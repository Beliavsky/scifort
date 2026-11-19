! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Independent implementation of the censored likelihood specified in
! https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.rv_continuous.fit.html
module scifort_censored
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_landau, only : landau_logpdf, landau_logcdf, landau_logsf
    use scifort_ksone, only : ksone_logsf
    use scifort_bessel_k, only : log_besselk_scaled
    use scifort_rel_breitwigner, only : rel_breitwigner_cdf, rel_breitwigner_sf
    use scifort_gausshyper, only : gausshyper_logkernel_integral
    use scifort_vonmises, only : vonmises_standard_logpdf
    use scifort_orthogonal_roots, only : roots_legendre
    use scifort_kinds, only : dp
    use scifort_constants, only : scifort_log_sqrt_two_pi, scifort_log_two, scifort_pi
    use scifort_empirical_distribution, only : censored_data
    use scifort_math, only : quiet_nan, negative_infinity, positive_infinity, valid_loc_scale, &
        expm1_safe, log1p_safe, log1pexp
    use scifort_error_fresnel, only : erfcx
    use scifort_normal, only : normal_logpdf, normal_logcdf, normal_logsf
    use scifort_exponential, only : exponential_logpdf, exponential_logcdf, exponential_logsf
    use scifort_log_gamma, only : log_gamma_one_plus, stirling_remainder, log_beta, &
        stirling_threshold, log_gamma_difference, log_gamma_ratio_scaled
    use scifort_incomplete_beta, only : incomplete_beta_xy, log_beta_kernel
    use scifort_gamma, only : gamma_logpdf
    use scifort_exponential_integrals, only : exp1
    use scifort_incomplete_gamma, only : log_gammainc, log_gammaincc, log_gamma_kernel
    use scifort_weibull, only : weibull_logpdf, weibull_logcdf, weibull_logsf
    implicit none
    private

    public :: censored_data
    public :: levy_stable_censored_loglikelihood, levy_stable_censored_nnlf
    public :: studentized_range_censored_loglikelihood, studentized_range_censored_nnlf
    public :: kstwo_censored_loglikelihood, kstwo_censored_nnlf
    public :: nct_censored_loglikelihood, nct_censored_nnlf
    public :: landau_censored_loglikelihood, landau_censored_nnlf
    public :: tukeylambda_censored_loglikelihood, tukeylambda_censored_nnlf
    public :: ksone_censored_loglikelihood, ksone_censored_nnlf
    public :: vonmises_censored_loglikelihood, vonmises_censored_nnlf
    public :: geninvgauss_censored_loglikelihood, geninvgauss_censored_nnlf
    public :: norminvgauss_censored_loglikelihood, norminvgauss_censored_nnlf
    public :: genhyperbolic_censored_loglikelihood, genhyperbolic_censored_nnlf
    public :: kstwobign_censored_loglikelihood, kstwobign_censored_nnlf
    public :: wrapcauchy_censored_loglikelihood, wrapcauchy_censored_nnlf
    public :: rel_breitwigner_censored_loglikelihood, rel_breitwigner_censored_nnlf
    public :: ncx2_censored_loglikelihood, ncx2_censored_nnlf
    public :: rice_censored_loglikelihood, rice_censored_nnlf
    public :: ncf_censored_loglikelihood, ncf_censored_nnlf
    public :: exponnorm_censored_loglikelihood, exponnorm_censored_nnlf
    public :: skewnorm_censored_loglikelihood, skewnorm_censored_nnlf
    public :: dpareto_lognorm_censored_loglikelihood, dpareto_lognorm_censored_nnlf
    public :: gausshyper_censored_loglikelihood, gausshyper_censored_nnlf
    public :: irwinhall_censored_loglikelihood, irwinhall_censored_nnlf
    public :: vonmises_line_censored_loglikelihood, vonmises_line_censored_nnlf
    public :: argus_censored_loglikelihood, argus_censored_nnlf
    public :: pearson3_censored_loglikelihood, pearson3_censored_nnlf
    public :: kappa4_censored_loglikelihood, kappa4_censored_nnlf
    public :: alpha_censored_loglikelihood, alpha_censored_nnlf
    public :: genextreme_censored_loglikelihood, genextreme_censored_nnlf
    public :: crystalball_censored_loglikelihood, crystalball_censored_nnlf
    public :: gengamma_censored_loglikelihood, gengamma_censored_nnlf
    public :: truncnorm_censored_loglikelihood, truncnorm_censored_nnlf
    public :: truncweibull_min_censored_loglikelihood, truncweibull_min_censored_nnlf
    public :: invgauss_censored_loglikelihood, invgauss_censored_nnlf
    public :: wald_censored_loglikelihood, wald_censored_nnlf
    public :: recipinvgauss_censored_loglikelihood, recipinvgauss_censored_nnlf
    public :: t_censored_loglikelihood, t_censored_nnlf
    public :: f_censored_loglikelihood, f_censored_nnlf
    public :: jf_skew_t_censored_loglikelihood, jf_skew_t_censored_nnlf
    public :: chi2_censored_loglikelihood, chi2_censored_nnlf
    public :: chi_censored_loglikelihood, chi_censored_nnlf
    public :: nakagami_censored_loglikelihood, nakagami_censored_nnlf
    public :: gennorm_censored_loglikelihood, gennorm_censored_nnlf
    public :: halfgennorm_censored_loglikelihood, halfgennorm_censored_nnlf
    public :: rdist_censored_loglikelihood, rdist_censored_nnlf
    public :: foldnorm_censored_loglikelihood, foldnorm_censored_nnlf
    public :: foldcauchy_censored_loglikelihood, foldcauchy_censored_nnlf
    public :: genpareto_censored_loglikelihood, genpareto_censored_nnlf
    public :: loguniform_censored_loglikelihood, loguniform_censored_nnlf
    public :: truncpareto_censored_loglikelihood, truncpareto_censored_nnlf
    public :: genhalflogistic_censored_loglikelihood, genhalflogistic_censored_nnlf
    public :: truncexpon_censored_loglikelihood, truncexpon_censored_nnlf
    public :: johnsonsb_censored_loglikelihood, johnsonsb_censored_nnlf
    public :: kappa3_censored_loglikelihood, kappa3_censored_nnlf
    public :: cosine_censored_loglikelihood, cosine_censored_nnlf
    public :: anglit_censored_loglikelihood, anglit_censored_nnlf
    public :: semicircular_censored_loglikelihood, semicircular_censored_nnlf
    public :: fatiguelife_censored_loglikelihood, fatiguelife_censored_nnlf
    public :: johnsonsu_censored_loglikelihood, johnsonsu_censored_nnlf
    public :: genexpon_censored_loglikelihood, genexpon_censored_nnlf
    public :: exponpow_censored_loglikelihood, exponpow_censored_nnlf
    public :: exponweib_censored_loglikelihood, exponweib_censored_nnlf
    public :: mielke_censored_loglikelihood, mielke_censored_nnlf
    public :: powernorm_censored_loglikelihood, powernorm_censored_nnlf
    public :: powerlognorm_censored_loglikelihood, powerlognorm_censored_nnlf
    public :: skewcauchy_censored_loglikelihood, skewcauchy_censored_nnlf
    public :: invgamma_censored_loglikelihood, invgamma_censored_nnlf
    public :: betaprime_censored_loglikelihood, betaprime_censored_nnlf
    public :: maxwell_censored_loglikelihood, maxwell_censored_nnlf
    public :: dgamma_censored_loglikelihood, dgamma_censored_nnlf
    public :: dweibull_censored_loglikelihood, dweibull_censored_nnlf
    public :: loggamma_censored_loglikelihood, loggamma_censored_nnlf
    public :: bradford_censored_loglikelihood, bradford_censored_nnlf
    public :: triang_censored_loglikelihood, triang_censored_nnlf
    public :: trapezoid_censored_loglikelihood, trapezoid_censored_nnlf
    public :: moyal_censored_loglikelihood, moyal_censored_nnlf
    public :: levy_censored_loglikelihood, levy_censored_nnlf
    public :: levy_l_censored_loglikelihood, levy_l_censored_nnlf
    public :: genlogistic_censored_loglikelihood, genlogistic_censored_nnlf
    public :: gompertz_censored_loglikelihood, gompertz_censored_nnlf
    public :: invweibull_censored_loglikelihood, invweibull_censored_nnlf
    public :: powerlaw_censored_loglikelihood, powerlaw_censored_nnlf
    public :: arcsine_censored_loglikelihood, arcsine_censored_nnlf
    public :: uniform_censored_loglikelihood, uniform_censored_nnlf
    public :: burr12_censored_loglikelihood, burr12_censored_nnlf
    public :: burr_censored_loglikelihood, burr_censored_nnlf
    public :: fisk_censored_loglikelihood, fisk_censored_nnlf
    public :: lomax_censored_loglikelihood, lomax_censored_nnlf
    public :: pareto_censored_loglikelihood, pareto_censored_nnlf
    public :: loglaplace_censored_loglikelihood, loglaplace_censored_nnlf

    integer, parameter, public :: censored_status_success = 0
    integer, parameter, public :: censored_status_invalid_input = 1

    abstract interface
        pure function censored_log_function(x, params) result(value)
            import dp
            implicit none
            real(dp), intent(in) :: x !! evaluation point, possibly infinite
            real(dp), intent(in) :: params(:) !! finite distribution parameters
            real(dp) :: value
        end function censored_log_function
        pure function censored_log_interval(a, b, params) result(value)
            import dp
            implicit none
            real(dp), intent(in) :: a !! lower endpoint, strictly less than b
            real(dp), intent(in) :: b !! upper endpoint, possibly infinite
            real(dp), intent(in) :: params(:) !! finite distribution parameters
            real(dp) :: value
        end function censored_log_interval
    end interface

    public :: make_censored_data, censored_data_valid, censored_size, censored_loglikelihood
    public :: censored_log_function, censored_log_interval
    public :: normal_censored_loglikelihood, normal_censored_nnlf
    public :: exponential_censored_loglikelihood, exponential_censored_nnlf
    public :: weibull_censored_loglikelihood, weibull_censored_nnlf
    public :: lognormal_censored_loglikelihood, lognormal_censored_nnlf
    public :: gibrat_censored_loglikelihood, gibrat_censored_nnlf
    public :: gamma_censored_loglikelihood, gamma_censored_nnlf
    public :: erlang_censored_loglikelihood, erlang_censored_nnlf
    public :: beta_censored_loglikelihood, beta_censored_nnlf
    public :: rayleigh_censored_loglikelihood, rayleigh_censored_nnlf
    public :: weibull_max_censored_loglikelihood, weibull_max_censored_nnlf
    public :: halfcauchy_censored_loglikelihood, halfcauchy_censored_nnlf
    public :: halfnorm_censored_loglikelihood, halfnorm_censored_nnlf
    public :: hypsecant_censored_loglikelihood, hypsecant_censored_nnlf
    public :: laplace_asymmetric_censored_loglikelihood, laplace_asymmetric_censored_nnlf
    public :: halflogistic_censored_loglikelihood, halflogistic_censored_nnlf
    public :: cauchy_censored_loglikelihood, cauchy_censored_nnlf
    public :: laplace_censored_loglikelihood, laplace_censored_nnlf
    public :: logistic_censored_loglikelihood, logistic_censored_nnlf
    public :: gumbel_r_censored_loglikelihood, gumbel_r_censored_nnlf
    public :: gumbel_l_censored_loglikelihood, gumbel_l_censored_nnlf

contains

    ! Independent defining integrals; see CODE_PROVENANCE.md.
    pure function ls_logs(z, p, mode) result(value)
        implicit none
        real(dp), intent(in) :: z !! central standardized coordinate, abs(z)<=6
        real(dp), intent(in) :: p(:) !! alpha, loc, scale, beta and GL16 rule
        integer, intent(in) :: mode !! zero density, one CDF, two SF
        real(dp) :: value, upper, h, v, t, u, phase, term, total, point
        integer :: j, i, panels
        value = quiet_nan(0.0_dp)
        if (.not. ieee_is_finite(z)) return
        if (abs(z) > 6.0_dp) return
        point = z
        if (mode == 2) point = -z
        upper = 48.0_dp**(1.0_dp/(4.0_dp*p(1)))
        panels = ceiling(upper/0.025_dp)
        h = upper/real(panels,dp)
        total = 0.0_dp
        do j = 1, panels
            do i = 1, 16
                v = (real(j-1,dp)+0.5_dp*(1.0_dp+p(4+i)))*h
                t = v**4
                u = v**(4.0_dp*p(1))
                phase = p(4)*tan(0.5_dp*scifort_pi*p(1))*u
                if (mode == 2) phase = -phase
                phase = point*t-phase
                if (mode == 0) then
                    term = 4.0_dp*v**3*exp(-u)*cos(phase)
                else
                    term = 4.0_dp*exp(-u)*sin(phase)/v
                end if
                total = total+p(20+i)*term
            end do
        end do
        total = total*0.5_dp*h/scifort_pi
        if (mode /= 0) total = 0.5_dp+total
        if (total <= 0.0_dp) return
        value = log(total)
        if (mode == 0) value = value-log(p(3))
    end function ls_logs

    ! Rotate the characteristic-function contour into the lower half-plane.
    ! Subtract its zero-scale term analytically before quadrature; see provenance.
    pure function ls_tail(lz, beta, p, density) result(value)
        implicit none
        real(dp), intent(in) :: lz !! logarithm of positive standardized distance, >=log(6)
        real(dp), intent(in) :: beta !! reflected S1 skewness in [-0.5,0.5]
        real(dp), intent(in) :: p(:) !! alpha, affine parameters and GL16 rule
        logical, intent(in) :: density !! density rather than upper-tail probability
        real(dp) :: value, theta, q, v, y, total, term
        complex(dp) :: rotation, coefficient, base, w, ratio, power, integrand
        integer :: j, i, n
        theta = scifort_pi/(8.0_dp*p(1))
        rotation = cmplx(cos(theta),-sin(theta),dp)
        coefficient = cmplx(1.0_dp,-beta*tan(0.5_dp*scifort_pi*p(1)),dp) &
            *cmplx(cos(p(1)*theta),-sin(p(1)*theta),dp)
        q = exp(-p(1)*lz)
        total = 0.0_dp
        ! y=v**4, 0<=v<=4. The damping is at least exp(-sin(pi/14)*y).
        do j = 1, 160
            do i = 1, 16
                v = (real(j-1,dp)+0.5_dp*(1.0_dp+p(4+i)))*0.025_dp
                y = v**4
                base = -coefficient*y**p(1)
                w = base*q
                if (abs(w) < 0.5_dp) then
                    ! Entire function (exp(w)-1)/w, including q=0 after underflow.
                    ratio = cmplx(1.0_dp,0.0_dp,dp)
                    power = ratio
                    do n = 1, 20
                        power = power*w/real(n+1,dp)
                        ratio = ratio+power
                    end do
                else
                    ratio = (exp(w)-1.0_dp)/w
                end if
                integrand = base*ratio*exp(cmplx(0.0_dp,-y,dp)*rotation)
                if (density) then
                    term = real(rotation*integrand,dp)*4.0_dp*v**3
                else
                    term = aimag(integrand)*4.0_dp/v
                end if
                total = total+p(20+i)*term
            end do
        end do
        total = total*0.0125_dp/scifort_pi
        value = quiet_nan(0.0_dp)
        if (total <= 0.0_dp) return
        value = log(total)-p(1)*lz
        if (density) value = value-lz-log(p(3))
    end function ls_tail

    pure function ls_value(x, p, mode) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite observation or censoring endpoint
        real(dp), intent(in) :: p(:) !! alpha, affine parameters, beta and GL16 rule
        integer, intent(in) :: mode !! zero density, one CDF, two SF
        real(dp) :: value, z, lz, beta
        z = location_scale_z(x,p(2:3))
        if (abs(z) <= 6.0_dp) then
            value = ls_logs(z,p,mode)
            return
        end if
        lz = log_distance(max(x,p(2)),min(x,p(2)))-log(p(3))
        beta = p(4)
        if (x < p(2)) beta = -beta
        value = ls_tail(lz,beta,p,mode == 0)
        if ((mode == 1 .and. x > p(2)) .or. (mode == 2 .and. x < p(2))) &
            value = ig_complement(value)
    end function ls_value

    pure function ls_pdf(x, p) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: p(:) !! stable shape, affine parameters and rule
        real(dp) :: value
        value = ls_value(x,p,0)
    end function ls_pdf

    pure function ls_cdf(x, p) result(value)
        implicit none
        real(dp), intent(in) :: x !! upper censoring endpoint
        real(dp), intent(in) :: p(:) !! stable shape, affine parameters and rule
        real(dp) :: value
        value = ls_value(x,p,1)
    end function ls_cdf

    pure function ls_sf(x, p) result(value)
        implicit none
        real(dp), intent(in) :: x !! lower censoring endpoint
        real(dp), intent(in) :: p(:) !! stable shape, affine parameters and rule
        real(dp) :: value
        value = ls_value(x,p,2)
    end function ls_sf

    pure function integrated_censored_interval(a, b, p, pdf, cdf, sf, bounded) result(value)
        implicit none
        real(dp), intent(in) :: a, b !! ordered finite interval endpoints
        real(dp), intent(in) :: p(:) !! shapes, affine parameters and GL16 rule
        procedure(censored_log_function) :: pdf, cdf, sf !! log density and tails
        logical, intent(in) :: bounded !! two-sided stable support rather than positive range
        real(dp) :: value, left, la, lb, qa, qb, width, u, x, term, za, zb
        integer :: j, i, panels
        value = quiet_nan(0.0_dp)
        left = a
        if (.not. bounded) left = max(a,p(2))
        if (b <= left) then
            value = negative_infinity(0.0_dp)
            return
        end if
        if (.not. bounded .and. left == p(2)) then
            value = cdf(b,p)
            return
        end if
        la = cdf(left,p)
        lb = cdf(b,p)
        qa = sf(left,p)
        qb = sf(b,p)
        if (any(ieee_is_nan([la,lb,qa,qb]))) return
        if (lb < -scifort_log_two) then
            if (la < lb-0.01_dp) then
                value = log_difference(lb,la)
                return
            end if
        else if (qb < qa-0.01_dp) then
            value = log_difference(qa,qb)
            return
        end if
        width = log_distance(b,left)
        panels = 1
        if (bounded) then
            ! Tail intervals reaching this fallback are relatively narrow.
            ! Only the central segment needs absolute-width subdivision.
            za = max(-6.0_dp,min(6.0_dp,location_scale_z(left,p(2:3))))
            zb = max(-6.0_dp,min(6.0_dp,location_scale_z(b,p(2:3))))
            panels = max(1,ceiling((zb-za)/0.25_dp))
        end if
        value = negative_infinity(0.0_dp)
        do j = 1, panels
            do i = 1, 16
                u = (real(j-1,dp)+0.5_dp*(1.0_dp+p(4+i)))/real(panels,dp)
                x = (1.0_dp-u)*left+u*b
                term = pdf(x,p)+log(p(20+i))
                if (ieee_is_nan(term)) then
                    value = term
                    return
                end if
                value = bounded_logadd(value,term)
            end do
        end do
        value = min(0.0_dp,value+width-scifort_log_two-log(real(panels,dp)))
    end function integrated_censored_interval

    pure function ls_interval(a, b, p) result(value)
        implicit none
        real(dp), intent(in) :: a, b !! ordered finite interval endpoints
        real(dp), intent(in) :: p(:) !! stable shape, affine parameters and rule
        real(dp) :: value
        value = integrated_censored_interval(a,b,p,ls_pdf,ls_cdf,ls_sf,.true.)
    end function ls_interval

    pure function levy_stable_censored_loglikelihood(data, alpha, beta, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: alpha !! stability shape; see documented supported subsets
        real(dp), intent(in) :: beta !! S1 skewness in [-1,1], subject to subset limits
        real(dp), intent(in), optional :: loc !! finite S1 location, default zero
        real(dp), intent(in), optional :: scale !! normal positive scale, default one
        real(dp) :: value, p(36)
        real(dp), allocatable :: nodes(:), weights(:)
        value = quiet_nan(0.0_dp)
        p = 0.0_dp
        p(1:4) = [alpha,0.0_dp,1.0_dp,beta]
        if (present(loc)) p(2) = loc
        if (present(scale)) p(3) = scale
        if (.not. all(ieee_is_finite(p))) return
        if (p(3) < tiny(1.0_dp) .or. p(3) > huge(1.0_dp)/100.0_dp) return
        if (abs(p(2))/1.0e4_dp > p(3)) return
        if (abs(beta) > 1.0_dp) return
        if (alpha == 2.0_dp) then
            value = normal_censored_loglikelihood(data,p(2),sqrt(2.0_dp)*p(3))
        else if (alpha == 1.0_dp .and. beta == 0.0_dp) then
            value = cauchy_censored_loglikelihood(data,p(2),p(3))
        else if (alpha == 0.5_dp .and. beta == 1.0_dp) then
            value = levy_censored_loglikelihood(data,p(2),p(3))
        else if (alpha == 0.5_dp .and. beta == -1.0_dp) then
            value = levy_l_censored_loglikelihood(data,p(2),p(3))
        else
            if (alpha < 1.25_dp .or. alpha > 1.75_dp .or. abs(beta) > 0.5_dp) return
            call roots_legendre(16,nodes,weights)
            p(5:20) = nodes
            p(21:36) = weights
            value = censored_loglikelihood(data,p,ls_pdf,ls_cdf,ls_sf,ls_interval)
        end if
    end function levy_stable_censored_loglikelihood

    pure function levy_stable_censored_nnlf(data, alpha, beta, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: alpha !! stability shape in supported subsets
        real(dp), intent(in) :: beta !! S1 skewness in supported subsets
        real(dp), intent(in), optional :: loc !! finite S1 location, default zero
        real(dp), intent(in), optional :: scale !! normal positive scale, default one
        real(dp) :: value
        value = -levy_stable_censored_loglikelihood(data,alpha,beta,loc,scale)
    end function levy_stable_censored_nnlf

    pure function sr_range(lr, k, mode, p) result(value)
        implicit none
        real(dp), intent(in) :: lr !! log positive normal range
        real(dp), intent(in) :: k !! number of means in [2,10]
        integer, intent(in) :: mode !! zero PDF, one CDF, two SF
        real(dp), intent(in) :: p(:) !! shapes, affine parameters and GL16 rule
        real(dp) :: value, r, center, z, ld, lq, ratio, hazard, term
        integer :: j, i
        r = exp(lr)
        if (k == 2.0_dp) then
            z = r/sqrt(2.0_dp)
            select case (mode)
            case (0)
                value = scifort_log_two/2+normal_logpdf(z)
            case (1)
                value = normal_small_interval(0.0_dp,lr-scifort_log_two/2)
                if (ieee_is_nan(value)) value = normal_interval(0.0_dp,z,[0.0_dp,1.0_dp])
                value = value+scifort_log_two
            case (2)
                value = scifort_log_two+normal_logsf(z)
            end select
            return
        end if
        if (mode == 1 .and. r > 40.0_dp) then
            value = 0.0_dp
            return
        end if
        center = 0.0_dp
        if (mode /= 1) center = -0.5_dp*r
        value = negative_infinity(0.0_dp)
        do j = 1, 24
            do i = 1, 16
                z = center-12.0_dp+real(j-1,dp)+0.5_dp*(1.0_dp+p(4+i))
                if (mode /= 2) then
                    ld = normal_small_interval(z,lr)
                    if (ieee_is_nan(ld)) ld = normal_interval(z,z+r,[0.0_dp,1.0_dp])
                end if
                select case (mode)
                case (0)
                    term = log(k*(k-1.0_dp))+normal_logpdf(z)+normal_logpdf(z+r)+(k-2.0_dp)*ld
                case (1)
                    term = log(k)+normal_logpdf(z)+(k-1.0_dp)*ld
                case (2)
                    lq = normal_logsf(z)
                    ratio = min(0.0_dp,normal_logsf(z+r)-lq)
                    if (ratio < -36.0_dp) then
                        hazard = log(k-1.0_dp)+ratio
                    else if (ratio == 0.0_dp) then
                        hazard = positive_infinity(0.0_dp)
                    else
                        hazard = log(k-1.0_dp)+log(-log(-expm1_safe(ratio)))
                    end if
                    term = log(k)+normal_logpdf(z)+(k-1.0_dp)*lq+log_exponential_mass(hazard)
                end select
                value = bounded_logadd(value,term+log(p(20+i)))
            end do
        end do
        value = value-scifort_log_two
    end function sr_range

    pure function sr_integral(lq, p, mode) result(value)
        implicit none
        real(dp), intent(in) :: lq !! log positive standardized variate
        real(dp), intent(in) :: p(:) !! k, location, scale, encoded df and GL16 rule
        integer, intent(in) :: mode !! zero PDF, one CDF, two SF
        real(dp) :: value, lo, hi, h, v, y, lr, term, df, constant
        integer :: j, i, panels
        if (p(4) == 0.0_dp) then
            if (lq > log(40.0_dp)) then
                ! Reject coordinates beyond the controlled infinite-df range.
                value = quiet_nan(0.0_dp)
                return
            end if
            value = sr_range(lq,p(1),mode,p)
            return
        end if
        df = p(4)
        constant = scifort_log_two+0.5_dp*df*log(0.5_dp*df)-log_gamma(0.5_dp*df)
        lo = -100.0_dp/df
        hi = 0.5_dp*log((df+200.0_dp)/df)
        if (mode /= 1) then
            lo = min(lq,0.0_dp)+lo
            hi = min(lq+hi,log(40.0_dp))
        end if
        panels = max(1,ceiling((hi-lo)/0.25_dp))
        h = (hi-lo)/real(panels,dp)
        value = negative_infinity(0.0_dp)
        do j = 1, panels
            do i = 1, 16
                v = lo+(real(j-1,dp)+0.5_dp*(1.0_dp+p(4+i)))*h
                y = v
                lr = lq+v
                if (mode /= 1) then
                    y = v-lq
                    lr = v
                end if
                if (y > 0.5_dp*log((df+200.0_dp)/df)) cycle
                if (lr > log(40.0_dp)) then
                    term = 0.0_dp
                else
                    term = sr_range(lr,p(1),mode,p)
                end if
                term = term+constant+df*y-0.5_dp*df*exp(2.0_dp*y)
                if (mode == 0) term = term+y
                value = bounded_logadd(value,term+log(p(20+i)))
            end do
        end do
        value = value+log(0.5_dp*h)
    end function sr_integral

    pure function sr_value(x, p, mode) result(value)
        implicit none
        real(dp), intent(in) :: x !! exact observation or censoring endpoint
        real(dp), intent(in) :: p(:) !! k, location, scale, encoded df and GL16 rule
        integer, intent(in) :: mode !! zero PDF, one CDF, two SF
        real(dp) :: value, z, first
        integer :: tail
        value = negative_infinity(0.0_dp)
        if (x <= p(2)) then
            if (mode == 2) value = 0.0_dp
            if (x == p(2) .and. mode == 0 .and. p(1) == 2.0_dp) then
                value = -0.5_dp*log(scifort_pi)-log(p(3))
                if (p(4) /= 0.0_dp) value = value+0.5_dp*log(2.0_dp/p(4))+ &
                    log_gamma(0.5_dp*(p(4)+1.0_dp))-log_gamma(0.5_dp*p(4))
            end if
            return
        end if
        z = log_distance(x,p(2))-log(p(3))
        if (mode == 0) then
            value = sr_integral(z,p,0)-log(p(3))
            return
        end if
        tail = 1
        if (z > scifort_log_two) tail = 2
        first = sr_integral(z,p,tail)
        value = first
        if (ieee_is_nan(value)) return
        if (first > -scifort_log_two) then
            tail = 3-tail
            value = sr_integral(z,p,tail)
            if (ieee_is_nan(value)) return
        end if
        value = min(0.0_dp,value)
        if (mode /= tail) value = ig_complement(value)
    end function sr_value

    pure function sr_pdf(x, p) result(v)
        implicit none
        real(dp), intent(in) :: x !! exact finite observation
        real(dp), intent(in) :: p(:) !! packed studentized-range parameters
        real(dp) :: v
        v = sr_value(x,p,0)
    end function sr_pdf

    pure function sr_cdf(x, p) result(v)
        implicit none
        real(dp), intent(in) :: x !! upper censoring endpoint
        real(dp), intent(in) :: p(:) !! packed studentized-range parameters
        real(dp) :: v
        v = sr_value(x,p,1)
    end function sr_cdf

    pure function sr_sf(x, p) result(v)
        implicit none
        real(dp), intent(in) :: x !! lower censoring endpoint
        real(dp), intent(in) :: p(:) !! packed studentized-range parameters
        real(dp) :: v
        v = sr_value(x,p,2)
    end function sr_sf

    pure function sr_interval(a, b, p) result(v)
        implicit none
        real(dp), intent(in) :: a, b !! ordered finite interval endpoints
        real(dp), intent(in) :: p(:) !! packed studentized-range parameters
        real(dp) :: v
        v = integrated_censored_interval(a,b,p,sr_pdf,sr_cdf,sr_sf,.false.)
    end function sr_interval

    pure function studentized_range_censored_loglikelihood(data, k, df, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: k !! number of means in [2,10]
        real(dp), intent(in) :: df !! degrees of freedom: [1,100] or >=100000 (asymptotic)
        real(dp), intent(in), optional :: loc !! lower support location, default zero
        real(dp), intent(in), optional :: scale !! normal positive scale, default one
        real(dp) :: value, p(36)
        real(dp), allocatable :: nodes(:), weights(:)
        value = quiet_nan(0.0_dp)
        if (ieee_is_nan(df)) return
        if (df < 1.0_dp) return
        if (df > 100.0_dp .and. df < 1.0e5_dp) return
        p = 0.0_dp
        p(1:4) = [k,0.0_dp,1.0_dp,0.0_dp]
        if (df <= 100.0_dp) p(4) = df
        if (present(loc)) p(2) = loc
        if (present(scale)) p(3) = scale
        if (.not. all(ieee_is_finite(p))) return
        if (k < 2.0_dp .or. k > 10.0_dp) return
        if (p(3) < tiny(1.0_dp) .or. p(3) > huge(1.0_dp)/100.0_dp) return
        if (abs(p(2))/1.0e4_dp > p(3)) return
        call roots_legendre(16,nodes,weights)
        p(5:20) = nodes
        p(21:36) = weights
        value = censored_loglikelihood(data,p,sr_pdf,sr_cdf,sr_sf,sr_interval)
    end function studentized_range_censored_loglikelihood

    pure function studentized_range_censored_nnlf(data, k, df, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: k !! number of means in [2,10]
        real(dp), intent(in) :: df !! degrees of freedom: [1,100] or >=100000 (asymptotic)
        real(dp), intent(in), optional :: loc !! lower support location, default zero
        real(dp), intent(in), optional :: scale !! normal positive scale, default one
        real(dp) :: value
        value = -studentized_range_censored_loglikelihood(data,k,df,loc,scale)
    end function studentized_range_censored_nnlf


    ! Uniform-order-statistic rectangle volumes and positive boundary fluxes.
    pure subroutine kw_volume(n, d, cdf, sf, pdf)
        implicit none
        integer, intent(in) :: n !! sample size in [2,32]
        real(dp), intent(in) :: d !! standardized statistic strictly above 1/n
        real(dp), intent(out) :: cdf, sf, pdf !! lower tail, upper tail and density
        real(dp) :: nodes(0:65), lower(32), upper(32), f(0:32,0:65), g(0:32,0:65)
        real(dp) :: v(0:32), w(0:32), powers(0:32), fac, y, delta, hold, rest
        integer :: m, i, j, k, r, countlo, counthi, at
        nodes(0:1) = [0.0_dp,1.0_dp]
        m = 1
        do i = 1, n
            lower(i) = real(i,dp)/real(n,dp)-d
            upper(i) = real(i-1,dp)/real(n,dp)+d
            if (lower(i) > 0.0_dp .and. lower(i) < 1.0_dp) then
                m = m+1
                nodes(m) = lower(i)
            end if
            if (upper(i) > 0.0_dp .and. upper(i) < 1.0_dp) then
                m = m+1
                nodes(m) = upper(i)
            end if
        end do
        do i = 1, m
            hold = nodes(i)
            j = i-1
            do while (j >= 0)
                if (nodes(j) <= hold) exit
                nodes(j+1) = nodes(j)
                j = j-1
            end do
            nodes(j+1) = hold
        end do
        at = 0
        do i = 1, m
            if (nodes(i) == nodes(at)) cycle
            at = at+1
            nodes(at) = nodes(i)
        end do
        m = at
        fac = exp(log_gamma(real(n+1,dp)))
        f = 0.0_dp
        g = 0.0_dp
        v = 0.0_dp
        v(0) = 1.0_dp
        sf = 0.0_dp
        do j = 1, m
            y = nodes(j)
            delta = y-nodes(j-1)
            powers(0) = 1.0_dp
            do k = 1, n
                powers(k) = powers(k-1)*delta/real(k,dp)
            end do
            w = 0.0_dp
            do r = 0, n
                do k = 0, r
                    w(r) = w(r)+v(r-k)*powers(k)
                end do
            end do
            f(:,j) = w
            countlo = 0
            counthi = n
            do i = 1, n
                if (upper(i) == y) countlo = max(countlo,i)
                if (lower(i) == y) counthi = min(counthi,i-1)
            end do
            do r = 0, n
                if (r >= countlo .and. r <= counthi) cycle
                rest = 1.0_dp
                do k = 1, n-r
                    rest = rest*(1.0_dp-y)/real(k,dp)
                end do
                sf = sf+fac*w(r)*rest
                w(r) = 0.0_dp
            end do
            v = w
        end do
        cdf = fac*v(n)
        v = 0.0_dp
        v(0) = 1.0_dp
        do j = m-1, 0, -1
            y = nodes(j)
            delta = nodes(j+1)-y
            powers(0) = 1.0_dp
            do k = 1, n
                powers(k) = powers(k-1)*delta/real(k,dp)
            end do
            w = 0.0_dp
            do r = 0, n
                do k = 0, r
                    w(r) = w(r)+v(r-k)*powers(k)
                end do
            end do
            g(:,j) = w
            countlo = 0
            counthi = n
            do i = 1, n
                if (lower(i) == y) countlo = max(countlo,n-i+1)
                if (upper(i) == y) counthi = min(counthi,n-i)
            end do
            do r = 0, n
                if (r < countlo .or. r > counthi) w(r) = 0.0_dp
            end do
            v = w
        end do
        pdf = 0.0_dp
        do i = 1, n
            do j = 1, m-1
                y = nodes(j)
                if (lower(i) == y) pdf = pdf+fac*f(i-1,j)*g(n-i,j)
                if (upper(i) == y) pdf = pdf+fac*f(i-1,j)*g(n-i,j)
            end do
        end do
    end subroutine kw_volume

    pure function kw_density(z, p) result(value)
        implicit none
        real(dp), intent(in) :: z !! interior standard coordinate; right density at 1/n
        real(dp), intent(in) :: p(:) !! n, location, scale, lower, upper, knot and rule
        real(dp) :: value, cdf, sf, pdf, d
        integer :: n
        n = int(p(1))
        if (n == 1) then
            value = scifort_log_two
        else if (z >= 1.0_dp-1.0_dp/p(1)) then
            value = log(2.0_dp*p(1))+(p(1)-1.0_dp)*log1p_safe(-z)
        else
            d = max(z,nearest(1.0_dp/p(1),1.0_dp))
            call kw_volume(n,d,cdf,sf,pdf)
            value = log(pdf)
        end if
    end function kw_density

    pure subroutine kw_logs(x, p, ld, lp, lq)
        implicit none
        real(dp), intent(in) :: x !! finite physical observation
        real(dp), intent(in) :: p(:) !! n, location, scale, lower, upper, knot and rule
        real(dp), intent(out) :: ld, lp, lq !! log density, CDF and SF
        real(dp) :: z, lw, lu, cdf, sf, pdf, nf
        ld = negative_infinity(0.0_dp)
        lp = ld
        lq = 0.0_dp
        if (x <= p(4)) return
        lp = 0.0_dp
        lq = negative_infinity(0.0_dp)
        if (x >= p(5)) return
        lu = log_distance(p(5),x)-log(p(3))
        nf = log_gamma(p(1)+1.0_dp)
        if (p(1) == 1.0_dp) then
            lp = scifort_log_two+log_distance(x,p(4))-log(p(3))
            lq = scifort_log_two+lu
            ld = scifort_log_two-log(p(3))
        else if (x < p(6)) then
            lw = scifort_log_two+log_distance(x,p(4))-log(p(3))
            lp = nf+p(1)*lw
            lq = ig_complement(lp)
            ld = log(2.0_dp*p(1))+nf+(p(1)-1.0_dp)*lw-log(p(3))
        else if (x == p(6)) then
            lp = nf-p(1)*log(p(1))
            lq = ig_complement(lp)
            ld = log(2.0_dp*p(1))+nf-(p(1)-1.0_dp)*log(p(1))-log(p(3))
        else if (lu <= -log(p(1))) then
            lq = scifort_log_two+p(1)*lu
            lp = ig_complement(lq)
            ld = log(2.0_dp*p(1))+(p(1)-1.0_dp)*lu-log(p(3))
        else
            z = max(location_scale_z(x,p(2:3)),nearest(1.0_dp/p(1),1.0_dp))
            call kw_volume(int(p(1)),z,cdf,sf,pdf)
            lp = log(cdf)
            lq = log(sf)
            ld = log(pdf)-log(p(3))
        end if
        if (lp <= lq) then
            lq = ig_complement(lp)
        else
            lp = ig_complement(lq)
        end if
    end subroutine kw_logs

    pure function kw_pdf(x, p) result(v)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: p(:) !! packed parameters
        real(dp) :: v, lp, lq
        call kw_logs(x,p,v,lp,lq)
    end function kw_pdf

    pure function kw_cdf(x, p) result(v)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: p(:) !! packed parameters
        real(dp) :: v, ld, lq
        call kw_logs(x,p,ld,v,lq)
    end function kw_cdf

    pure function kw_sf(x, p) result(v)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: p(:) !! packed parameters
        real(dp) :: v, ld, lp
        call kw_logs(x,p,ld,lp,v)
    end function kw_sf

    pure function kw_interval(a, b, p) result(v)
        implicit none
        real(dp), intent(in) :: a, b !! finite ordered physical endpoints
        real(dp), intent(in) :: p(:) !! packed parameters
        real(dp) :: v, left, right, la, lb, qa, qb, unused, lw, d, next, cut, z, base, term
        integer :: i, k
        left = max(a,p(4))
        right = min(b,p(5))
        v = negative_infinity(0.0_dp)
        if (right <= left) return
        if (left == p(4)) then
            v = kw_cdf(right,p)
            return
        else if (right == p(5)) then
            v = kw_sf(left,p)
            return
        end if
        lw = log_distance(right,left)
        if (p(1) == 1.0_dp) then
            v = scifort_log_two+lw-log(p(3))
            return
        else if (right <= p(6)) then
            d = log(p(1))+burr_logsoftplus(lw-log_distance(left,p(4)))
            v = kw_cdf(right,p)+log_exponential_mass(d)
            return
        else if (log_distance(p(5),left)-log(p(3)) <= -log(p(1))) then
            d = log(p(1))+burr_logsoftplus(lw-log_distance(p(5),right))
            v = kw_sf(left,p)+log_exponential_mass(d)
            return
        end if
        call kw_logs(left,p,unused,la,qa)
        call kw_logs(right,p,unused,lb,qb)
        if (lb < -scifort_log_two) then
            if (la < lb-0.01_dp) then
                v = log_difference(lb,la)
                return
            end if
        else if (qb < qa-0.01_dp) then
            v = log_difference(qa,qb)
            return
        end if
        ! Each density piece has degree at most n-1<=31; GL16 integrates it.
        do k = 2, 2*int(p(1))
            cut = p(2)+p(3)*(real(k,dp)/(2.0_dp*p(1)))
            if (k == 2) cut = p(6)
            if (k == 2*int(p(1))) cut = p(5)
            next = min(right,cut)
            if (next <= left) cycle
            if (next <= p(6)) then
                d = log(p(1))+burr_logsoftplus(log_distance(next,left)-log_distance(left,p(4)))
                term = kw_cdf(next,p)+log_exponential_mass(d)
            else
                lw = log_distance(next,left)-log(p(3))
                base = location_scale_z(left,p(2:3))
                if (left == p(6)) base = 1.0_dp/p(1)
                term = negative_infinity(0.0_dp)
                do i = 1, 16
                    z = base+0.5_dp*(1.0_dp+p(6+i))*exp(lw)
                    d = kw_density(z,p)+log(p(22+i))
                    term = bounded_logadd(term,d)
                end do
                term = term+lw-scifort_log_two
            end if
            v = bounded_logadd(v,term)
            left = next
            if (left == right) exit
        end do
        v = min(0.0_dp,v)
    end function kw_interval

    pure function kstwo_censored_loglikelihood(data, n, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: n !! integer sample size in [1,32]
        real(dp), intent(in), optional :: loc !! finite location with abs(loc)/scale<=100
        real(dp), intent(in), optional :: scale !! normal positive scale, default one
        real(dp) :: value, p(38)
        real(dp), allocatable :: nodes(:), weights(:)
        value = quiet_nan(0.0_dp)
        p = 0.0_dp
        p(1:3) = [n,0.0_dp,1.0_dp]
        if (present(loc)) p(2) = loc
        if (present(scale)) p(3) = scale
        if (.not. all(ieee_is_finite(p))) return
        if (n < 1.0_dp .or. n > 32.0_dp) return
        if (n /= real(nint(n),dp)) return
        if (p(3) < tiny(1.0_dp) .or. p(3) > huge(1.0_dp)/2.0_dp) return
        if (abs(p(2))/100.0_dp > p(3)) return
        if (p(2) > huge(1.0_dp)-p(3)) return
        p(4) = p(2)+p(3)/(2.0_dp*n)
        p(5) = p(2)+p(3)
        p(6) = p(2)+p(3)/n
        if (p(4) >= p(6) .or. p(6) > p(5)) return
        call roots_legendre(16,nodes,weights)
        p(7:22) = nodes
        p(23:38) = weights
        value = censored_loglikelihood(data,p,kw_pdf,kw_cdf,kw_sf,kw_interval)
    end function kstwo_censored_loglikelihood

    pure function kstwo_censored_nnlf(data, n, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: n !! integer sample size in [1,32]
        real(dp), intent(in), optional :: loc !! finite location with abs(loc)/scale<=100
        real(dp), intent(in), optional :: scale !! normal positive scale, default one
        real(dp) :: value
        value = -kstwo_censored_loglikelihood(data,n,loc,scale)
    end function kstwo_censored_nnlf


    ! Positive normal/chi mixture in log(sqrt(V/df)); see CODE_PROVENANCE.md.
    pure function nt_integral(lz, mode, p) result(value)
        implicit none
        real(dp), intent(in) :: lz !! log positive standardized magnitude after reflection
        integer, intent(in) :: mode !! zero density, one upper tail, two lower tail
        real(dp), intent(in) :: p(:) !! family, location, scale, df, nc, normalizer and rule
        real(dp) :: value, lo, hi, h, v, term, arg, lr, df, nc, normal
        integer :: n, i, j
        df = p(4)
        nc = p(5)
        lo = min(0.0_dp,-lz)-(100.0_dp+nc*nc)/df
        hi = 0.5_dp*log((df+200.0_dp+nc*nc)/df)
        if (mode /= 2) hi = min(hi,log(abs(nc)+40.0_dp)-lz)
        n = max(1,ceiling((hi-lo)/0.25_dp))
        h = (hi-lo)/real(n,dp)
        value = negative_infinity(0.0_dp)
        do j = 1, n
            do i = 1, 16
                v = lo+(real(j-1,dp)+0.5_dp*(1.0_dp+p(6+i)))*h
                lr = lz+v
                if (lr > log(abs(nc)+40.0_dp)) then
                    if (mode /= 2) cycle
                    normal = 0.0_dp
                else
                    arg = exp(lr)-nc
                    select case (mode)
                    case (0)
                        normal = -0.5_dp*arg*arg-scifort_log_sqrt_two_pi
                    case (1)
                        normal = normal_logsf(arg)
                    case default
                        normal = normal_logcdf(arg)
                    end select
                end if
                term = df*v-0.5_dp*df*exp(2.0_dp*v)+normal+log(p(22+i))
                if (mode == 0) term = term+v
                value = bounded_logadd(value,term)
            end do
        end do
        value = value+p(6)+log(0.5_dp*h)
        if (mode /= 0) value = min(0.0_dp,value)
    end function nt_integral

    pure subroutine sp_coordinate(x, p, z, lz)
        implicit none
        real(dp), intent(in) :: x !! finite physical observation
        real(dp), intent(in) :: p(:) !! packed specialist distribution parameters
        real(dp), intent(out) :: z !! standardized value, possibly infinite
        real(dp), intent(out) :: lz !! log absolute standardized value
        lz = negative_infinity(0.0_dp)
        if (x > p(2)) lz = log_distance(x,p(2))-log(p(3))
        if (x < p(2)) lz = log_distance(p(2),x)-log(p(3))
        z = location_scale_z(x,p(2:3))
    end subroutine sp_coordinate

    pure subroutine sp_logs(x, p, ld, lp, lq)
        implicit none
        real(dp), intent(in) :: x !! finite physical observation
        real(dp), intent(in) :: p(:) !! family, location, scale, shapes and quadrature rule
        real(dp), intent(out) :: ld !! log density in physical units
        real(dp), intent(out) :: lp !! log CDF
        real(dp), intent(out) :: lq !! log SF
        real(dp) :: z, lz, pp(38), small, large
        call sp_coordinate(x,p,z,lz)
        if (p(1) == 1.0_dp) then
            pp = p
            if (x < p(2)) pp(5) = -p(5)
            small = nt_integral(lz,1,pp)
            if (small <= -scifort_log_two) then
                large = ig_complement(small)
            else
                large = nt_integral(lz,2,pp)
                small = ig_complement(large)
            end if
            if (x < p(2)) then
                lp = small
                lq = large
            else
                lp = large
                lq = small
            end if
            ld = nt_integral(lz,0,pp)-log(p(3))
        else
            if (z <= -1.0_dp) then
                call ld_left(z,p,ld,lp)
                ld = ld-log(p(3))
                lq = ig_complement(lp)
            else if (lz >= log(1.0e100_dp) .and. x > p(2)) then
                lq = log(2.0_dp/scifort_pi)-lz
                lp = ig_complement(lq)
                ld = lq-lz-log(p(3))
            else
                ld = landau_logpdf(z)-log(p(3))
                lp = landau_logcdf(z)
                lq = landau_logsf(z)
                if (lp <= -scifort_log_two) then
                    lq = ig_complement(lp)
                else
                    lp = ig_complement(lq)
                end if
            end if
        end if
    end subroutine sp_logs

    ! Bromwich inversion on the real saddle of E[exp(-s X)]=exp((2/pi)*s*log(s)).
    pure subroutine ld_left(z, p, ld, lp)
        implicit none
        real(dp), intent(in) :: z !! standardized observation at most -1
        real(dp), intent(in) :: p(:) !! packed Gauss-Legendre nodes and weights
        real(dp), intent(out) :: ld !! standard log density
        real(dp), intent(out) :: lp !! standard log CDF
        real(dp) :: ls, sigma, root, hi, h, v, u, u2, phase, decay, wd, sd, sp
        integer :: n, i, j
        ld = negative_infinity(0.0_dp)
        lp = ld
        if (z < -500.0_dp) return
        ls = -0.5_dp*scifort_pi*z-log(0.5_dp*scifort_pi)-1.0_dp
        if (ls >= log(huge(1.0_dp))) return
        sigma = exp(ls)
        if (ls > log(1.0e12_dp)) then
            ! The first omitted relative correction is O(1/sigma).
            ld = -sigma+0.5_dp*ls+log(0.5_dp*scifort_pi)-scifort_log_sqrt_two_pi
            lp = -sigma-0.5_dp*ls-scifort_log_sqrt_two_pi
            return
        end if
        root = sqrt(sigma)
        hi = max(16.0_dp,80.0_dp/root)
        n = ceiling(hi/0.5_dp)
        h = hi/real(n,dp)
        sd = 0.0_dp
        sp = 0.0_dp
        do j = 1, n
            do i = 1, 16
                v = (real(j-1,dp)+0.5_dp*(1.0_dp+p(6+i)))*h
                u = v/root
                u2 = u*u
                decay = sigma*(0.5_dp*log1p_safe(u2)-u*atan(u))
                if (u < 0.01_dp) then
                    phase = sigma*u*u2*(1.0_dp/6.0_dp-u2*(1.0_dp/20.0_dp &
                        -u2*(1.0_dp/42.0_dp-u2/72.0_dp)))
                else
                    phase = sigma*(atan(u)+0.5_dp*u*log1p_safe(u2)-u)
                end if
                wd = p(22+i)*exp(decay)
                sd = sd+wd*cos(phase)
                sp = sp+wd*(cos(phase)+u*sin(phase))/(1.0_dp+u2)
            end do
        end do
        ld = -sigma+0.5_dp*ls-scifort_log_two+log(0.5_dp*h*sd)
        lp = -sigma-0.5_dp*ls-log(scifort_pi)+log(0.5_dp*h*sp)
    end subroutine ld_left

    pure function sp_pdf(x, p) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: p(:) !! packed specialist distribution parameters
        real(dp) :: value, lp, lq, z, lz, pp(38)
        if (p(1) == 1.0_dp) then
            call sp_coordinate(x,p,z,lz)
            pp = p
            if (x < p(2)) pp(5) = -p(5)
            value = nt_integral(lz,0,pp)-log(p(3))
        else
            call sp_logs(x,p,value,lp,lq)
        end if
    end function sp_pdf

    pure function sp_cdf(x, p) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: p(:) !! packed specialist distribution parameters
        real(dp) :: value, ld, lq
        call sp_logs(x,p,ld,value,lq)
    end function sp_cdf

    pure function sp_sf(x, p) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: p(:) !! packed specialist distribution parameters
        real(dp) :: value, ld, lp
        call sp_logs(x,p,ld,lp,value)
    end function sp_sf

    pure function sp_interval(a, b, p) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper interval endpoint, greater than a
        real(dp), intent(in) :: p(:) !! packed specialist distribution parameters
        real(dp) :: value, la, lb, qa, qb, da, db, lw, z, lz, width, t, pp(38), ld
        integer :: i
        call sp_logs(a,p,da,la,qa)
        call sp_logs(b,p,db,lb,qb)
        if (lb <= -scifort_log_two) then
            if (la < lb-0.01_dp .or. .not. ieee_is_finite(lb)) then
                value = log_difference(lb,la)
                if (.not. ieee_is_finite(lb)) value = lb
                return
            end if
        else if (qa <= -scifort_log_two) then
            if (qb < qa-0.01_dp) then
                value = log_difference(qa,qb)
                return
            end if
        else if (la < lb-0.01_dp) then
            value = log_difference(lb,la)
            return
        end if
        ! Positive density integration keeps adjacent physical endpoints distinct.
        lw = log_distance(b,a)-log(p(3))
        call sp_coordinate(a,p,z,lz)
        if (lz >= log(1.0e100_dp)) then
            ! Beyond this magnitude power-law corrections are negligible. Use
            ! physical distance ratios even when quadrature nodes would overflow.
            t = 1.0_dp
            if (p(1) == 1.0_dp) t = p(4)
            if (a > p(2)) then
                ld = log_distance(b,a)-log_distance(a,p(2))
                value = qa+log_exponential_mass(log(t)+burr_logsoftplus(ld))
            else if (b < p(2) .and. p(1) == 1.0_dp) then
                ld = log_distance(b,a)-log_distance(p(2),b)
                value = lb+log_exponential_mass(log(t)+burr_logsoftplus(ld))
            else
                value = negative_infinity(0.0_dp)
            end if
            return
        end if
        if (lw < log(tiny(1.0_dp))) then
            value = da+log_distance(b,a)
            return
        end if
        width = exp(lw)
        pp = p
        pp(2:3) = [0.0_dp,1.0_dp]
        value = negative_infinity(0.0_dp)
        do i = 1, 16
            t = z+0.5_dp*(1.0_dp+p(6+i))*width
            ld = sp_pdf(t,pp)
            value = bounded_logadd(value,ld+log(p(22+i)))
        end do
        value = min(0.0_dp,value+lw-scifort_log_two)
    end function sp_interval

    pure function nct_censored_loglikelihood(data, df, nc, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: df !! degrees of freedom in [0.5,100]
        real(dp), intent(in) :: nc !! noncentrality in [-10,10]
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(38)
        real(dp), allocatable :: nodes(:), weights(:)
        value = quiet_nan(0.0_dp)
        params = 0.0_dp
        params(1:6) = [1.0_dp,0.0_dp,1.0_dp,df,nc,0.0_dp]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        if (.not. all(ieee_is_finite(params))) return
        if (params(3) <= 0.0_dp) return
        if (df < 0.5_dp .or. df > 100.0_dp .or. abs(nc) > 10.0_dp) return
        params(6) = scifort_log_two+0.5_dp*df*log(0.5_dp*df)-log_gamma(0.5_dp*df)
        call roots_legendre(16,nodes,weights)
        params(7:22) = nodes
        params(23:38) = weights
        value = censored_loglikelihood(data,params,sp_pdf,sp_cdf,sp_sf,sp_interval)
    end function nct_censored_loglikelihood

    pure function nct_censored_nnlf(data, df, nc, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: df !! degrees of freedom in [0.5,100]
        real(dp), intent(in) :: nc !! noncentrality in [-10,10]
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -nct_censored_loglikelihood(data,df,nc,loc,scale)
    end function nct_censored_nnlf

    pure function landau_censored_loglikelihood(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(38)
        real(dp), allocatable :: nodes(:), weights(:)
        value = quiet_nan(0.0_dp)
        params = 0.0_dp
        params(1:6) = [2.0_dp,0.0_dp,1.0_dp,0.0_dp,0.0_dp,0.0_dp]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        if (.not. all(ieee_is_finite(params))) return
        if (params(3) <= 0.0_dp) return
        call roots_legendre(16,nodes,weights)
        params(7:22) = nodes
        params(23:38) = weights
        value = censored_loglikelihood(data,params,sp_pdf,sp_cdf,sp_sf,sp_interval)
    end function landau_censored_loglikelihood

    pure function landau_censored_nnlf(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -landau_censored_loglikelihood(data,loc,scale)
    end function landau_censored_nnlf


    ! Log-probability inversion of the Tukey quantile; see CODE_PROVENANCE.md.
    pure function tl_target(lq, lam, endpoint) result(v)
        implicit none
        real(dp), intent(in) :: lq !! log smaller probability, at most -log(2)
        real(dp), intent(in) :: lam !! nonzero shape with abs(lam) in [1e-6,5]
        logical, intent(in) :: endpoint !! positive-shape endpoint-gap representation
        real(dp) :: v, lp, delta, small
        lp = ig_complement(lq)
        if (endpoint) then
            small = lq
            if (lq > -36.0_dp) small = log(-lp)
            v = bounded_logadd(log_exponential_mass(log(lam)+small),lam*lq)-log(lam)
        else
            delta = max(0.0_dp,lp-lq)
            if (delta == 0.0_dp) then
                v = negative_infinity(0.0_dp)
            else if (lam > 0.0_dp) then
                v = lam*lp+ig_complement(-lam*delta)-log(lam)
            else
                v = lam*lq+ig_complement(lam*delta)-log(-lam)
            end if
        end if
    end function tl_target

    pure subroutine tl_logs(target, p, endpoint, lp, lq, ld)
        implicit none
        real(dp), intent(in) :: target !! log endpoint gap or log abs(z), selected by endpoint
        real(dp), intent(in) :: p(:) !! shape, location, scale, finite lower/upper bounds
        logical, intent(in) :: endpoint !! target is a positive-shape endpoint gap
        real(dp), intent(out) :: lp, lq, ld !! log large/small probabilities and standardized density
        real(dp) :: lo, hi, mid, q, lam
        integer :: k
        lam = p(1)
        if (lam == 0.0_dp) then
            if (target > log(huge(1.0_dp))) then
                lq = negative_infinity(0.0_dp)
            else
                lq = -log1pexp(exp(target))
            end if
        else
            lo = -1.0_dp
            hi = -scifort_log_two
            do k = 1, 40
                q = tl_target(lo,lam,endpoint)
                if (endpoint) then
                    if (q <= target) exit
                else
                    if (q >= target) exit
                end if
                lo = 2.0_dp*lo
            end do
            do k = 1, 100
                mid = lo+(hi-lo)/2.0_dp
                if (mid == lo .or. mid == hi) exit
                q = tl_target(mid,lam,endpoint)
                if ((endpoint .and. q < target) .or. &
                    (.not. endpoint .and. q > target)) then
                    lo = mid
                else
                    hi = mid
                end if
            end do
            lq = lo+(hi-lo)/2.0_dp
        end if
        lp = ig_complement(lq)
        ld = -bounded_logadd((lam-1.0_dp)*lp,(lam-1.0_dp)*lq)
    end subroutine tl_logs

    pure subroutine tl_at(x, p, lp, lq, ld)
        implicit none
        real(dp), intent(in) :: x !! finite observation or censoring limit
        real(dp), intent(in) :: p(:) !! shape, location, scale, finite lower/upper bounds
        real(dp), intent(out) :: lp, lq, ld !! physical log CDF, SF and density
        real(dp) :: target, temp
        logical :: endpoint
        lp = negative_infinity(0.0_dp)
        lq = 0.0_dp
        ld = lp
        if (p(1) > 0.0_dp) then
            if (x < p(4) .or. x > p(5)) then
                if (x > p(5)) then
                    lp = 0.0_dp
                    lq = negative_infinity(0.0_dp)
                end if
                return
            end if
            if (x == p(4) .or. x == p(5)) then
                if (x == p(5)) then
                    lp = 0.0_dp
                    lq = negative_infinity(0.0_dp)
                end if
                if (p(1) >= 1.0_dp) ld = -log(p(3))
                if (p(1) == 1.0_dp) ld = -scifort_log_two-log(p(3))
                return
            end if
        end if
        if (x == p(2)) then
            lp = -scifort_log_two
            lq = lp
            ld = (p(1)-2.0_dp)*scifort_log_two-log(p(3))
            return
        end if
        target = log_distance(max(x,p(2)),min(x,p(2)))-log(p(3))
        endpoint = .false.
        if (p(1) > 0.0_dp) then
            endpoint = target > -log(p(1))-scifort_log_two
            if (endpoint) then
                if (x < p(2)) then
                    target = log_distance(x,p(4))-log(p(3))
                else
                    target = log_distance(p(5),x)-log(p(3))
                end if
            end if
        end if
        call tl_logs(target,p,endpoint,lp,lq,ld)
        if (x < p(2)) then
            temp = lp
            lp = lq
            lq = temp
        end if
        ld = ld-log(p(3))
    end subroutine tl_at

    pure function tl_pdf(x, p) result(v)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: p(:) !! kernel parameters
        real(dp) :: v, lp, lq
        call tl_at(x,p,lp,lq,v)
    end function tl_pdf

    pure function tl_cdf(x, p) result(v)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: p(:) !! kernel parameters
        real(dp) :: v, lq, ld
        call tl_at(x,p,v,lq,ld)
    end function tl_cdf

    pure function tl_sf(x, p) result(v)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: p(:) !! kernel parameters
        real(dp) :: v, lp, ld
        call tl_at(x,p,lp,v,ld)
    end function tl_sf

    pure function tl_interval(a, b, p) result(v)
        implicit none
        real(dp), intent(in) :: a, b !! finite strictly ordered physical endpoints
        real(dp), intent(in) :: p(:) !! shape, location, scale, lower/upper bounds
        real(dp) :: v, left, right, lw, ld, lc, target, offset, t, weight, lp, lq
        real(dp) :: nodes(3), weights(3), ratio, z
        logical :: endpoint
        integer :: k
        left = a
        right = b
        if (p(1) > 0.0_dp) then
            left = max(a,p(4))
            right = min(b,p(5))
            if (right <= left) then
                v = negative_infinity(0.0_dp)
                return
            end if
            if (left == p(4)) then
                v = tl_cdf(right,p)
                return
            else if (right == p(5)) then
                v = tl_sf(left,p)
                return
            end if
        end if
        lw = log_distance(right,left)
        ratio = lw-log(p(3))
        if (p(1) > 0.0_dp) then
            ratio = max(ratio,lw-log_distance(left,p(4)),lw-log_distance(p(5),right))
        else if (left > p(2)) then
            ratio = lw-bounded_logadd(log(p(3)),log_distance(left,p(2)))
        else if (right < p(2)) then
            ratio = lw-bounded_logadd(log(p(3)),log_distance(p(2),right))
        end if
        if (ratio < log(1.0e-4_dp)) then
            nodes = [0.5_dp,(1.0_dp-sqrt(3.0_dp/5.0_dp))/2.0_dp, &
                (1.0_dp+sqrt(3.0_dp/5.0_dp))/2.0_dp]
            weights = [4.0_dp/9.0_dp,5.0_dp/18.0_dp,5.0_dp/18.0_dp]
            v = negative_infinity(0.0_dp)
            do k = 1, 3
                t = nodes(k)
                weight = weights(k)
                offset = lw+log(t)
                endpoint = .false.
                if (left >= p(2) .or. right <= p(2)) then
                    if (left >= p(2)) then
                        ld = negative_infinity(0.0_dp)
                        if (left > p(2)) ld = log_distance(left,p(2))
                        target = bounded_logadd(ld,offset)-log(p(3))
                    else
                        ld = log_distance(p(2),left)
                        target = log_difference(ld,offset)-log(p(3))
                    end if
                    if (p(1) > 0.0_dp) then
                        endpoint = target > -log(p(1))-scifort_log_two
                        if (endpoint) then
                            if (left >= p(2)) then
                                ld = log_distance(p(5),left)
                                target = log_difference(ld,offset)-log(p(3))
                            else
                                ld = log_distance(left,p(4))
                                target = bounded_logadd(ld,offset)-log(p(3))
                            end if
                        end if
                    end if
                else
                    z = location_scale_z(left,p(2:3))+t*exp(lw-log(p(3)))
                    target = bounded_log(abs(z))
                end if
                call tl_logs(target,p,endpoint,lp,lq,lc)
                v = bounded_logadd(v,log(weight)+lc)
            end do
            v = v+lw-log(p(3))
        else
            v = tail_interval(left,right,p,tl_cdf,tl_sf)
            if (v > -scifort_log_two) then
                lp = tl_cdf(left,p)
                lq = tl_sf(right,p)
                v = ig_complement(bounded_logadd(lp,lq))
            end if
        end if
    end function tl_interval

    ! Derivative of the finite positive Smirnov survival sum.
    pure function ko_regular(z, n) result(v)
        implicit none
        real(dp), intent(in) :: z !! standardized interior coordinate in [1/n,1-1/n]
        integer, intent(in) :: n !! sample size in [2,100]
        real(dp) :: v, total, correction, term, a, b, g, logt, next, adjusted
        integer :: j, m
        total = 0.0_dp
        correction = 0.0_dp
        m = min(n-2,floor(real(n,dp)*(1.0_dp-z)))
        do j = 0, m
            a = z+real(j,dp)/real(n,dp)
            b = 1.0_dp-z-real(j,dp)/real(n,dp)
            if (b <= 0.0_dp) cycle
            if (j == 0) then
                logt = real(n,dp)*log(b)
                g = -real(n,dp)/b
            else
                logt = log(z)+log_gamma(real(n+1,dp))-log_gamma(real(j+1,dp)) &
                    -log_gamma(real(n-j+1,dp))+real(j-1,dp)*log(a) &
                    +real(n-j,dp)*log(b)
                g = 1.0_dp/z+real(j-1,dp)/a-real(n-j,dp)/b
            end if
            term = -exp(logt)*g
            adjusted = term-correction
            next = total+adjusted
            correction = (next-total)-adjusted
            total = next
        end do
        v = bounded_log(max(0.0_dp,total))
    end function ko_regular

    pure function ko_density(z, lu, p, side) result(v)
        implicit none
        real(dp), intent(in) :: z !! standardized coordinate in (0,1)
        real(dp), intent(in) :: lu !! log standardized distance to the upper endpoint
        real(dp), intent(in) :: p(:) !! n, location, scale, upper endpoint
        integer, intent(in) :: side !! at first knot: -1 left, +1 right, 0 midpoint
        real(dp) :: v, low, high
        integer :: n
        n = int(p(1))
        if (n == 1) then
            v = 0.0_dp
        else if (z < 1.0_dp/p(1)) then
            v = real(n-2,dp)*log1p_safe(z)+log1p_safe(p(1)*z)
        else if (z == 1.0_dp/p(1)) then
            low = real(n-2,dp)*log1p_safe(z)+log1p_safe(p(1)*z)
            high = ko_regular(z,n)
            v = low
            if (side > 0) v = high
            if (side == 0) v = bounded_logadd(low,high)-scifort_log_two
        else if (z >= 1.0_dp-1.0_dp/p(1)) then
            v = log(p(1))+real(n-1,dp)*lu
        else
            v = ko_regular(z,n)
        end if
    end function ko_density

    pure subroutine ko_tails(x, p, lp, lq)
        implicit none
        real(dp), intent(in) :: x !! finite physical censoring limit
        real(dp), intent(in) :: p(:) !! n, location, scale, upper endpoint
        real(dp), intent(out) :: lp, lq !! complementary log probabilities
        real(dp) :: lz, lu, z
        lp = negative_infinity(0.0_dp)
        lq = 0.0_dp
        if (x <= p(2)) return
        if (x >= p(4)) then
            lp = 0.0_dp
            lq = negative_infinity(0.0_dp)
            return
        end if
        lz = log_distance(x,p(2))-log(p(3))
        lu = log_distance(p(4),x)-log(p(3))
        z = exp(lz)
        if (p(1) == 1.0_dp) then
            if (lz < lu) then
                lp = lz
                lq = ig_complement(lp)
            else
                lq = lu
                lp = ig_complement(lq)
            end if
        else if (lz <= -log(p(1))) then
            lp = lz+(p(1)-1.0_dp)*log1p_safe(z)
            lq = ig_complement(lp)
        else if (lu <= -log(p(1))) then
            lq = p(1)*lu
            lp = ig_complement(lq)
        else
            lq = ksone_logsf(z,p(1))
            lp = ig_complement(lq)
        end if
    end subroutine ko_tails

    pure function ko_pdf(x, p) result(v)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: p(:) !! n, location, scale, upper endpoint
        real(dp) :: v, z, lu, knot
        integer :: side
        v = negative_infinity(0.0_dp)
        if (x < p(2) .or. x > p(4)) return
        if (p(1) == 1.0_dp) then
            v = -log(p(3))
            return
        end if
        if (x == p(2) .or. x == p(4)) return
        z = location_scale_z(x,p(2:3))
        lu = log_distance(p(4),x)-log(p(3))
        knot = p(2)+p(3)/p(1)
        side = -1
        if (x > knot) side = 1
        if (x == knot) then
            side = 0
            z = 1.0_dp/p(1)
        end if
        v = ko_density(z,lu,p,side)-log(p(3))
    end function ko_pdf

    pure function ko_cdf(x, p) result(v)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: p(:) !! kernel parameters
        real(dp) :: v, lq
        call ko_tails(x,p,v,lq)
    end function ko_cdf

    pure function ko_sf(x, p) result(v)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: p(:) !! kernel parameters
        real(dp) :: v, lp
        call ko_tails(x,p,lp,v)
    end function ko_sf

    pure function ko_small(a, b, p, side) result(v)
        implicit none
        real(dp), intent(in) :: a, b !! finite ordered endpoints inside support
        real(dp), intent(in) :: p(:) !! kernel parameters
        integer, intent(in) :: side !! selected side at the first density knot
        real(dp) :: v, z, za, lw, lu, offset, t, term, weight, nodes(3)
        integer :: k
        nodes = [0.5_dp,(1.0_dp-sqrt(3.0_dp/5.0_dp))/2.0_dp, &
            (1.0_dp+sqrt(3.0_dp/5.0_dp))/2.0_dp]
        za = location_scale_z(a,p(2:3))
        lw = log_distance(b,a)
        v = negative_infinity(0.0_dp)
        do k = 1, 3
            t = nodes(k)
            weight = 5.0_dp/18.0_dp
            if (k == 1) weight = 4.0_dp/9.0_dp
            offset = lw+log(t)
            z = za+t*exp(lw-log(p(3)))
            lu = log_difference(log_distance(p(4),a),offset)-log(p(3))
            if (abs(z-1.0_dp/p(1)) < 64.0_dp*epsilon(1.0_dp)) z = 1.0_dp/p(1)
            term = ko_density(z,lu,p,side)
            v = bounded_logadd(v,log(weight)+term)
        end do
        v = v+lw-log(p(3))
    end function ko_small

    pure function ko_interval(a, b, p) result(v)
        implicit none
        real(dp), intent(in) :: a, b !! finite strictly ordered physical endpoints
        real(dp), intent(in) :: p(:) !! n, location, scale, upper endpoint
        real(dp) :: v, left, right, la, lb, ua, ub, lw, delta, knot, lp, lq
        integer :: side
        left = max(a,p(2))
        right = min(b,p(4))
        v = negative_infinity(0.0_dp)
        if (right <= left) return
        if (left == p(2)) then
            v = ko_cdf(right,p)
            return
        else if (right == p(4)) then
            v = ko_sf(left,p)
            return
        end if
        lw = log_distance(right,left)
        if (p(1) == 1.0_dp) then
            v = lw-log(p(3))
            return
        end if
        la = log_distance(left,p(2))
        lb = log_distance(right,p(2))
        ua = log_distance(p(4),left)
        ub = log_distance(p(4),right)
        if (lb-log(p(3)) <= -log(p(1))) then
            delta = bounded_logadd(burr_logsoftplus(lw-la),log(p(1)-1.0_dp) &
                +burr_logsoftplus(lw-log(p(3))-log1p_safe(exp(la-log(p(3))))))
            v = ko_cdf(right,p)+log_exponential_mass(delta)
        else if (ua-log(p(3)) <= -log(p(1))) then
            delta = log(p(1))+burr_logsoftplus(lw-ub)
            v = ko_sf(left,p)+log_exponential_mass(delta)
        else if (lw-log(p(3)) < log(1.0e-5_dp)) then
            knot = p(2)+p(3)/p(1)
            if (left < knot .and. right > knot) then
                v = bounded_logadd(ko_small(left,knot,p,-1),ko_small(knot,right,p,1))
            else
                side = -1
                if (left >= knot) side = 1
                v = ko_small(left,right,p,side)
            end if
        else
            v = tail_interval(left,right,p,ko_cdf,ko_sf)
            if (v > -scifort_log_two) then
                lp = ko_cdf(left,p)
                lq = ko_sf(right,p)
                v = ig_complement(bounded_logadd(lp,lq))
            end if
        end if
    end function ko_interval

    pure function tukeylambda_censored_loglikelihood(data, lam, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: lam !! shape in [-5,5], zero or abs(lam) >= 1e-6
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(5), extent
        value = quiet_nan(0.0_dp)
        params = 0.0_dp
        params(:3) = [lam,0.0_dp,1.0_dp]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        if (.not. all(ieee_is_finite(params(:3)))) return
        if (params(3) <= 0.0_dp) return
        if (abs(lam) > 5.0_dp) return
        if (lam /= 0.0_dp .and. abs(lam) < 1.0e-6_dp) return
        if (lam > 0.0_dp) then
            if (params(3) < tiny(1.0_dp)) return
            if (params(2) /= 0.0_dp) then
            if (log(abs(params(2)))-log(params(3)) > log(10.0_dp)) return
        end if
            if (log(params(3))-log(lam) > log(huge(1.0_dp)/4.0_dp)) return
            extent = params(3)/lam
            if (abs(params(2)) > huge(1.0_dp)-extent) return
            params(4) = params(2)-extent
            params(5) = params(2)+extent
            if (.not. all(ieee_is_finite(params(4:5)))) return
            if (params(4) >= params(2) .or. params(5) <= params(2)) return
        end if
        value = censored_loglikelihood(data,params,tl_pdf,tl_cdf,tl_sf,tl_interval)
    end function tukeylambda_censored_loglikelihood

    pure function tukeylambda_censored_nnlf(data, lam, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: lam !! shape in [-5,5], zero or abs(lam) >= 1e-6
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -tukeylambda_censored_loglikelihood(data,lam,loc,scale)
    end function tukeylambda_censored_nnlf

    pure function ksone_censored_loglikelihood(data, n, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: n !! integer sample size in [1,100]
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(4), extent
        value = quiet_nan(0.0_dp)
        params = 0.0_dp
        params(:3) = [n,0.0_dp,1.0_dp]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        if (.not. all(ieee_is_finite(params(:3)))) return
        if (params(3) <= 0.0_dp) return
        if (n < 1.0_dp .or. n > 100.0_dp) return
        if (n /= real(nint(n),dp)) return
        if (params(3) < tiny(1.0_dp) .or. params(3) > huge(1.0_dp)/4.0_dp) return
        if (abs(params(2))/params(3) > 10.0_dp) return
        if (abs(params(2)) > huge(1.0_dp)-params(3)) return
        params(4) = params(2)+params(3)
        if (.not. ieee_is_finite(params(4))) return
        if (params(4) <= params(2)) return
        value = censored_loglikelihood(data,params,ko_pdf,ko_cdf,ko_sf,ko_interval)
    end function ksone_censored_loglikelihood

    pure function ksone_censored_nnlf(data, n, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: n !! integer sample size in [1,100]
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -ksone_censored_loglikelihood(data,n,loc,scale)
    end function ksone_censored_nnlf

    pure function vonmises_censored_loglikelihood(data, kappa, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: kappa !! concentration in [0,100]
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        ! SciFort extension: ordinary censoring on the fixed affine [-pi,pi] branch.
        value = vonmises_line_censored_loglikelihood(data,kappa,loc,scale)
    end function vonmises_censored_loglikelihood

    pure function vonmises_censored_nnlf(data, kappa, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: kappa !! concentration in [0,100]
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -vonmises_censored_loglikelihood(data,kappa,loc,scale)
    end function vonmises_censored_nnlf


    ! Positive transformed-density quadrature from the cited mathematical PDFs.
    pure function bk_density(z, p) result(v)
        implicit none
        real(dp), intent(in) :: z !! log standardized x for GIG, standardized x otherwise
        real(dp), intent(in) :: p(:) !! kernel, loc, scale, shapes, normalizer, GL rule
        real(dp) :: v, cost, r, lr, arg, lk, gap, factor
        v = negative_infinity(0.0_dp)
        if (.not. ieee_is_finite(z)) return
        if (p(1) == 1.0_dp) then
            arg = log(p(5)/2.0_dp)+abs(z)
            if (arg > log(huge(1.0_dp))) return
            factor = exp(abs(z)/2.0_dp)
            cost = (p(5)/2.0_dp)*factor
            if (cost > 1.0_dp) then
                if (factor > huge(1.0_dp)/cost) return
            end if
            cost = cost*factor+exp(log(p(5)/2.0_dp)-abs(z))
            v = p(7)+(p(4)-1.0_dp)*z-cost
        else
            r = hypot(1.0_dp,z)
            lr = log(r)
            arg = log(p(5))+lr
            if (arg > log(1.0e100_dp)) then
                lk = 0.5_dp*(log(scifort_pi/2.0_dp)-arg)
            else
                lk = log_besselk_scaled(p(4)-0.5_dp,exp(arg))
            end if
            gap = p(5)-p(6)
            if (z < 0.0_dp) gap = p(5)+p(6)
            if (gap > 1.0_dp) then
                if (abs(z) > huge(1.0_dp)/gap) return
            end if
            cost = gap*abs(z)+(p(5)/r)/(1.0_dp+abs(z)/r)
            v = p(7)+lk+(p(4)-0.5_dp)*lr-cost
        end if
    end function bk_density

    pure function bk_coordinate(x, p) result(z)
        implicit none
        real(dp), intent(in) :: x !! physical coordinate inside the support
        real(dp), intent(in) :: p(:) !! kernel parameters
        real(dp) :: z
        if (p(1) == 1.0_dp) then
            z = log_distance(x,p(2))-log(p(3))
        else
            z = location_scale_z(x,p(2:3))
        end if
    end function bk_coordinate

    pure function bk_pdf(x, p) result(v)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: p(:) !! kernel parameters
        real(dp) :: v
        v = negative_infinity(0.0_dp)
        if (p(1) == 1.0_dp .and. x <= p(2)) return
        v = bk_density(bk_coordinate(x,p),p)-log(p(3))
    end function bk_pdf

    pure function bk_transformed(y, p) result(v)
        implicit none
        real(dp), intent(in) :: y !! log or asinh standardized coordinate
        real(dp), intent(in) :: p(:) !! kernel parameters
        real(dp) :: v
        if (p(1) == 1.0_dp) then
            v = bk_density(y,p)+y
        else
            v = bk_density(sinh(y),p)+log(cosh(y))
        end if
    end function bk_transformed

    pure function bk_integral(y, direction, p) result(v)
        implicit none
        real(dp), intent(in) :: y !! finite transformed tail endpoint
        real(dp), intent(in) :: direction !! -1 for lower, +1 for upper tail
        real(dp), intent(in) :: p(:) !! kernel parameters including the GL rule
        real(dp) :: v, far, near, mid, base, width, node, half, term
        integer :: k, j, segments
        base = bk_transformed(y,p)
        far = direction*max(12.0_dp,abs(y)+1.0_dp)
        ! Resolve the endpoint-relative tail even when its width is very small.
        block
            near = y
            do k = 1, 48
                mid = near+(far-near)/2.0_dp
                if (bk_transformed(mid,p) < base-50.0_dp) then
                    far = mid
                else
                    near = mid
                end if
            end do
        end block
        width = abs(far-y)
        segments = max(16,ceiling(width/0.125_dp))
        half = width/(2.0_dp*real(segments,dp))
        v = negative_infinity(0.0_dp)
        do j = 1, segments
            mid = y+direction*(real(2*j-1,dp)*half)
            do k = 1, 16
                node = mid+half*p(7+k)
                term = log(half*p(23+k))+bk_transformed(node,p)
                v = bounded_logadd(v,term)
            end do
        end do
        v = min(0.0_dp,v)
    end function bk_integral

    pure subroutine bk_tails(x, p, lp, lq)
        implicit none
        real(dp), intent(in) :: x !! finite censoring limit
        real(dp), intent(in) :: p(:) !! kernel parameters
        real(dp), intent(out) :: lp, lq !! logarithms of complementary probabilities
        real(dp) :: z, y, ld, gap
        lp = negative_infinity(0.0_dp)
        lq = 0.0_dp
        if (p(1) == 1.0_dp .and. x <= p(2)) return
        z = bk_coordinate(x,p)
        if (p(1) == 1.0_dp) then
            y = z
            if (abs(y) > log(1.0e8_dp)) then
                ld = bk_density(y,p)
                if (y < 0.0_dp) then
                    lp = ld+2.0_dp*y-log(p(5)/2.0_dp)
                    lq = ig_complement(lp)
                else
                    lq = ld-log(p(5)/2.0_dp)
                    lp = ig_complement(lq)
                end if
                return
            end if
        else
            if (abs(z) > 1.0e8_dp) then
                ld = bk_density(z,p)
                gap = p(5)-p(6)
                if (z < 0.0_dp) gap = p(5)+p(6)
                if (z < 0.0_dp) then
                    lp = ld-log(gap)
                    lq = ig_complement(lp)
                else
                    lq = ld-log(gap)
                    lp = ig_complement(lq)
                end if
                return
            end if
            y = asinh(z)
        end if
        if (y < -6.0_dp) then
            lp = bk_integral(y,-1.0_dp,p)
            lq = ig_complement(lp)
        else if (y > 6.0_dp) then
            lq = bk_integral(y,1.0_dp,p)
            lp = ig_complement(lq)
        else
            lp = bk_integral(y,-1.0_dp,p)
            lq = bk_integral(y,1.0_dp,p)
            if (lp < lq) then
                lq = ig_complement(lp)
            else
                lp = ig_complement(lq)
            end if
        end if
    end subroutine bk_tails

    pure function bk_cdf(x, p) result(v)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: p(:) !! kernel parameters
        real(dp) :: v, lq
        call bk_tails(x,p,v,lq)
    end function bk_cdf

    pure function bk_sf(x, p) result(v)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: p(:) !! kernel parameters
        real(dp) :: v, lp
        call bk_tails(x,p,lp,v)
    end function bk_sf

    pure function bk_interval(a, b, p) result(v)
        implicit none
        real(dp), intent(in) :: a, b !! finite strictly ordered physical endpoints
        real(dp), intent(in) :: p(:) !! kernel parameters
        real(dp) :: v, z, zb, lw, slope, t, offset, term, weight, lp, lq
        integer :: k
        v = negative_infinity(0.0_dp)
        if (p(1) == 1.0_dp) then
            if (b <= p(2)) return
            if (a <= p(2)) then
                v = bk_cdf(b,p)
                return
            end if
        end if
        z = bk_coordinate(a,p)
        lw = log_distance(b,a)-log(p(3))
        zb = bk_coordinate(b,p)
        ! Retain the physical increment when endpoint tail logs round together.
        if (p(1) == 1.0_dp) then
            if (z > log(1.0e8_dp)) then
                call bk_tails(a,p,lp,lq)
                v = lq+log_exponential_mass(lw+log(p(5)/2.0_dp))
                return
            else if (zb < -log(1.0e8_dp)) then
                call bk_tails(b,p,lp,lq)
                v = lp+log_exponential_mass(log(p(5)/2.0_dp)+lw-z-zb)
                return
            end if
        else
            if (z > 1.0e8_dp) then
                call bk_tails(a,p,lp,lq)
                v = lq+log_exponential_mass(lw+log(p(5)-p(6)))
                return
            else if (zb < -1.0e8_dp) then
                call bk_tails(b,p,lp,lq)
                v = lp+log_exponential_mass(lw+log(p(5)+p(6)))
                return
            end if
        end if
        slope = log(100.0_dp)
        if (p(1) == 1.0_dp) slope = log(100.0_dp)+max(0.0_dp,-2.0_dp*z)
        if (lw+slope < log(0.01_dp)) then
            v = negative_infinity(0.0_dp)
            do k = 1, 3
                t = 0.5_dp
                weight = 4.0_dp/9.0_dp
                if (k /= 1) then
                    t = (1.0_dp+real(2*k-5,dp)*sqrt(3.0_dp/5.0_dp))/2.0_dp
                    weight = 5.0_dp/18.0_dp
                end if
                if (p(1) == 1.0_dp) then
                    offset = bounded_logadd(z,lw+log(t))
                else
                    offset = z+t*exp(lw)
                end if
                term = bk_density(offset,p)+log(weight)
                v = bounded_logadd(v,term)
            end do
            v = v+lw
        else
            v = tail_interval(a,b,p,bk_cdf,bk_sf)
            if (v > -scifort_log_two) then
                lp = bk_cdf(a,p)
                lq = bk_sf(b,p)
                v = ig_complement(bounded_logadd(lp,lq))
            end if
        end if
    end function bk_interval

    pure function geninvgauss_censored_loglikelihood(data, p, b, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: p !! order shape in [-5,5]
        real(dp), intent(in) :: b !! positive shape in [0.5,10]
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(39), delta
        real(dp), allocatable :: nodes(:), weights(:)
        integer :: status
        value = quiet_nan(0.0_dp)
        params = 0.0_dp
        params(:6) = [1.0_dp,0.0_dp,1.0_dp,p,b,0.0_dp]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        if (.not. all(ieee_is_finite(params(:6)))) return
        if (params(3) <= 0.0_dp) return
        if (abs(params(4)) > 5.0_dp) return
        if (params(5) < 0.5_dp .or. params(5) > 10.0_dp) return
        if (abs(params(6)) > 0.8_dp*params(5)) return
        if (params(1) == 1.0_dp) then
            params(7) = -scifort_log_two-log_besselk_scaled(params(4),params(5))+params(5)
        else
            delta = sqrt((params(5)-params(6))*(params(5)+params(6)))
            params(7) = params(4)*log(delta)-scifort_log_sqrt_two_pi &
                -(params(4)-0.5_dp)*log(params(5))-log_besselk_scaled(params(4),delta)+delta
        end if
        call roots_legendre(16,nodes,weights,status=status)
        if (status /= 0) return
        params(8:23) = nodes
        params(24:39) = weights
        value = censored_loglikelihood(data,params,bk_pdf,bk_cdf,bk_sf,bk_interval)
    end function geninvgauss_censored_loglikelihood

    pure function geninvgauss_censored_nnlf(data, p, b, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: p !! order shape in [-5,5]
        real(dp), intent(in) :: b !! positive shape in [0.5,10]
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -geninvgauss_censored_loglikelihood(data,p, b,loc,scale)
    end function geninvgauss_censored_nnlf

    pure function norminvgauss_censored_loglikelihood(data, a, b, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: a !! tail shape in [0.5,10]
        real(dp), intent(in) :: b !! skew shape with abs(b) <= 0.8*a
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(39), delta
        real(dp), allocatable :: nodes(:), weights(:)
        integer :: status
        value = quiet_nan(0.0_dp)
        params = 0.0_dp
        params(:6) = [2.0_dp,0.0_dp,1.0_dp,-0.5_dp,a,b]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        if (.not. all(ieee_is_finite(params(:6)))) return
        if (params(3) <= 0.0_dp) return
        if (abs(params(4)) > 5.0_dp) return
        if (params(5) < 0.5_dp .or. params(5) > 10.0_dp) return
        if (abs(params(6)) > 0.8_dp*params(5)) return
        if (params(1) == 1.0_dp) then
            params(7) = -scifort_log_two-log_besselk_scaled(params(4),params(5))+params(5)
        else
            delta = sqrt((params(5)-params(6))*(params(5)+params(6)))
            params(7) = params(4)*log(delta)-scifort_log_sqrt_two_pi &
                -(params(4)-0.5_dp)*log(params(5))-log_besselk_scaled(params(4),delta)+delta
        end if
        call roots_legendre(16,nodes,weights,status=status)
        if (status /= 0) return
        params(8:23) = nodes
        params(24:39) = weights
        value = censored_loglikelihood(data,params,bk_pdf,bk_cdf,bk_sf,bk_interval)
    end function norminvgauss_censored_loglikelihood

    pure function norminvgauss_censored_nnlf(data, a, b, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: a !! tail shape in [0.5,10]
        real(dp), intent(in) :: b !! skew shape with abs(b) <= 0.8*a
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -norminvgauss_censored_loglikelihood(data,a, b,loc,scale)
    end function norminvgauss_censored_nnlf

    pure function genhyperbolic_censored_loglikelihood(data, p, a, b, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: p !! order shape in [-5,5]
        real(dp), intent(in) :: a !! tail shape in [0.5,10]
        real(dp), intent(in) :: b !! skew shape with abs(b) <= 0.8*a
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(39), delta
        real(dp), allocatable :: nodes(:), weights(:)
        integer :: status
        value = quiet_nan(0.0_dp)
        params = 0.0_dp
        params(:6) = [2.0_dp,0.0_dp,1.0_dp,p,a,b]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        if (.not. all(ieee_is_finite(params(:6)))) return
        if (params(3) <= 0.0_dp) return
        if (abs(params(4)) > 5.0_dp) return
        if (params(5) < 0.5_dp .or. params(5) > 10.0_dp) return
        if (abs(params(6)) > 0.8_dp*params(5)) return
        if (params(1) == 1.0_dp) then
            params(7) = -scifort_log_two-log_besselk_scaled(params(4),params(5))+params(5)
        else
            delta = sqrt((params(5)-params(6))*(params(5)+params(6)))
            params(7) = params(4)*log(delta)-scifort_log_sqrt_two_pi &
                -(params(4)-0.5_dp)*log(params(5))-log_besselk_scaled(params(4),delta)+delta
        end if
        call roots_legendre(16,nodes,weights,status=status)
        if (status /= 0) return
        params(8:23) = nodes
        params(24:39) = weights
        value = censored_loglikelihood(data,params,bk_pdf,bk_cdf,bk_sf,bk_interval)
    end function genhyperbolic_censored_loglikelihood

    pure function genhyperbolic_censored_nnlf(data, p, a, b, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: p !! order shape in [-5,5]
        real(dp), intent(in) :: a !! tail shape in [0.5,10]
        real(dp), intent(in) :: b !! skew shape with abs(b) <= 0.8*a
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -genhyperbolic_censored_loglikelihood(data,p, a, b,loc,scale)
    end function genhyperbolic_censored_nnlf


    ! Log-domain kernels from the mathematical specifications in CODE_PROVENANCE.md.
    pure subroutine kt_logs(lz, lp, lq, ld)
        implicit none
        real(dp), intent(in) :: lz !! log positive standardized coordinate
        real(dp), intent(out) :: lp, lq, ld !! log CDF, SF, and standardized density
        real(dp) :: a, la, total, weighted, term, m, signv, power
        integer :: k
        lp = negative_infinity(0.0_dp)
        lq = 0.0_dp
        ld = lp
        if (lz < log(0.82_dp)) then
            la = log(scifort_pi**2/8.0_dp)-2.0_dp*lz
            if (la > log(huge(1.0_dp))) return
            a = exp(la)
            total = 1.0_dp
            weighted = 1.0_dp
            if (a < 100.0_dp) then
                do k = 2, 16
                    m = real(2*k-1,dp)**2
                    term = exp(-(m-1.0_dp)*a)
                    total = total+term
                    weighted = weighted+term*(m-exp(-scifort_log_two-la)) &
                        /(1.0_dp-exp(-scifort_log_two-la))
                end do
            end if
            lp = scifort_log_sqrt_two_pi-lz-a+log(total)
            ld = scifort_log_sqrt_two_pi-2.0_dp*lz-a+scifort_log_two+la &
                +log1p_safe(-exp(-scifort_log_two-la))+log(weighted)
            lq = ig_complement(lp)
        else
            la = scifort_log_two+2.0_dp*lz
            lp = 0.0_dp
            lq = negative_infinity(0.0_dp)
            ld = lq
            if (la > log(huge(1.0_dp))) return
            a = exp(la)
            total = 1.0_dp
            weighted = 1.0_dp
            if (a < 100.0_dp) then
                do k = 2, 16
                    power = real(k,dp)**2
                    signv = 1.0_dp
                    if (mod(k,2) == 0) signv = -1.0_dp
                    term = signv*exp(-(power-1.0_dp)*a)
                    total = total+term
                    weighted = weighted+power*term
                end do
            end if
            lq = scifort_log_two-a+log(total)
            ld = log(8.0_dp)+lz-a+log(weighted)
            lp = ig_complement(lq)
        end if
    end subroutine kt_logs

    pure function kt_pdf(x, p) result(v)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: p(:) !! finite location and positive scale
        real(dp) :: v, lp, lq
        v = negative_infinity(0.0_dp)
        if (x <= p(1)) return
        call kt_logs(log_distance(x,p(1))-log(p(2)),lp,lq,v)
        v = v-log(p(2))
    end function kt_pdf

    pure function kt_cdf(x, p) result(v)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: p(:) !! finite location and positive scale
        real(dp) :: v, lq, ld
        v = negative_infinity(0.0_dp)
        if (x <= p(1)) return
        call kt_logs(log_distance(x,p(1))-log(p(2)),v,lq,ld)
    end function kt_cdf

    pure function kt_sf(x, p) result(v)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: p(:) !! finite location and positive scale
        real(dp) :: v, lp, ld
        v = 0.0_dp
        if (x <= p(1)) return
        call kt_logs(log_distance(x,p(1))-log(p(2)),lp,v,ld)
    end function kt_sf

    pure function kt_interval(a, b, p) result(v)
        implicit none
        real(dp), intent(in) :: a !! finite lower endpoint
        real(dp), intent(in) :: b !! finite upper endpoint greater than a
        real(dp), intent(in) :: p(:) !! finite location and positive scale
        real(dp) :: v, lza, lzb, lr, la, delta, width, lp, lq
        v = negative_infinity(0.0_dp)
        if (b <= p(1)) return
        if (a <= p(1)) then
            v = kt_cdf(b,p)
            return
        end if
        lza = log_distance(a,p(1))-log(p(2))
        lzb = log_distance(b,p(1))-log(p(2))
        width = log_distance(b,a)
        lr = burr_logsoftplus(width-log_distance(a,p(1)))
        if (lzb <= log(0.1_dp)) then
            ! F(z)=sqrt(2*pi)/z*exp(-pi*pi/(8*z*z)); remaining terms < exp(-980).
            la = log(scifort_pi**2/8.0_dp)-2.0_dp*lza
            delta = la+log_exponential_mass(scifort_log_two+lr)
            delta = delta+ig_complement(lr-delta)
            v = kt_cdf(b,p)+log_exponential_mass(delta)
        else if (lza >= log(4.0_dp)) then
            ! The omitted SF terms are below exp(-96) relative to the first term.
            delta = scifort_log_two+2.0_dp*lza+log_expm1_from_log(scifort_log_two+lr)
            v = kt_sf(a,p)+log_exponential_mass(delta)
        else
            la = max(scifort_log_two+2.0_dp*lzb, &
                log(scifort_pi**2/8.0_dp)-2.0_dp*lza)
            if (width-log_distance(a,p(1))+max(0.0_dp,la)+log(8.0_dp) < log(0.01_dp)) then
                v = tt_small_mass(a,b,p,1)
            else
                v = tail_interval(a,b,p,kt_cdf,kt_sf)
                if (v > -scifort_log_two) then
                    lp = kt_cdf(a,p)
                    lq = kt_sf(b,p)
                    v = ig_complement(bounded_logadd(lp,lq))
                end if
            end if
        end if
    end function kt_interval

    pure function wc_small(ld, c) result(v)
        implicit none
        real(dp), intent(in) :: ld !! log angle distance in (0,pi]
        real(dp), intent(in) :: c !! concentration in (0,0.99]
        real(dp) :: v, t
        t = ld-scifort_log_two
        if (ld > -20.0_dp) t = log(tan(exp(min(ld,log(scifort_pi)))/2.0_dp))
        t = t+log1p_safe(c)-log1p_safe(-c)
        if (t < -20.0_dp) then
            v = t-log(scifort_pi)
        else if (t > 20.0_dp) then
            v = log(scifort_pi/2.0_dp)+log1p_safe(-2.0_dp*atan(exp(-t))/scifort_pi) &
                -log(scifort_pi)
        else
            v = log(atan(exp(t)))-log(scifort_pi)
        end if
    end function wc_small

    pure function wc_log_density(ld, p) result(v)
        implicit none
        real(dp), intent(in) :: ld !! log standardized distance to the closer support endpoint
        real(dp), intent(in) :: p(:) !! concentration, location, scale, upper support
        real(dp) :: v, ls, denominator
        ls = ld-scifort_log_two
        if (ld > -20.0_dp) ls = log(sin(exp(ld)/2.0_dp))
        denominator = bounded_logadd(2.0_dp*log1p_safe(-p(1)), &
            log(4.0_dp)+log(p(1))+2.0_dp*ls)
        v = log1p_safe(-p(1))+log1p_safe(p(1))-log(2.0_dp*scifort_pi) &
            -denominator-log(p(3))
    end function wc_log_density

    pure function wc_pdf(x, p) result(v)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: p(:) !! concentration, location, scale, finite upper support
        real(dp) :: v, ld
        v = negative_infinity(0.0_dp)
        if (x < p(2) .or. x > p(4)) return
        ld = negative_infinity(0.0_dp)
        if (x > p(2) .and. x < p(4)) ld = &
            min(log_distance(x,p(2)),log_distance(p(4),x))-log(p(3))
        v = wc_log_density(ld,p)
    end function wc_pdf

    pure function wc_cdf(x, p) result(v)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: p(:) !! concentration, location, scale, upper support
        real(dp) :: v, left, right
        v = negative_infinity(0.0_dp)
        if (x <= p(2)) return
        v = 0.0_dp
        if (x >= p(4)) return
        left = log_distance(x,p(2))-log(p(3))
        right = log_distance(p(4),x)-log(p(3))
        if (left <= right) then
            v = wc_small(left,p(1))
        else
            v = ig_complement(wc_small(right,p(1)))
        end if
    end function wc_cdf

    pure function wc_sf(x, p) result(v)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: p(:) !! concentration, location, scale, upper support
        real(dp) :: v, left, right
        v = 0.0_dp
        if (x <= p(2)) return
        v = negative_infinity(0.0_dp)
        if (x >= p(4)) return
        left = log_distance(x,p(2))-log(p(3))
        right = log_distance(p(4),x)-log(p(3))
        if (right <= left) then
            v = wc_small(right,p(1))
        else
            v = ig_complement(wc_small(left,p(1)))
        end if
    end function wc_sf

    pure function wc_interval(a, b, p) result(v)
        implicit none
        real(dp), intent(in) :: a !! finite lower endpoint
        real(dp), intent(in) :: b !! finite upper endpoint greater than a
        real(dp), intent(in) :: p(:) !! concentration, location, scale, upper support
        real(dp) :: v, lo, hi, lp, lq
        v = negative_infinity(0.0_dp)
        lo = max(a,p(2))
        hi = min(b,p(4))
        if (hi <= lo) return
        if (lo == p(2)) then
            v = wc_cdf(hi,p)
        else if (hi == p(4)) then
            v = wc_sf(lo,p)
        else if (log_distance(hi,lo)-log(p(3))+log(2.0_dp)-log1p_safe(-p(1)) &
                < log(0.001_dp)) then
            v = tt_small_mass(lo,hi,p,2)
        else
            v = tail_interval(lo,hi,p,wc_cdf,wc_sf)
            if (v > -scifort_log_two) then
                lp = wc_cdf(lo,p)
                lq = wc_sf(hi,p)
                v = ig_complement(bounded_logadd(lp,lq))
            end if
        end if
    end function wc_interval

    pure function rb_upper(lz, p, lr) result(v)
        implicit none
        real(dp), intent(in) :: lz !! log positive lower standardized coordinate
        real(dp), intent(in) :: p(:) !! shape, location, scale, log density normalization
        real(dp), intent(in), optional :: lr !! log of log(upper/lower) for interval masses
        real(dp) :: v, aa, bb, previous, current, next, total, factor, base
        integer :: k
        aa = 2.0_dp*exp(2.0_dp*(log(p(1))-lz))
        bb = exp(4.0_dp*(log(p(1))-lz))+exp(2.0_dp*log(p(1))-4.0_dp*lz)
        total = 1.0_dp
        previous = 0.0_dp
        current = 1.0_dp
        base = 0.0_dp
        if (present(lr)) base = log_exponential_mass(log(3.0_dp)+lr)
        do k = 1, 64
            next = aa*current-bb*previous
            factor = 3.0_dp/real(2*k+3,dp)
            if (present(lr)) factor = factor*exp( &
                log_exponential_mass(log(real(2*k+3,dp))+lr)-base)
            total = total+factor*next
            previous = current
            current = next
        end do
        v = p(4)+2.0_dp*log(p(1))-3.0_dp*lz-log(3.0_dp)+log(total)+base
    end function rb_upper

    pure function rb_log_density(lz, p) result(v)
        implicit none
        real(dp), intent(in) :: lz !! log positive standardized coordinate, or negative infinity
        real(dp), intent(in) :: p(:) !! shape, location, scale, log normalization
        real(dp) :: v, q, lq
        if (lz > log(2.0_dp*p(1))) then
            lq = 2.0_dp*lz-log(p(1))+log1p_safe(-exp(2.0_dp*(log(p(1))-lz)))
        else
            q = (exp(lz)-p(1))*(exp(lz)+p(1))/p(1)
            lq = bounded_log(abs(q))
        end if
        v = p(4)-log1pexp(2.0_dp*lq)-log(p(3))
    end function rb_log_density

    pure function rb_pdf(x, p) result(v)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: p(:) !! shape, location, scale, log normalization
        real(dp) :: v, lz
        v = negative_infinity(0.0_dp)
        if (x < p(2)) return
        lz = negative_infinity(0.0_dp)
        if (x > p(2)) lz = log_distance(x,p(2))-log(p(3))
        v = rb_log_density(lz,p)
    end function rb_pdf

    pure subroutine rb_tails(x, p, lp, lq)
        implicit none
        real(dp), intent(in) :: x !! finite censoring limit
        real(dp), intent(in) :: p(:) !! shape, location, scale, log normalization
        real(dp), intent(out) :: lp, lq !! log lower and upper probabilities
        real(dp) :: lz, z, lower, upper
        lp = negative_infinity(0.0_dp)
        lq = 0.0_dp
        if (x <= p(2)) return
        lz = log_distance(x,p(2))-log(p(3))
        if (lz < log(1.0e-10_dp*min(p(1),1.0_dp))) then
            lp = p(4)-log1p_safe(p(1)**2)+lz
            lq = ig_complement(lp)
        else if (lz >= log(2.0_dp*max(p(1),sqrt(p(1))))) then
            lq = rb_upper(lz,p)
            lp = ig_complement(lq)
        else
            z = exp(lz)
            lower = rel_breitwigner_cdf(z,p(1))
            upper = rel_breitwigner_sf(z,p(1))
            if (lower < upper) then
                lp = log(lower)
                lq = ig_complement(lp)
            else
                lq = log(upper)
                lp = ig_complement(lq)
            end if
        end if
    end subroutine rb_tails

    pure function rb_cdf(x, p) result(v)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: p(:) !! shape, location, scale, log normalization
        real(dp) :: v, lq
        call rb_tails(x,p,v,lq)
    end function rb_cdf

    pure function rb_sf(x, p) result(v)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: p(:) !! shape, location, scale, log normalization
        real(dp) :: v, lp
        call rb_tails(x,p,lp,v)
    end function rb_sf

    pure function rb_interval(a, b, p) result(v)
        implicit none
        real(dp), intent(in) :: a !! finite lower endpoint
        real(dp), intent(in) :: b !! finite upper endpoint greater than a
        real(dp), intent(in) :: p(:) !! shape, location, scale, log normalization
        real(dp) :: v, la, lb, lw, lr, lp, lq
        v = negative_infinity(0.0_dp)
        if (b <= p(2)) return
        if (a <= p(2)) then
            v = rb_cdf(b,p)
            return
        end if
        la = log_distance(a,p(2))-log(p(3))
        lb = log_distance(b,p(2))-log(p(3))
        lw = log_distance(b,a)
        lr = burr_logsoftplus(lw-log_distance(a,p(2)))
        if (la >= log(2.0_dp*max(p(1),sqrt(p(1))))) then
            v = rb_upper(la,p,lr)
        else if (lb < log(1.0e-10_dp*min(p(1),1.0_dp))) then
            v = p(4)-log1p_safe(p(1)**2)+lw-log(p(3))
        else if (lw-log(p(3))+scifort_log_two+max(0.0_dp,lb)-log(p(1)) < log(0.001_dp)) then
            v = tt_small_mass(a,b,p,3)
        else
            v = tail_interval(a,b,p,rb_cdf,rb_sf)
            if (v > -scifort_log_two) then
                lp = rb_cdf(a,p)
                lq = rb_sf(b,p)
                v = ig_complement(bounded_logadd(lp,lq))
            end if
        end if
    end function rb_interval

    pure function tt_small_mass(a, b, p, family) result(v)
        implicit none
        real(dp), intent(in) :: a, b !! finite ordered interval endpoints inside the support
        real(dp), intent(in) :: p(:) !! family-specific parameters
        integer, intent(in) :: family !! Kolmogorov=1, wrapped Cauchy=2, Breit-Wigner=3
        real(dp) :: v, t(3), term, weight, lw, dl, dr, ls, offset, lp, lq, lz
        integer :: k
        t = [0.5_dp,0.5_dp-sqrt(3.0_dp/5.0_dp)/2.0_dp,0.5_dp+sqrt(3.0_dp/5.0_dp)/2.0_dp]
        lw = log_distance(b,a)
        if (family == 1) then
            dl = log_distance(a,p(1))
            ls = log(p(2))
        else
            dl = log_distance(a,p(2))
            ls = log(p(3))
        end if
        dr = 0.0_dp
        if (family == 2) dr = log_distance(p(4),a)
        v = negative_infinity(0.0_dp)
        do k = 1, 3
            ! Add increments to distances, not to potentially enormous physical locations.
            offset = lw+log(t(k))
            term = bounded_logadd(dl,offset)-ls
            weight = 5.0_dp/18.0_dp
            if (k == 1) weight = 4.0_dp/9.0_dp
            select case (family)
            case (1)
                lz = term
                call kt_logs(lz,lp,lq,term)
                term = term-ls
            case (2)
                term = min(term,log_difference(dr,offset)-ls)
                term = wc_log_density(term,p)
            case (3)
                term = rb_log_density(term,p)
            end select
            v = bounded_logadd(v,term+log(weight))
        end do
        v = v+lw
    end function tt_small_mass

    pure function kstwobign_censored_loglikelihood(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, p(2)
        value = quiet_nan(0.0_dp)
        p = [0.0_dp,1.0_dp]
        if (present(loc)) p(1) = loc
        if (present(scale)) p(2) = scale
        if (.not. all(ieee_is_finite(p))) return
        if (p(2) <= 0.0_dp) return
        value = censored_loglikelihood(data,p,kt_pdf,kt_cdf,kt_sf,kt_interval)
    end function kstwobign_censored_loglikelihood

    pure function kstwobign_censored_nnlf(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -kstwobign_censored_loglikelihood(data, loc, scale)
    end function kstwobign_censored_nnlf

    pure function wrapcauchy_censored_loglikelihood(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! concentration in (0,0.99]
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, p(4)
        value = quiet_nan(0.0_dp)
        p = [c,0.0_dp,1.0_dp,0.0_dp]
        if (present(loc)) p(2) = loc
        if (present(scale)) p(3) = scale
        if (.not. all(ieee_is_finite(p))) return
        if (p(3) <= 0.0_dp) return
        if (c <= 0.0_dp .or. c > 0.99_dp) return
        if (p(3) < tiny(1.0_dp)) return
        if (abs(p(2))/max(1.0_dp,p(3)) > 10.0_dp) return
        if (p(3) < 1.0_dp) then
            if (abs(p(2)) > 10.0_dp*p(3)) return
        end if
        if (p(3) > huge(1.0_dp)/(4.0_dp*scifort_pi)) return
        if (p(2) > huge(1.0_dp)-2.0_dp*scifort_pi*p(3)) return
        p(4) = p(2)+2.0_dp*scifort_pi*p(3)
        if (p(4) <= p(2)) return
        value = censored_loglikelihood(data,p,wc_pdf,wc_cdf,wc_sf,wc_interval)
    end function wrapcauchy_censored_loglikelihood

    pure function wrapcauchy_censored_nnlf(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! concentration in (0,0.99]
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -wrapcauchy_censored_loglikelihood(data, c, loc, scale)
    end function wrapcauchy_censored_nnlf

    pure function rel_breitwigner_censored_loglikelihood(data, rho, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: rho !! resonance-to-width ratio in [0.1,100]
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, p(4), u
        value = quiet_nan(0.0_dp)
        p = [rho,0.0_dp,1.0_dp,0.0_dp]
        if (present(loc)) p(2) = loc
        if (present(scale)) p(3) = scale
        if (.not. all(ieee_is_finite(p))) return
        if (p(3) <= 0.0_dp) return
        if (rho < 0.1_dp .or. rho > 100.0_dp) return
        u = sqrt(1.0_dp+1.0_dp/(rho*rho))
        p(4) = log(2.0_dp/scifort_pi)+0.5_dp*log(2.0_dp*u*u/(1.0_dp+u))
        value = censored_loglikelihood(data,p,rb_pdf,rb_cdf,rb_sf,rb_interval)
    end function rel_breitwigner_censored_loglikelihood

    pure function rel_breitwigner_censored_nnlf(data, rho, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: rho !! resonance-to-width ratio in [0.1,100]
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -rel_breitwigner_censored_loglikelihood(data, rho, loc, scale)
    end function rel_breitwigner_censored_nnlf


    ! Positive Poisson mixtures derived from the specifications cited in provenance.
    pure function pm_component(x, y, p, j, mode) result(v)
        implicit none
        real(dp), intent(in) :: x !! observation or lower endpoint
        real(dp), intent(in) :: y !! upper endpoint for interval mode
        real(dp), intent(in) :: p(:) !! family, location, scale, gamma/beta shapes, Poisson mean, ratio
        integer, intent(in) :: j !! nonnegative Poisson index
        integer, intent(in) :: mode !! zero density, one CDF, two SF, three interval
        real(dp) :: v, a, ly, ld, jac, bp(4), lo, lt, lh, threshold
        a = p(4)+real(j,dp)
        if (p(1) == 3.0_dp) then
            bp = [a,p(5),p(2),p(3)*p(7)]
            lt = log_distance(x,p(2))-log(bp(4))
            lh = lt
            if (mode == 3) lh = log_distance(y,p(2))-log(bp(4))
            threshold = log(epsilon(1.0_dp)/16.0_dp)-log(a+p(5))
            if (lh < threshold) then
                ! Integrate the endpoint power before converting subnormal coordinates.
                v = beta_small_tail(a,p(5),lt)
                select case (mode)
                case (0)
                    v = v+log(a)-log_distance(x,p(2))
                case (2)
                    v = ig_complement(v)
                case (3)
                    ld = burr_logsoftplus(log_distance(y,x)-log_distance(x,p(2)))
                    v = beta_small_tail(a,p(5),lh)+log_exponential_mass(log(a)+ld)
                end select
                return
            else if (lt > -threshold) then
                v = beta_small_tail(p(5),a,-lt)
                select case (mode)
                case (0)
                    v = v+log(p(5))-log_distance(x,p(2))
                case (1)
                    v = ig_complement(v)
                case (3)
                    ld = burr_logsoftplus(log_distance(y,x)-log_distance(x,p(2)))
                    v = v+log_exponential_mass(log(p(5))+ld)
                end select
                return
            end if
            select case (mode)
            case (0)
                v = betaprime_lpdf(x,bp)
            case (1)
                v = betaprime_lcdf(x,bp)
            case (2)
                v = betaprime_lsf(x,bp)
            case (3)
                v = betaprime_interval(x,y,bp)
            end select
            return
        end if
        ly = negative_infinity(0.0_dp)
        jac = 0.0_dp
        if (x > p(2)) then
            ld = log_distance(x,p(2))-log(p(3))
            if (p(1) == 1.0_dp) then
                ly = ld-scifort_log_two
                jac = -log(p(3))-scifort_log_two
            else
                ly = 2.0_dp*ld-scifort_log_two
                jac = ld-log(p(3))
            end if
        end if
        select case (mode)
        case (0)
            v = loggamma_lpdf(ly,[a,0.0_dp,1.0_dp])-ly+jac
        case (1)
            v = loggamma_lcdf(ly,[a,0.0_dp,1.0_dp])
        case (2)
            v = loggamma_lsf(ly,[a,0.0_dp,1.0_dp])
        case (3)
            lo = max(x,p(2))
            ld = log_distance(y,lo)-log(p(3))-scifort_log_two
            if (p(1) == 2.0_dp) ld = ld+bounded_logadd( &
                log_distance(y,p(2)),log_distance(lo,p(2)))-log(p(3))
            v = gk_gamma_mass(ly,ld,a)
        end select
    end function pm_component

    pure function pm_sum(x, y, p, mode) result(v)
        implicit none
        real(dp), intent(in) :: x !! observation or lower endpoint
        real(dp), intent(in) :: y !! upper endpoint for interval mode
        real(dp), intent(in) :: p(:) !! family, location, scale, shapes, Poisson mean, ratio
        integer, intent(in) :: mode !! zero density, one CDF, two SF, three interval
        real(dp) :: v, lw, term, lr, a, endpoint, ly, residual
        integer :: j
        v = negative_infinity(0.0_dp)
        lw = -p(6)
        do j = 0, 10000
            term = lw+pm_component(x,y,p,j,mode)
            if (ieee_is_nan(term)) then
                v = term
                return
            end if
            v = bounded_logadd(v,term)
            if (p(6) == 0.0_dp) return
            a = p(4)+real(j,dp)
            lr = log(p(6))-log(real(j+1,dp))
            if (p(1) == 3.0_dp) then
                ! The beta density ratio bounds every partial integral as well.
                lr = lr+log1p_safe(p(5)/a)
            else if (mode /= 1) then
                endpoint = x
                if (mode == 3) endpoint = y
                ly = log_distance(endpoint,p(2))-log(p(3))-scifort_log_two
                if (p(1) == 2.0_dp) ly = 2.0_dp*(ly+scifort_log_two)-scifort_log_two
                if (mode == 2) then
                    ! For a>=1, Q(a+1,t)/Q(a,t) <= 1+t/a.
                    if (a < 1.0_dp) then
                        lw = lw+log(p(6))-log(real(j+1,dp))
                        cycle
                    end if
                    lr = lr+log1pexp(ly-log(a))
                else
                    lr = lr+ly-log(a)
                end if
            end if
            if (lr < 0.0_dp) then
                residual = term+lr-ig_complement(lr)
                if (residual <= v+log(2.0e-16_dp)) return
                if (.not. ieee_is_finite(v) .and. term == v) return
            end if
            lw = lw+log(p(6))-log(real(j+1,dp))
        end do
        ! Do not report a truncated series as a valid tail probability.
        v = quiet_nan(0.0_dp)
    end function pm_sum

    pure function pm_lpdf(x, p) result(v)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: p(:) !! prepared mixture parameters
        real(dp) :: v
        v = negative_infinity(0.0_dp)
        if (x < p(2)) return
        if (x == p(2)) then
            if (p(1) == 2.0_dp) return
            if (p(4) < 1.0_dp) v = positive_infinity(0.0_dp)
            if (p(4) == 1.0_dp) then
                v = -p(6)-log(p(3))-scifort_log_two
                if (p(1) == 3.0_dp) v = -p(6)+log(p(5))-log(p(3))-log(p(7))
            end if
            return
        end if
        v = pm_sum(x,x,p,0)
    end function pm_lpdf

    pure function pm_tail(x, p, upper) result(v)
        implicit none
        real(dp), intent(in) :: x !! finite censoring limit
        real(dp), intent(in) :: p(:) !! prepared mixture parameters
        logical, intent(in) :: upper !! select survival probability
        real(dp) :: v, opposite
        integer :: mode
        v = negative_infinity(0.0_dp)
        if (upper) v = 0.0_dp
        if (x <= p(2)) return
        mode = 1
        if (upper) mode = 2
        v = pm_sum(x,x,p,mode)
        if (ieee_is_nan(v)) return
        if (v > -scifort_log_two) then
            opposite = pm_sum(x,x,p,3-mode)
            v = opposite
            if (ieee_is_nan(v)) return
            v = ig_complement(opposite)
        end if
    end function pm_tail

    pure function pm_lcdf(x, p) result(v)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: p(:) !! prepared mixture parameters
        real(dp) :: v
        v = pm_tail(x,p,.false.)
    end function pm_lcdf

    pure function pm_lsf(x, p) result(v)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: p(:) !! prepared mixture parameters
        real(dp) :: v
        v = pm_tail(x,p,.true.)
    end function pm_lsf

    pure function pm_interval(x, y, p) result(v)
        implicit none
        real(dp), intent(in) :: x !! finite lower interval endpoint
        real(dp), intent(in) :: y !! finite upper endpoint greater than x
        real(dp), intent(in) :: p(:) !! prepared mixture parameters
        real(dp) :: v, lp, lq
        v = negative_infinity(0.0_dp)
        if (y <= p(2)) return
        if (x <= p(2)) then
            v = pm_lcdf(y,p)
        else
            v = pm_sum(x,y,p,3)
            if (ieee_is_nan(v)) return
            if (v > -scifort_log_two) then
                lp = pm_lcdf(x,p)
                lq = pm_lsf(y,p)
                if (ieee_is_nan(lp) .or. ieee_is_nan(lq)) then
                    v = quiet_nan(0.0_dp)
                    return
                end if
                v = ig_complement(bounded_logadd(lp,lq))
            end if
        end if
    end function pm_interval

    pure function ncx2_censored_loglikelihood(data, df, nc, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: df !! degrees of freedom in [0.1,2000]
        real(dp), intent(in) :: nc !! noncentrality in [0,100]
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, p(7)
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite([df, nc]))) return
        if (df < 0.1_dp .or. df > 2000.0_dp) return
        if (nc < 0.0_dp .or. nc > 100.0_dp) return
        p = [1.0_dp,0.0_dp,1.0_dp,df/2.0_dp,0.0_dp,nc/2.0_dp,1.0_dp]
        if (present(loc)) p(2) = loc
        if (present(scale)) p(3) = scale
        if (.not. all(ieee_is_finite(p))) return
        if (p(3) <= 0.0_dp) return
        value = censored_loglikelihood(data,p,pm_lpdf,pm_lcdf,pm_lsf,pm_interval)
    end function ncx2_censored_loglikelihood

    pure function ncx2_censored_nnlf(data, df, nc, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: df !! degrees of freedom in [0.1,2000]
        real(dp), intent(in) :: nc !! noncentrality in [0,100]
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -ncx2_censored_loglikelihood(data,df, nc,loc,scale)
    end function ncx2_censored_nnlf

    pure function rice_censored_loglikelihood(data, b, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: b !! Rice shape in [0,10]
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, p(7)
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite([b]))) return
        if (b < 0.0_dp .or. b > 10.0_dp) return
        p = [2.0_dp,0.0_dp,1.0_dp,1.0_dp,0.0_dp,b*b/2.0_dp,1.0_dp]
        if (present(loc)) p(2) = loc
        if (present(scale)) p(3) = scale
        if (.not. all(ieee_is_finite(p))) return
        if (p(3) <= 0.0_dp) return
        value = censored_loglikelihood(data,p,pm_lpdf,pm_lcdf,pm_lsf,pm_interval)
    end function rice_censored_loglikelihood

    pure function rice_censored_nnlf(data, b, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: b !! Rice shape in [0,10]
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -rice_censored_loglikelihood(data,b,loc,scale)
    end function rice_censored_nnlf

    pure function ncf_censored_loglikelihood(data, dfn, dfd, nc, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: dfn !! degrees of freedom in [0.1,2000]
        real(dp), intent(in) :: dfd !! degrees of freedom in [0.1,2000]
        real(dp), intent(in) :: nc !! noncentrality in [0,100]
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, p(7)
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite([dfn, dfd, nc]))) return
        if (dfn < 0.1_dp .or. dfn > 2000.0_dp) return
        if (dfd < 0.1_dp .or. dfd > 2000.0_dp) return
        if (nc < 0.0_dp .or. nc > 100.0_dp) return
        p = [3.0_dp,0.0_dp,1.0_dp,dfn/2.0_dp,dfd/2.0_dp,nc/2.0_dp,dfd/dfn]
        if (present(loc)) p(2) = loc
        if (present(scale)) p(3) = scale
        if (.not. all(ieee_is_finite(p))) return
        if (p(3) <= 0.0_dp) return
        if (p(7) > 1.0_dp) then
            if (p(3) > huge(1.0_dp)/p(7)) return
        end if
        if (p(3) < tiny(1.0_dp)/p(7)) return
        value = censored_loglikelihood(data,p,pm_lpdf,pm_lcdf,pm_lsf,pm_interval)
    end function ncf_censored_loglikelihood

    pure function ncf_censored_nnlf(data, dfn, dfd, nc, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: dfn !! degrees of freedom in [0.1,2000]
        real(dp), intent(in) :: dfd !! degrees of freedom in [0.1,2000]
        real(dp), intent(in) :: nc !! noncentrality in [0,100]
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -ncf_censored_loglikelihood(data,dfn, dfd, nc,loc,scale)
    end function ncf_censored_nnlf


    ! Parameters: family, loc, scale, up to four shapes, GL16 nodes and weights.
    pure subroutine nm_prepare(par, valid)
        implicit none
        real(dp), intent(inout) :: par(39) !! packed family parameters and local quadrature rule
        logical, intent(out) :: valid !! true for supported finite parameters
        real(dp), allocatable :: nodes(:), weights(:)
        integer :: status
        valid = .false.
        if (.not. all(ieee_is_finite(par))) return
        if (par(3) <= 0.0_dp) return
        select case (int(par(1)))
        case (1)
            if (par(4) < 0.01_dp .or. par(4) > 100.0_dp) return
        case (2)
            if (abs(par(4)) > 10.0_dp) return
            call roots_legendre(16,nodes,weights,status=status)
            if (status /= 0) return
            par(8:23) = nodes
            par(24:39) = weights
        case (3)
            if (abs(par(4)) > 100.0_dp) return
            if (par(5) < 0.01_dp .or. par(5) > 10.0_dp) return
            if (minval(par(6:7)) < 0.01_dp .or. maxval(par(6:7)) > 100.0_dp) return
        end select
        valid = .true.
    end subroutine nm_prepare

    pure function nm_normal_pdf(z) result(value)
        implicit none
        real(dp), intent(in) :: z !! real standard coordinate, possibly infinite
        real(dp) :: value
        value = negative_infinity(z)
        if (abs(z) > sqrt(2.0_dp)*sqrt(huge(1.0_dp))) return
        value = normal_logpdf(z)
    end function nm_normal_pdf

    pure function nm_emg_h(z, rate) result(value)
        implicit none
        real(dp), intent(in) :: z !! standardized real normal coordinate
        real(dp), intent(in) :: rate !! exponential rate in [1e-4,1000]
        real(dp) :: value, q
        value = negative_infinity(z)
        if (.not. ieee_is_finite(z)) return
        if (z <= rate) then
            q = rate-z
            if (q == 0.0_dp) then
                value = nm_normal_pdf(z)+0.5_dp*log(scifort_pi/2.0_dp)
            else
                value = nm_normal_pdf(z)+ig_logmills(log(q))
            end if
        else
            if (rate > 1.0_dp) then
                if (z > huge(1.0_dp)/rate) return
            end if
            value = rate*(rate/2.0_dp-z)+normal_logcdf(z-rate)
        end if
    end function nm_emg_h

    pure function nm_emg_pdf(z, rate) result(value)
        implicit none
        real(dp), intent(in) :: z !! standardized real normal coordinate
        real(dp), intent(in) :: rate !! positive exponential rate
        real(dp) :: value
        value = log(rate)+nm_emg_h(z,rate)
    end function nm_emg_pdf

    pure function nm_emg_cdf(z, rate) result(value)
        implicit none
        real(dp), intent(in) :: z !! standardized real normal coordinate
        real(dp), intent(in) :: rate !! positive exponential rate
        real(dp) :: value, upper
        if (.not. ieee_is_finite(z)) then
            value = 0.0_dp
            if (z < 0.0_dp) value = negative_infinity(z)
            return
        end if
        if (z <= -12.0_dp) then
            value = nm_normal_pdf(z)+ig_mills_difference(log(-z),log(rate-z),log(rate))
        else if (z > 0.0_dp) then
            upper = min(0.0_dp,bounded_logadd(normal_logsf(z),nm_emg_h(z,rate)))
            if (upper < -scifort_log_two) then
                value = ig_complement(upper)
            else
                value = log_difference(normal_logcdf(z),nm_emg_h(z,rate))
            end if
        else
            value = log_difference(normal_logcdf(z),nm_emg_h(z,rate))
        end if
    end function nm_emg_cdf

    pure function nm_emg_sf(z, rate) result(value)
        implicit none
        real(dp), intent(in) :: z !! standardized real normal coordinate
        real(dp), intent(in) :: rate !! positive exponential rate
        real(dp) :: value
        if (z <= 0.0_dp) then
            value = ig_complement(nm_emg_cdf(z,rate))
        else
            value = min(0.0_dp,bounded_logadd(normal_logsf(z),nm_emg_h(z,rate)))
            if (value > -scifort_log_two) value = ig_complement(nm_emg_cdf(z,rate))
        end if
    end function nm_emg_sf

    pure function nm_emg_mass(za, zb, lw, rate) result(value)
        implicit none
        real(dp), intent(in) :: za, zb !! ordered standard coordinates, possibly rounded equal
        real(dp), intent(in) :: lw !! log of positive standard interval width
        real(dp), intent(in) :: rate !! positive exponential rate
        real(dp) :: value, lbound, width, lp, lq, node, total
        integer :: j
        if (za > rate+40.0_dp) then
            value = nm_emg_h(za,rate)+log_exponential_mass(log(rate)+lw)
            return
        end if
        if (ieee_is_finite(za) .and. ieee_is_finite(zb)) then
            lbound = log(max(1.0_dp,abs(za),abs(zb)))+log(2.0_dp+rate)
            if (lw+lbound < log(0.001_dp)) then
                width = exp(lw)
                total = negative_infinity(0.0_dp)
                do j = -1, 1
                    node = (1.0_dp+real(j,dp)*sqrt(3.0_dp/5.0_dp))/2.0_dp
                    lp = log(5.0_dp/18.0_dp)
                    if (j == 0) lp = log(4.0_dp/9.0_dp)
                    total = bounded_logadd(total,lp+nm_emg_pdf(za+node*width,rate))
                end do
                value = lw+total
                return
            end if
        end if
        lp = nm_emg_cdf(zb,rate)
        lq = nm_emg_sf(za,rate)
        if (lp <= lq) then
            value = log_difference(lp,nm_emg_cdf(za,rate))
        else
            value = log_difference(lq,nm_emg_sf(zb,rate))
        end if
    end function nm_emg_mass

    pure function nm_skew_pdf(z, shape) result(value)
        implicit none
        real(dp), intent(in) :: z !! real standardized observation
        real(dp), intent(in) :: shape !! skewness in [-10,10]
        real(dp) :: value, lp, lq, az
        value = negative_infinity(z)
        if (.not. ieee_is_finite(z)) return
        lp = nm_normal_pdf(z)
        if (.not. ieee_is_finite(lp)) return
        if (abs(shape) > 1.0_dp) then
            if (abs(z) > huge(1.0_dp)/abs(shape)) return
        end if
        az = shape*z
        lq = normal_logcdf(az)
        if (.not. ieee_is_finite(lq)) return
        if (lp < -huge(1.0_dp)-lq) return
        value = scifort_log_two+lp+lq
    end function nm_skew_pdf

    pure function nm_skew_positive(z, shape, par) result(value)
        implicit none
        real(dp), intent(in) :: z !! nonpositive standardized threshold
        real(dp), intent(in) :: shape !! nonnegative skewness in [0,10]
        real(dp), intent(in) :: par(:) !! local GL16 rule in elements 8:39
        real(dp) :: value, base, y, ay, rate, v, w, relative, total, logmills
        integer :: i, j
        base = nm_skew_pdf(z,shape)
        value = base
        if (.not. ieee_is_finite(base)) return
        y = -z
        ay = shape*y
        rate = 1.0_dp+shape+(1.0_dp+shape*shape)*y
        logmills = 0.0_dp
        if (ay >= 12.0_dp) logmills = ig_logmills(log(ay))
        total = negative_infinity(0.0_dp)
        do i = 0, 23
            do j = 1, 16
                v = real(2*i+1,dp)+par(7+j)
                w = v/rate
                if (ay >= 12.0_dp) then
                    relative = -(1.0_dp+shape*shape)*y*w &
                        -0.5_dp*(1.0_dp+shape*shape)*w*w &
                        +ig_logmills(log(ay+shape*w))-logmills
                else
                    relative = -y*w-0.5_dp*w*w &
                        +normal_logcdf(-ay-shape*w)-normal_logcdf(-ay)
                end if
                total = bounded_logadd(total,log(par(23+j))+relative)
            end do
        end do
        value = min(0.0_dp,base-log(rate)+total)
    end function nm_skew_positive

    pure function nm_skew_lower(z, shape, par) result(value)
        implicit none
        real(dp), intent(in) :: z !! nonpositive standard coordinate
        real(dp), intent(in) :: shape !! skewness in [-10,10]
        real(dp), intent(in) :: par(:) !! prepared family parameters
        real(dp) :: value
        if (shape == 0.0_dp) then
            value = normal_logcdf(z)
        else if (shape > 0.0_dp) then
            value = nm_skew_positive(z,shape,par)
        else
            value = log_difference(scifort_log_two+normal_logcdf(z), &
                nm_skew_positive(z,-shape,par))
        end if
    end function nm_skew_lower

    pure function nm_skew_tail(z, shape, par) result(value)
        implicit none
        real(dp), intent(in) :: z !! real standard coordinate
        real(dp), intent(in) :: shape !! skewness in [-10,10]
        real(dp), intent(in) :: par(:) !! prepared family parameters
        real(dp) :: value
        if (z <= 0.0_dp) then
            value = nm_skew_lower(z,shape,par)
        else
            value = ig_complement(nm_skew_lower(-z,-shape,par))
        end if
    end function nm_skew_tail

    pure function nm_dpl_tail(z, par, upper) result(value)
        implicit none
        real(dp), intent(in) :: z !! standardized lognormal coordinate
        real(dp), intent(in) :: par(:) !! u, s, a, b in elements 4:7
        logical, intent(in) :: upper !! true for log survival probability
        real(dp) :: value, plus, minus, wa, wb
        wa = log(par(7))-log(par(6)+par(7))
        wb = log(par(6))-log(par(6)+par(7))
        if (upper) then
            plus = nm_emg_sf(z,par(6)*par(5))
            minus = nm_emg_cdf(-z,par(7)*par(5))
        else
            plus = nm_emg_cdf(z,par(6)*par(5))
            minus = nm_emg_sf(-z,par(7)*par(5))
        end if
        value = min(0.0_dp,bounded_logadd(wa+plus,wb+minus))
        if (value > -scifort_log_two) then
            if (upper) then
                plus = nm_emg_cdf(z,par(6)*par(5))
                minus = nm_emg_sf(-z,par(7)*par(5))
            else
                plus = nm_emg_sf(z,par(6)*par(5))
                minus = nm_emg_cdf(-z,par(7)*par(5))
            end if
            value = ig_complement(min(0.0_dp,bounded_logadd(wa+plus,wb+minus)))
        end if
    end function nm_dpl_tail

    pure function nm_lpdf(x, par) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: par(:) !! prepared family parameters
        real(dp) :: value, z, ld, plus, minus, norm
        if (int(par(1)) == 3) then
            value = negative_infinity(x)
            if (x <= par(2)) return
            ld = log_distance(x,par(2))
            z = (ld-log(par(3))-par(4))/par(5)
            norm = log(par(6)+par(7))
            plus = log(par(7))-norm+nm_emg_pdf(z,par(6)*par(5))
            minus = log(par(6))-norm+nm_emg_pdf(-z,par(7)*par(5))
            value = bounded_logadd(plus,minus)-log(par(5))-ld
        else
            z = location_scale_z(x,par(2:3))
            if (int(par(1)) == 1) then
                value = nm_emg_pdf(z,1.0_dp/par(4))-log(par(3))
            else
                value = nm_skew_pdf(z,par(4))-log(par(3))
            end if
        end if
    end function nm_lpdf

    pure function nm_lcdf(x, par) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite left-censoring threshold
        real(dp), intent(in) :: par(:) !! prepared family parameters
        real(dp) :: value, z
        if (int(par(1)) == 3) then
            value = negative_infinity(x)
            if (x <= par(2)) return
            z = (log_distance(x,par(2))-log(par(3))-par(4))/par(5)
            value = nm_dpl_tail(z,par,.false.)
        else
            z = location_scale_z(x,par(2:3))
            if (int(par(1)) == 1) then
                value = nm_emg_cdf(z,1.0_dp/par(4))
            else
                value = nm_skew_tail(z,par(4),par)
            end if
        end if
    end function nm_lcdf

    pure function nm_lsf(x, par) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite right-censoring threshold
        real(dp), intent(in) :: par(:) !! prepared family parameters
        real(dp) :: value, z
        if (int(par(1)) == 3) then
            value = 0.0_dp
            if (x <= par(2)) return
            z = (log_distance(x,par(2))-log(par(3))-par(4))/par(5)
            value = nm_dpl_tail(z,par,.true.)
        else
            z = location_scale_z(x,par(2:3))
            if (int(par(1)) == 1) then
                value = nm_emg_sf(z,1.0_dp/par(4))
            else
                value = nm_skew_tail(-z,-par(4),par)
            end if
        end if
    end function nm_lsf

    pure function nm_interval(a, b, par) result(value)
        implicit none
        real(dp), intent(in) :: a, b !! ordered finite physical interval endpoints
        real(dp), intent(in) :: par(:) !! prepared family parameters
        real(dp) :: value, za, zb, lw, plus, minus, norm, ld, width, total, node, weight
        integer :: j
        if (int(par(1)) == 3) then
            value = negative_infinity(a)
            if (b <= par(2)) return
            if (a <= par(2)) then
                value = nm_lcdf(b,par)
                return
            end if
            ld = log_distance(a,par(2))
            za = (ld-log(par(3))-par(4))/par(5)
            zb = (log_distance(b,par(2))-log(par(3))-par(4))/par(5)
            lw = burr_logsoftplus(log_distance(b,a)-ld)-log(par(5))
            norm = log(par(6)+par(7))
            plus = log(par(7))-norm+nm_emg_mass(za,zb,lw,par(6)*par(5))
            minus = log(par(6))-norm+nm_emg_mass(-zb,-za,lw,par(7)*par(5))
            value = bounded_logadd(plus,minus)
            return
        end if
        za = location_scale_z(a,par(2:3))
        zb = location_scale_z(b,par(2:3))
        lw = log_distance(b,a)-log(par(3))
        if (int(par(1)) == 1) then
            value = nm_emg_mass(za,zb,lw,1.0_dp/par(4))
            return
        end if
        if (ieee_is_finite(za) .and. ieee_is_finite(zb)) then
            ld = log(max(1.0_dp,abs(za),abs(zb)))+log(2.0_dp+abs(par(4))+par(4)**2)
            if (lw+ld < log(0.001_dp)) then
                width = exp(lw)
                total = negative_infinity(0.0_dp)
                do j = -1, 1
                    node = (1.0_dp+real(j,dp)*sqrt(3.0_dp/5.0_dp))/2.0_dp
                    weight = log(5.0_dp/18.0_dp)
                    if (j == 0) weight = log(4.0_dp/9.0_dp)
                    total = bounded_logadd(total,weight+nm_skew_pdf(za+node*width,par(4)))
                end do
                value = lw+total
                return
            end if
        end if
        value = tail_interval(a,b,par,nm_lcdf,nm_lsf)
    end function nm_interval

    pure function exponnorm_censored_loglikelihood(data, k, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: k !! normal-exponential scale ratio in [0.01,100]
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(39)
        logical :: valid
        params = 0.0_dp
        params(1:3) = [1.0_dp,0.0_dp,1.0_dp]
        params(4) = k
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        value = quiet_nan(0.0_dp)
        call nm_prepare(params,valid)
        if (.not. valid) return
        value = censored_loglikelihood(data,params,nm_lpdf,nm_lcdf,nm_lsf,nm_interval)
    end function exponnorm_censored_loglikelihood

    pure function exponnorm_censored_nnlf(data, k, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: k !! normal-exponential scale ratio in [0.01,100]
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -exponnorm_censored_loglikelihood(data,k,loc,scale)
    end function exponnorm_censored_nnlf

    pure function skewnorm_censored_loglikelihood(data, a, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: a !! skewness in [-10,10]
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(39)
        logical :: valid
        params = 0.0_dp
        params(1:3) = [2.0_dp,0.0_dp,1.0_dp]
        params(4) = a
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        value = quiet_nan(0.0_dp)
        call nm_prepare(params,valid)
        if (.not. valid) return
        value = censored_loglikelihood(data,params,nm_lpdf,nm_lcdf,nm_lsf,nm_interval)
    end function skewnorm_censored_loglikelihood

    pure function skewnorm_censored_nnlf(data, a, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: a !! skewness in [-10,10]
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -skewnorm_censored_loglikelihood(data,a,loc,scale)
    end function skewnorm_censored_nnlf

    pure function dpareto_lognorm_censored_loglikelihood(data, u, s, a, b, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: u !! lognormal location shape in [-100,100]
        real(dp), intent(in) :: s !! lognormal scale shape in [0.01,10]
        real(dp), intent(in) :: a !! upper-Pareto shape in [0.01,100]
        real(dp), intent(in) :: b !! lower-Pareto shape in [0.01,100]
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(39)
        logical :: valid
        params = 0.0_dp
        params(1:3) = [3.0_dp,0.0_dp,1.0_dp]
        params(4) = u
        params(5) = s
        params(6) = a
        params(7) = b
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        value = quiet_nan(0.0_dp)
        call nm_prepare(params,valid)
        if (.not. valid) return
        value = censored_loglikelihood(data,params,nm_lpdf,nm_lcdf,nm_lsf,nm_interval)
    end function dpareto_lognorm_censored_loglikelihood

    pure function dpareto_lognorm_censored_nnlf(data, u, s, a, b, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: u !! lognormal location shape in [-100,100]
        real(dp), intent(in) :: s !! lognormal scale shape in [0.01,10]
        real(dp), intent(in) :: a !! upper-Pareto shape in [0.01,100]
        real(dp), intent(in) :: b !! lower-Pareto shape in [0.01,100]
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -dpareto_lognorm_censored_loglikelihood(data,u, s, a, b,loc,scale)
    end function dpareto_lognorm_censored_nnlf


    ! Bounded censoring: family, loc, scale, shapes, normalization, GL16 rule.
    pure subroutine bd_prepare(par, valid)
        implicit none
        real(dp), intent(inout) :: par(40) !! packed distribution parameters and local rule
        logical, intent(out) :: valid !! true when the validated numerical contract holds
        real(dp), allocatable :: nodes(:), weights(:)
        real(dp) :: extent
        integer :: status
        valid = .false.
        if (.not. all(ieee_is_finite(par))) return
        if (par(3) < tiny(1.0_dp)) return
        ! Keep affine endpoints finite and avoid severe location cancellation.
        if (abs(par(2))/max(1.0_dp,par(3)) > 1.0e4_dp) return
        if (par(3) < 1.0_dp) then
            if (abs(par(2)) > 1.0e4_dp*par(3)) return
        end if
        select case (int(par(1)))
        case (1)
            if (minval(par(4:5)) < 0.1_dp .or. maxval(par(4:5)) > 20.0_dp) return
            if (abs(par(6)) > 10.0_dp .or. par(7) < -0.5_dp .or. par(7) > 10.0_dp) return
            extent = 1.0_dp
            par(8) = gausshyper_logkernel_integral(0.0_dp,1.0_dp,par(4),par(5),par(6),par(7))
        case (2)
            if (par(4) < 1.0_dp .or. par(4) > 32.0_dp) return
            if (par(4) /= real(nint(par(4)),dp)) return
            extent = par(4)
        case (3)
            if (par(4) < 0.0_dp .or. par(4) > 100.0_dp) return
            extent = scifort_pi
            par(8) = -vonmises_standard_logpdf(0.0_dp,par(4))
        end select
        if (par(3) > huge(1.0_dp)/(2.0_dp*extent)) return
        if (abs(par(2)) > huge(1.0_dp)-par(3)*extent) return
        call roots_legendre(16,nodes,weights,status=status)
        if (status /= 0) return
        par(9:24) = nodes
        par(25:40) = weights
        valid = all(ieee_is_finite(par))
    end subroutine bd_prepare

    pure subroutine bd_support(par, lo, hi)
        implicit none
        real(dp), intent(in) :: par(:) !! prepared bounded family parameters
        real(dp), intent(out) :: lo, hi !! standardized support endpoints
        lo = 0.0_dp
        hi = 1.0_dp
        if (int(par(1)) == 2) hi = par(4)
        if (int(par(1)) == 3) then
            lo = -scifort_pi
            hi = scifort_pi
        end if
    end subroutine bd_support

    pure function bd_spline(x, order) result(value)
        implicit none
        real(dp), intent(in) :: x !! standard coordinate
        integer, intent(in) :: order !! cardinal B-spline order in [1,33]
        real(dp) :: value, old(0:34), new(0:34), t, neginf
        integer :: k, r
        neginf = negative_infinity(0.0_dp)
        value = neginf
        if (x < 0.0_dp .or. x > real(order,dp)) return
        if (order == 1) then
            value = 0.0_dp
            return
        end if
        if (x == 0.0_dp .or. x == real(order,dp)) return
        ! The first/last piece is a single power, even when its density underflows.
        t = min(x,real(order,dp)-x)
        if (t <= 1.0_dp) then
            value = real(order-1,dp)*log(t)-log_gamma(real(order,dp))
            return
        end if
        old = neginf
        do k = 0, order
            t = x-real(k,dp)
            if (t >= 0.0_dp .and. t < 1.0_dp) old(k) = 0.0_dp
        end do
        do r = 2, order
            new = neginf
            do k = 0, order+1-r
                t = x-real(k,dp)
                if (old(k) /= neginf .and. t > 0.0_dp) then
                    new(k) = log(t)+old(k)-log(real(r-1,dp))
                end if
                if (old(k+1) /= neginf .and. t < real(r,dp)) then
                    new(k) = bounded_logadd(new(k), &
                        log(real(r,dp)-t)+old(k+1)-log(real(r-1,dp)))
                end if
            end do
            old = new
        end do
        value = old(0)
    end function bd_spline

    pure function bd_standard_pdf(x, par) result(value)
        implicit none
        real(dp), intent(in) :: x !! standard coordinate inside the support
        real(dp), intent(in) :: par(:) !! prepared bounded family parameters
        real(dp) :: value
        value = 0.0_dp
        select case (int(par(1)))
        case (1)
            if (x == 0.0_dp .and. par(4) /= 1.0_dp) then
                value = sign(positive_infinity(x),1.0_dp-par(4))
                return
            end if
            if (x == 1.0_dp .and. par(5) /= 1.0_dp) then
                value = sign(positive_infinity(x),1.0_dp-par(5))
                return
            end if
            if (par(4) /= 1.0_dp) value = (par(4)-1.0_dp)*log(x)
            if (par(5) /= 1.0_dp) value = value+(par(5)-1.0_dp)*log1p_safe(-x)
            value = value-par(6)*log1p_safe(par(7)*x)-par(8)
        case (2)
            value = bd_spline(x,int(par(4)))
        case (3)
            ! -2*sin(x/2)**2 avoids cancellation in cos(x)-1 near zero.
            value = -2.0_dp*par(4)*sin(x/2.0_dp)**2-par(8)
        end select
    end function bd_standard_pdf

    pure function bd_lpdf(x, par) result(value)
        implicit none
        real(dp), intent(in) :: x !! physical observation
        real(dp), intent(in) :: par(:) !! prepared bounded family parameters
        real(dp) :: value, z, lo, hi
        z = location_scale_z(x,par(2:3))
        call bd_support(par,lo,hi)
        value = negative_infinity(x)
        if (z < lo .or. z > hi) return
        value = bd_standard_pdf(z,par)-log(par(3))
    end function bd_lpdf

    pure function bd_gauss(lo, logwidth, par) result(value)
        implicit none
        real(dp), intent(in) :: lo !! lower standard coordinate
        real(dp), intent(in) :: logwidth !! logarithm of positive standard width
        real(dp), intent(in) :: par(:) !! prepared bounded family parameters and GL16 rule
        real(dp) :: value, width, x
        integer :: j
        width = exp(logwidth)
        value = negative_infinity(0.0_dp)
        do j = 1, 16
            x = lo+width*((1.0_dp+par(8+j))/2.0_dp)
            value = bounded_logadd(value,log(par(24+j))+bd_standard_pdf(x,par))
        end do
        value = value+logwidth-scifort_log_two
    end function bd_gauss

    pure function bd_mass(lo, hi, par) result(value)
        implicit none
        real(dp), intent(in) :: lo, hi !! ordered standard endpoints clipped to support
        real(dp), intent(in) :: par(:) !! prepared bounded family parameters
        real(dp) :: value, a, b, h, part, left, right, ratio
        integer :: i, panels
        value = negative_infinity(0.0_dp)
        if (lo >= hi) return
        select case (int(par(1)))
        case (1)
            value = gausshyper_logkernel_integral(lo,hi,par(4),par(5),par(6),par(7))-par(8)
        case (2)
            ! A degree-(n-1) polynomial on each unit cell: GL16 is exact for n<=32.
            left = lo
            right = hi
            if (lo > par(4)/2.0_dp) then
                left = par(4)-hi
                right = par(4)-lo
            end if
            if (right <= 1.0_dp) then
                ! Integrate the outer power analytically, including subnormal widths.
                value = par(4)*log(right)-log_gamma(par(4)+1.0_dp)
                if (left > 0.0_dp) then
                    ratio = (right-left)/right
                    if (ratio < 1.0_dp) value = value &
                        +log(-expm1_safe(par(4)*log1p_safe(-ratio)))
                end if
                return
            end if
            a = left
            do while (a < right)
                b = min(right,real(floor(a)+1,dp))
                part = bd_gauss(a,log(b-a),par)
                value = bounded_logadd(value,part)
                a = b
            end do
        case (3)
            panels = max(1,ceiling(4.0_dp*(hi-lo)))
            h = (hi-lo)/real(panels,dp)
            do i = 0, panels-1
                a = lo+real(i,dp)*h
                value = bounded_logadd(value,bd_gauss(a,log(h),par))
            end do
        end select
        value = min(0.0_dp,value)
    end function bd_mass

    pure function bd_lcdf(x, par) result(value)
        implicit none
        real(dp), intent(in) :: x !! physical left-censoring threshold
        real(dp), intent(in) :: par(:) !! prepared bounded family parameters
        real(dp) :: value, z, lo, hi, tail
        z = location_scale_z(x,par(2:3))
        call bd_support(par,lo,hi)
        value = negative_infinity(x)
        if (z <= lo) return
        value = 0.0_dp
        if (z >= hi) return
        value = bd_mass(lo,z,par)
        if (value > -scifort_log_two) then
            tail = bd_mass(z,hi,par)
            value = ig_complement(tail)
        end if
    end function bd_lcdf

    pure function bd_lsf(x, par) result(value)
        implicit none
        real(dp), intent(in) :: x !! physical right-censoring threshold
        real(dp), intent(in) :: par(:) !! prepared bounded family parameters
        real(dp) :: value, z, lo, hi, tail
        z = location_scale_z(x,par(2:3))
        call bd_support(par,lo,hi)
        value = 0.0_dp
        if (z <= lo) return
        value = negative_infinity(x)
        if (z >= hi) return
        value = bd_mass(z,hi,par)
        if (value > -scifort_log_two) then
            tail = bd_mass(lo,z,par)
            value = ig_complement(tail)
        end if
    end function bd_lsf

    pure function bd_interval(a, b, par) result(value)
        implicit none
        real(dp), intent(in) :: a, b !! ordered finite physical interval endpoints
        real(dp), intent(in) :: par(:) !! prepared bounded family parameters
        real(dp) :: value, za, zb, lo, hi, lw, margin
        za = location_scale_z(a,par(2:3))
        zb = location_scale_z(b,par(2:3))
        call bd_support(par,lo,hi)
        value = negative_infinity(a)
        if (za >= hi .or. zb <= lo) return
        if (za <= lo .and. zb >= hi) then
            value = 0.0_dp
            return
        end if
        if (za > lo .and. zb < hi) then
            ! Retain physical width even when the two standardized coordinates coincide.
            lw = log_distance(b,a)-log(par(3))
            margin = min(za-lo,hi-zb)
            if (lw < log(margin)-log(1.0e5_dp)) then
                value = bd_gauss(za,lw,par)
                return
            end if
        end if
        value = bd_mass(max(lo,za),min(hi,zb),par)
    end function bd_interval

    pure function gausshyper_censored_loglikelihood(data, a, b, c, z, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: a !! first beta shape in [0.1,20]
        real(dp), intent(in) :: b !! second beta shape in [0.1,20]
        real(dp), intent(in) :: c !! tilt exponent in [-10,10]
        real(dp), intent(in) :: z !! tilt parameter in [-0.5,10]
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive normal scale, default one
        real(dp) :: value, params(40)
        logical :: valid
        params = 0.0_dp
        params(1:3) = [1.0_dp,0.0_dp,1.0_dp]
        params(4) = a
        params(5) = b
        params(6) = c
        params(7) = z
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        value = quiet_nan(0.0_dp)
        call bd_prepare(params,valid)
        if (.not. valid) return
        value = censored_loglikelihood(data,params,bd_lpdf,bd_lcdf,bd_lsf,bd_interval)
    end function gausshyper_censored_loglikelihood

    pure function gausshyper_censored_nnlf(data, a, b, c, z, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: a !! first beta shape in [0.1,20]
        real(dp), intent(in) :: b !! second beta shape in [0.1,20]
        real(dp), intent(in) :: c !! tilt exponent in [-10,10]
        real(dp), intent(in) :: z !! tilt parameter in [-0.5,10]
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive normal scale, default one
        real(dp) :: value
        value = -gausshyper_censored_loglikelihood(data,a, b, c, z,loc,scale)
    end function gausshyper_censored_nnlf

    pure function irwinhall_censored_loglikelihood(data, n, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: n !! integer number of uniforms in [1,32]
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive normal scale, default one
        real(dp) :: value, params(40)
        logical :: valid
        params = 0.0_dp
        params(1:3) = [2.0_dp,0.0_dp,1.0_dp]
        params(4) = n
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        value = quiet_nan(0.0_dp)
        call bd_prepare(params,valid)
        if (.not. valid) return
        value = censored_loglikelihood(data,params,bd_lpdf,bd_lcdf,bd_lsf,bd_interval)
    end function irwinhall_censored_loglikelihood

    pure function irwinhall_censored_nnlf(data, n, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: n !! integer number of uniforms in [1,32]
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive normal scale, default one
        real(dp) :: value
        value = -irwinhall_censored_loglikelihood(data,n,loc,scale)
    end function irwinhall_censored_nnlf

    pure function vonmises_line_censored_loglikelihood(data, kappa, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: kappa !! concentration in [0,100]
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive normal scale, default one
        real(dp) :: value, params(40)
        logical :: valid
        params = 0.0_dp
        params(1:3) = [3.0_dp,0.0_dp,1.0_dp]
        params(4) = kappa
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        value = quiet_nan(0.0_dp)
        call bd_prepare(params,valid)
        if (.not. valid) return
        value = censored_loglikelihood(data,params,bd_lpdf,bd_lcdf,bd_lsf,bd_interval)
    end function vonmises_line_censored_loglikelihood

    pure function vonmises_line_censored_nnlf(data, kappa, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: kappa !! concentration in [0,100]
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive normal scale, default one
        real(dp) :: value
        value = -vonmises_line_censored_loglikelihood(data,kappa,loc,scale)
    end function vonmises_line_censored_nnlf


    ! Independent gamma and kappa transformations; see CODE_PROVENANCE.md.
    pure function gk_gamma_mass(ly, ld, shape) result(value)
        implicit none
        real(dp), intent(in) :: ly !! log lower gamma coordinate, or negative infinity
        real(dp), intent(in) :: ld !! log positive gamma-coordinate interval width
        real(dp), intent(in) :: shape !! gamma shape in [1e-6,1e6]
        real(dp) :: value, lq, lh, q, t(3), d(3), r(3), correction, lp, lsf
        if (.not. ieee_is_finite(ly)) then
            value = loggamma_lcdf(ld,[shape,0.0_dp,1.0_dp])
            return
        end if
        lq = ld-ly
        if (lq+max(0.0_dp,log(shape),ly) <= log(0.001_dp)) then
            q = exp(lq)
            t = [0.5_dp,0.5_dp-sqrt(3.0_dp/5.0_dp)/2.0_dp, &
                0.5_dp+sqrt(3.0_dp/5.0_dp)/2.0_dp]
            d = log1p_safe(t*q)
            r = -d-exp(ld)*t
            if (q > 0.0_dp) then
                r = r+exp(log(shape)+lq)*(d/q)
            else
                r = r+exp(log(shape)+lq)*t
            end if
            correction = (4.0_dp/9.0_dp)*expm1_safe(r(1)) &
                +(5.0_dp/18.0_dp)*(expm1_safe(r(2))+expm1_safe(r(3)))
            value = loggamma_lpdf(ly,[shape,0.0_dp,1.0_dp])-ly+ld &
                +log1p_safe(correction)
            return
        end if
        if (ly > log(max(1.0_dp,shape))+log(1.0e16_dp) .and. ld > log(1.0e8_dp)) then
            value = loggamma_lsf(ly,[shape,0.0_dp,1.0_dp])
            return
        end if
        lh = bounded_logadd(ly,ld)
        lp = loggamma_lcdf(lh,[shape,0.0_dp,1.0_dp])
        lsf = loggamma_lsf(ly,[shape,0.0_dp,1.0_dp])
        if (lp <= lsf) then
            value = log_difference(lp,loggamma_lcdf(ly,[shape,0.0_dp,1.0_dp]))
        else
            value = log_difference(lsf,loggamma_lsf(lh,[shape,0.0_dp,1.0_dp]))
        end if
    end function gk_gamma_mass

    pure function ag_logy(z, params) result(value)
        implicit none
        real(dp), intent(in) :: z !! standardized point strictly between zero and one
        real(dp), intent(in) :: params(:) !! chi, loc, scale, log A, log gamma normalizer
        real(dp) :: value
        value = params(4)+log1p_safe(-z)+log1p_safe(z)
    end function ag_logy

    pure function ag_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! chi, loc, scale, log A, log gamma normalizer
        real(dp) :: value, z, ly
        value = negative_infinity(0.0_dp)
        z = location_scale_z(x,params(2:3))
        if (z <= 0.0_dp .or. z >= 1.0_dp) return
        ly = ag_logy(z,params)
        value = loggamma_lpdf(ly,[1.5_dp,0.0_dp,1.0_dp])-ly &
            +2.0_dp*log(params(1))+log(z)-log(params(3))-params(5)
    end function ag_lpdf

    pure function ag_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! chi, loc, scale, log A, log gamma normalizer
        real(dp) :: value, z, ly, other
        z = location_scale_z(x,params(2:3))
        value = negative_infinity(0.0_dp)
        if (z <= 0.0_dp) return
        value = 0.0_dp
        if (z >= 1.0_dp) return
        ly = ag_logy(z,params)
        other = min(0.0_dp,loggamma_lcdf(ly,[1.5_dp,0.0_dp,1.0_dp])-params(5))
        if (other < -scifort_log_two) then
            value = ig_complement(other)
        else
            value = min(0.0_dp,gk_gamma_mass(ly,params(4)+2.0_dp*log(z),1.5_dp)-params(5))
        end if
    end function ag_lcdf

    pure function ag_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! chi, loc, scale, log A, log gamma normalizer
        real(dp) :: value, z, other
        z = location_scale_z(x,params(2:3))
        value = 0.0_dp
        if (z <= 0.0_dp) return
        value = negative_infinity(0.0_dp)
        if (z >= 1.0_dp) return
        value = min(0.0_dp,loggamma_lcdf(ag_logy(z,params),[1.5_dp,0.0_dp,1.0_dp])-params(5))
        if (value > -scifort_log_two) then
            other = gk_gamma_mass(ag_logy(z,params),params(4)+2.0_dp*log(z),1.5_dp)-params(5)
            value = ig_complement(min(0.0_dp,other))
        end if
    end function ag_lsf

    pure function ag_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper endpoint, greater than a
        real(dp), intent(in) :: params(:) !! chi, loc, scale, log A, log gamma normalizer
        real(dp) :: value, za, zb, ly, ld
        za = location_scale_z(a,params(2:3))
        zb = location_scale_z(b,params(2:3))
        value = negative_infinity(0.0_dp)
        if (zb <= 0.0_dp .or. za >= 1.0_dp) return
        if (za <= 0.0_dp) then
            value = ag_lcdf(b,params)
        else if (zb >= 1.0_dp) then
            value = ag_lsf(a,params)
        else
            ly = ag_logy(zb,params)
            ld = params(4)+log(za+zb)+log_distance(b,a)-log(params(3))
            value = min(0.0_dp,gk_gamma_mass(ly,ld,1.5_dp)-params(5))
        end if
    end function ag_interval

    pure subroutine pe_coordinate(x, params, ly, region)
        implicit none
        real(dp), intent(in) :: x !! finite physical observation or censoring limit
        real(dp), intent(in) :: params(:) !! skew, loc, scale, beta, gamma shape
        real(dp), intent(out) :: ly !! log positive gamma coordinate in the interior
        integer, intent(out) :: region !! -1 outside support, 0 at endpoint, 1 interior
        real(dp) :: z, w
        z = location_scale_z(x,params(2:3))
        w = z+params(4)
        region = -1
        ly = negative_infinity(0.0_dp)
        if (w == 0.0_dp) then
            region = 0
            return
        end if
        if (sign(1.0_dp,w) /= sign(1.0_dp,params(4))) return
        region = 1
        if (ieee_is_finite(w)) then
            ly = log(abs(params(4)))+log(abs(w))
        else
            ly = log(abs(params(4))) &
                +log_distance(max(x,params(2)),min(x,params(2)))-log(params(3))
        end if
    end subroutine pe_coordinate

    pure function pe_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! skew, loc, scale, beta, gamma shape
        real(dp) :: value, ly, z
        integer :: region
        if (params(4) == 0.0_dp) then
            z = location_scale_z(x,params(2:3))
            value = negative_infinity(0.0_dp)
            if (abs(z) > sqrt(huge(1.0_dp))*sqrt(2.0_dp)) return
            value = normal_logpdf(z)-log(params(3))
            return
        end if
        call pe_coordinate(x,params,ly,region)
        value = negative_infinity(0.0_dp)
        if (region < 0) return
        if (region == 0) then
            if (params(5) < 1.0_dp) value = positive_infinity(0.0_dp)
            if (params(5) == 1.0_dp) value = log(abs(params(4)))-log(params(3))
            return
        end if
        value = loggamma_lpdf(ly,[params(5),0.0_dp,1.0_dp])-ly &
            +log(abs(params(4)))-log(params(3))
    end function pe_lpdf

    pure function pe_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! skew, loc, scale, beta, gamma shape
        real(dp) :: value, ly
        integer :: region
        if (params(4) == 0.0_dp) then
            value = normal_logcdf(location_scale_z(x,params(2:3)))
            return
        end if
        call pe_coordinate(x,params,ly,region)
        value = negative_infinity(0.0_dp)
        if (params(4) < 0.0_dp) value = 0.0_dp
        if (region <= 0) return
        if (params(4) > 0.0_dp) then
            value = loggamma_lcdf(ly,[params(5),0.0_dp,1.0_dp])
        else
            value = loggamma_lsf(ly,[params(5),0.0_dp,1.0_dp])
        end if
    end function pe_lcdf

    pure function pe_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! skew, loc, scale, beta, gamma shape
        real(dp) :: value, ly
        integer :: region
        if (params(4) == 0.0_dp) then
            value = normal_logsf(location_scale_z(x,params(2:3)))
            return
        end if
        call pe_coordinate(x,params,ly,region)
        value = 0.0_dp
        if (params(4) < 0.0_dp) value = negative_infinity(0.0_dp)
        if (region <= 0) return
        if (params(4) > 0.0_dp) then
            value = loggamma_lsf(ly,[params(5),0.0_dp,1.0_dp])
        else
            value = loggamma_lcdf(ly,[params(5),0.0_dp,1.0_dp])
        end if
    end function pe_lsf

    pure function pe_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper endpoint, greater than a
        real(dp), intent(in) :: params(:) !! skew, loc, scale, beta, gamma shape
        real(dp) :: value, ly, other, ld, za, zb, lw
        integer :: ra, rb
        if (params(4) == 0.0_dp) then
            za = location_scale_z(a,params(2:3))
            zb = location_scale_z(b,params(2:3))
            lw = log_distance(b,a)-log(params(3))
            if (za > 12.0_dp) then
                if (log(za)+lw > log(1.0e8_dp)) then
                    value = pe_lsf(a,params)
                    return
                end if
            else if (zb < -12.0_dp) then
                if (log(-zb)+lw > log(1.0e8_dp)) then
                    value = pe_lcdf(b,params)
                    return
                end if
            end if
            value = normal_small_interval(za,lw)
            if (ieee_is_nan(value)) value = tail_interval(a,b,params,pe_lcdf,pe_lsf)
            return
        end if
        call pe_coordinate(a,params,ly,ra)
        call pe_coordinate(b,params,other,rb)
        if (ra <= 0 .or. rb <= 0) then
            value = tail_interval(a,b,params,pe_lcdf,pe_lsf)
            return
        end if
        if (params(4) < 0.0_dp) ly = other
        ld = log(abs(params(4)))+log_distance(b,a)-log(params(3))
        value = gk_gamma_mass(ly,ld,params(5))
    end function pe_interval

    pure subroutine k4_terms(x, params, lt, lu, lv, lh, region)
        implicit none
        real(dp), intent(in) :: x !! finite observation or censoring limit
        real(dp), intent(in) :: params(:) !! h, loc, scale, k, lower bound, upper bound
        real(dp), intent(out) :: lt !! log inner transform T
        real(dp), intent(out) :: lu !! log inner support factor
        real(dp), intent(out) :: lv !! log outer support factor
        real(dp), intent(out) :: lh !! log reverse cumulative hazard
        integer, intent(out) :: region !! -1 below support, 0 interior, 1 above support
        real(dp) :: z, outer, lw, ld
        z = location_scale_z(x,params(2:3))
        region = -1
        lt = 0.0_dp
        lu = 0.0_dp
        lv = 0.0_dp
        lh = 0.0_dp
        if (params(1) > 0.0_dp .or. params(4) < 0.0_dp) then
            if (z <= params(5)) return
        end if
        region = 1
        if (params(4) > 0.0_dp) then
            if (z >= params(6)) return
        end if
        call gev_terms(x,[params(4),params(2),params(3)],lt,lu,region)
        if (region /= 0) return
        lh = lt
        if (params(1) == 0.0_dp) return
        outer = log(abs(params(1)))+lt
        if (params(1) < 0.0_dp) then
            lv = log1pexp(outer)
            lh = burr_logsoftplus(outer)-log(-params(1))
        else
            if (outer > -0.01_dp) then
                ! Retain distance from the finite outer endpoint when h*T rounds to one.
                lw = log_distance(z,params(5))+params(4)*log(params(1))
                ld = lw
                if (params(4) < 0.0_dp) then
                    ld = burr_logsoftplus(log(-params(4))+lw)-log(-params(4))
                else if (params(4) > 0.0_dp) then
                    ld = log(-ig_complement(log(params(4))+lw))-log(params(4))
                end if
                lv = log_exponential_mass(ld)
            else
                lv = ig_complement(outer)
            end if
            lh = log(-lv)-log(params(1))
        end if
    end subroutine k4_terms

    pure function k4_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! h, loc, scale, k, lower bound, upper bound
        real(dp) :: value, lt, lu, lv, lh, z, hk
        integer :: region
        call k4_terms(x,params,lt,lu,lv,lh,region)
        value = negative_infinity(0.0_dp)
        if (region /= 0) then
            z = location_scale_z(x,params(2:3))
            if (z == params(5)) then
                if (params(1) > 0.0_dp) then
                    if (params(1) == 1.0_dp) value = -log(params(3))
                    if (params(1) > 1.0_dp) value = positive_infinity(0.0_dp)
                else if (params(1) < 0.0_dp .and. params(4) < 0.0_dp) then
                    hk = params(1)*params(4)
                    if (hk == 1.0_dp) value = (1.0_dp/params(1)-1.0_dp)*log(-params(1)) &
                        -log(params(3))
                    if (hk > 1.0_dp) value = positive_infinity(0.0_dp)
                end if
            else if (z == params(6)) then
                if (params(4) == 1.0_dp) value = -log(params(3))
                if (params(4) > 1.0_dp) value = positive_infinity(0.0_dp)
            end if
            return
        end if
        if (lh > log(huge(1.0_dp))) return
        if (params(1) == 0.0_dp) then
            value = negative_hazard(lh)+lt-lu-log(params(3))
        else if (params(1) < 0.0_dp) then
            value = negative_hazard(lh)-log(-params(1)) &
                -log1pexp(-log(-params(1))-lt)-lu-log(params(3))
        else
            ! -H-lv=(1/h-1)*lv avoids subtracting large opposite terms.
            value = (1.0_dp/params(1)-1.0_dp)*lv+lt-lu-log(params(3))
        end if
    end function k4_lpdf

    pure function k4_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! h, loc, scale, k, lower bound, upper bound
        real(dp) :: value, lt, lu, lv, lh
        integer :: region
        call k4_terms(x,params,lt,lu,lv,lh,region)
        value = negative_infinity(0.0_dp)
        if (region < 0) return
        value = 0.0_dp
        if (region > 0) return
        value = negative_hazard(lh)
    end function k4_lcdf

    pure function k4_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! h, loc, scale, k, lower bound, upper bound
        real(dp) :: value, lt, lu, lv, lh
        integer :: region
        call k4_terms(x,params,lt,lu,lv,lh,region)
        value = 0.0_dp
        if (region < 0) return
        value = negative_infinity(0.0_dp)
        if (region > 0) return
        value = log_exponential_mass(lh)
    end function k4_lsf

    pure function k4_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper endpoint, greater than a
        real(dp), intent(in) :: params(:) !! h, loc, scale, k, lower bound, upper bound
        real(dp) :: value, ta, tb, ua, ub, va, vb, ha, hb, lw, lr, ld
        integer :: ra, rb
        call k4_terms(a,params,ta,ua,va,ha,ra)
        call k4_terms(b,params,tb,ub,vb,hb,rb)
        if (ra /= 0 .or. rb /= 0) then
            value = tail_interval(a,b,params,k4_lcdf,k4_lsf)
            return
        end if
        lw = log_distance(b,a)-log(params(3))
        if (params(4) == 0.0_dp) then
            ld = ta+log_exponential_mass(lw)
        else
            if (params(4) > 0.0_dp) then
                lr = log(params(4))+lw-ub
            else
                lr = log(-params(4))+lw-ua
            end if
            lr = burr_logsoftplus(lr)-log(abs(params(4)))
            ld = ta+log_exponential_mass(lr)
        end if
        if (params(1) /= 0.0_dp) then
            if (params(1) > 0.0_dp) then
                lr = log(params(1))+ld-va
            else
                lr = log(-params(1))+ld-vb
            end if
            ld = burr_logsoftplus(lr)-log(abs(params(1)))
        end if
        value = negative_hazard(hb)+log_exponential_mass(ld)
    end function k4_interval

    pure function argus_censored_loglikelihood(data, chi, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: chi !! positive ARGUS shape in [1e-6,100]
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(5)
        value = quiet_nan(0.0_dp)
        params = [chi,0.0_dp,1.0_dp,0.0_dp,0.0_dp]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        if (.not. all(ieee_is_finite(params))) return
        if (params(3) <= 0.0_dp) return
        if (chi < 1.0e-6_dp .or. chi > 100.0_dp) return
        params(4) = 2.0_dp*log(chi)-scifort_log_two
        params(5) = loggamma_lcdf(params(4),[1.5_dp,0.0_dp,1.0_dp])
        value = censored_loglikelihood(data,params,ag_lpdf,ag_lcdf, &
            ag_lsf,ag_interval)
    end function argus_censored_loglikelihood

    pure function argus_censored_nnlf(data, chi, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: chi !! positive ARGUS shape in [1e-6,100]
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -argus_censored_loglikelihood(data,chi,loc,scale)
    end function argus_censored_nnlf

    pure function pearson3_censored_loglikelihood(data, skew, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: skew !! skewness; normal if abs(skew)<1.6e-5, otherwise 0.002<=abs(skew)<=2000
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(5)
        value = quiet_nan(0.0_dp)
        params = [skew,0.0_dp,1.0_dp,0.0_dp,0.0_dp]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        if (.not. all(ieee_is_finite(params))) return
        if (params(3) <= 0.0_dp) return
        if (abs(skew) >= 1.6e-5_dp) then
            if (abs(skew) < 0.002_dp .or. abs(skew) > 2000.0_dp) return
            params(4) = 2.0_dp/skew
            params(5) = params(4)**2
        end if
        value = censored_loglikelihood(data,params,pe_lpdf,pe_lcdf, &
            pe_lsf,pe_interval)
    end function pearson3_censored_loglikelihood

    pure function pearson3_censored_nnlf(data, skew, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: skew !! skewness; normal if abs(skew)<1.6e-5, otherwise 0.002<=abs(skew)<=2000
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -pearson3_censored_loglikelihood(data,skew,loc,scale)
    end function pearson3_censored_nnlf

    pure function kappa4_censored_loglikelihood(data, h, k, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: h !! first kappa shape in [-10,10], zero or absolute value at least 1e-6
        real(dp), intent(in) :: k !! second kappa shape in [-10,10], zero or absolute value at least 1e-6
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(6)
        value = quiet_nan(0.0_dp)
        params = [h,0.0_dp,1.0_dp,k,0.0_dp,0.0_dp]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        if (.not. all(ieee_is_finite(params))) return
        if (params(3) <= 0.0_dp) return
        if (abs(h) > 10.0_dp .or. abs(k) > 10.0_dp) return
        if (h /= 0.0_dp .and. abs(h) < 1.0e-6_dp) return
        if (k /= 0.0_dp .and. abs(k) < 1.0e-6_dp) return
        if (h > 0.0_dp) then
            if (k*log(h) > log(10.0_dp)) return
        end if
        params(5) = 0.0_dp
        params(6) = 0.0_dp
        if (h > 0.0_dp) then
            if (k == 0.0_dp) then
                params(5) = log(h)
            else
                params(5) = -expm1_safe(-k*log(h))/k
            end if
        else if (k < 0.0_dp) then
            params(5) = 1.0_dp/k
        end if
        if (k > 0.0_dp) params(6) = 1.0_dp/k
        if (h > 0.0_dp .and. k > 0.0_dp) then
            if (params(5) >= params(6)) return
        end if
        value = censored_loglikelihood(data,params,k4_lpdf,k4_lcdf, &
            k4_lsf,k4_interval)
    end function kappa4_censored_loglikelihood

    pure function kappa4_censored_nnlf(data, h, k, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: h !! first kappa shape in [-10,10], zero or absolute value at least 1e-6
        real(dp), intent(in) :: k !! second kappa shape in [-10,10], zero or absolute value at least 1e-6
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -kappa4_censored_loglikelihood(data,h, k,loc,scale)
    end function kappa4_censored_nnlf


    ! Independent transformations of the cited densities; see CODE_PROVENANCE.md.
    pure function alpha_coordinate(x, params) result(w)
        implicit none
        real(dp), intent(in) :: x !! finite point above location
        real(dp), intent(in) :: params(:) !! shape, loc, scale
        real(dp) :: w, lr
        lr = log(params(3))-log_distance(x,params(2))
        w = negative_infinity(0.0_dp)
        if (lr > log(huge(1.0_dp))) return
        w = params(1)-exp(lr)
    end function alpha_coordinate

    pure function alpha_lpdf_c(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! shape, loc, scale
        real(dp) :: value, w, ld
        value = negative_infinity(0.0_dp)
        if (x <= params(2)) return
        w = alpha_coordinate(x,params)
        if (abs(w) > sqrt(huge(1.0_dp))*sqrt(2.0_dp)) return
        ld = log_distance(x,params(2))
        value = (-0.5_dp*w)*w-scifort_log_sqrt_two_pi-normal_logcdf(params(1)) &
            +log(params(3))-2.0_dp*ld
    end function alpha_lpdf_c

    pure function alpha_lsf_c(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite censoring limit
        real(dp), intent(in) :: params(:) !! shape, loc, scale
        real(dp) :: value, w, lr
        value = 0.0_dp
        if (x <= params(2)) return
        w = alpha_coordinate(x,params)
        lr = log(params(3))-log_distance(x,params(2))
        value = normal_small_interval(w,lr)
        if (ieee_is_nan(value)) then
            if (w < -sqrt(huge(1.0_dp))*sqrt(2.0_dp)) then
                value = normal_logcdf(params(1))
            else
                if (w >= 0.0_dp) then
                    value = log_difference(normal_logsf(w),normal_logsf(params(1)))
                else
                    value = log_difference(normal_logcdf(params(1)),normal_logcdf(w))
                end if
            end if
        end if
        value = min(0.0_dp,value-normal_logcdf(params(1)))
    end function alpha_lsf_c

    pure function alpha_lcdf_c(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite censoring limit
        real(dp), intent(in) :: params(:) !! shape, loc, scale
        real(dp) :: value, w, q
        value = negative_infinity(0.0_dp)
        if (x <= params(2)) return
        q = alpha_lsf_c(x,params)
        if (q < -scifort_log_two) then
            value = ig_complement(q)
        else
            w = alpha_coordinate(x,params)
            if (w < -sqrt(huge(1.0_dp))*sqrt(2.0_dp)) return
            value = min(0.0_dp,normal_logcdf(w)-normal_logcdf(params(1)))
        end if
    end function alpha_lcdf_c

    pure function alpha_interval_c(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper endpoint, greater than a
        real(dp), intent(in) :: params(:) !! shape, loc, scale
        real(dp) :: value, w, lw, upper_w
        value = negative_infinity(0.0_dp)
        if (b <= params(2)) return
        if (a <= params(2)) then
            value = alpha_lcdf_c(b,params)
            return
        end if
        w = alpha_coordinate(a,params)
        lw = log(params(3))+log_distance(b,a)-log_distance(a,params(2)) &
            -log_distance(b,params(2))
        upper_w = alpha_coordinate(b,params)
        if (upper_w < -12.0_dp) then
            ! If the normal log density changes by more than 1e8, the lower
            ! endpoint contribution is negligible even when rounded logs coincide.
            if (log(-upper_w)+lw > log(1.0e8_dp)) then
                value = alpha_lcdf_c(b,params)
                return
            end if
        end if
        value = normal_small_interval(w,lw)
        if (ieee_is_nan(value)) then
            value = tail_interval(a,b,params,alpha_lcdf_c,alpha_lsf_c)
        else
            value = value-normal_logcdf(params(1))
        end if
    end function alpha_interval_c

    pure subroutine gev_terms(x, params, lt, lu, region)
        implicit none
        real(dp), intent(in) :: x !! finite observation or censoring limit
        real(dp), intent(in) :: params(:) !! shape, loc, scale
        real(dp), intent(out) :: lt !! log transformed reverse hazard in the interior
        real(dp), intent(out) :: lu !! log positive support factor in the interior
        integer, intent(out) :: region !! -1 below support, 0 interior, 1 above support
        real(dp) :: z, cz, lc
        z = location_scale_z(x,params(2:3))
        region = 0
        lu = 0.0_dp
        lt = -z
        if (params(1) == 0.0_dp) return
        if (z == 0.0_dp) return
        if (ieee_is_finite(z)) then
            lc = log(abs(params(1)))+log(abs(z))
        else
            lc = log(abs(params(1)))+log_distance(max(x,params(2)),min(x,params(2))) &
                -log(params(3))
        end if
        if (.not. ieee_is_finite(z) .or. lc > log(huge(1.0_dp))-1.0_dp) then
            if (sign(1.0_dp,z) == sign(1.0_dp,params(1))) then
                region = int(sign(1.0_dp,params(1)))
                return
            end if
            lu = log1pexp(lc)
        else
            cz = params(1)*z
            if (cz >= 1.0_dp) then
                region = int(sign(1.0_dp,params(1)))
                return
            end if
            lu = log1p_safe(-cz)
        end if
        lt = lu/params(1)
    end subroutine gev_terms

    pure function gev_lpdf_c(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! shape, loc, scale
        real(dp) :: value, lt, lu, z
        integer :: region
        call gev_terms(x,params,lt,lu,region)
        value = negative_infinity(0.0_dp)
        if (region /= 0) then
            ! Preserve the scalar API's density convention at an exact upper endpoint.
            if (params(1) > 0.0_dp) then
                z = location_scale_z(x,params(2:3))
                if (z == 1.0_dp/params(1)) then
                    if (params(1) == 1.0_dp) value = -log(params(3))
                    if (params(1) > 1.0_dp) value = positive_infinity(0.0_dp)
                end if
            end if
            return
        end if
        if (lt > log(huge(1.0_dp))) return
        value = negative_hazard(lt)+lt-lu-log(params(3))
    end function gev_lpdf_c

    pure function gev_lcdf_c(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite censoring limit
        real(dp), intent(in) :: params(:) !! shape, loc, scale
        real(dp) :: value, lt, lu
        integer :: region
        call gev_terms(x,params,lt,lu,region)
        value = negative_infinity(0.0_dp)
        if (region < 0) return
        value = 0.0_dp
        if (region > 0) return
        value = negative_hazard(lt)
    end function gev_lcdf_c

    pure function gev_lsf_c(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite censoring limit
        real(dp), intent(in) :: params(:) !! shape, loc, scale
        real(dp) :: value, lt, lu
        integer :: region
        call gev_terms(x,params,lt,lu,region)
        value = 0.0_dp
        if (region < 0) return
        value = negative_infinity(0.0_dp)
        if (region > 0) return
        value = log_exponential_mass(lt)
    end function gev_lsf_c

    pure function gev_interval_c(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper endpoint, greater than a
        real(dp), intent(in) :: params(:) !! shape, loc, scale
        real(dp) :: value, la, lb, ua, ub, lw, lr, ld, delta
        integer :: ra, rb
        call gev_terms(a,params,la,ua,ra)
        call gev_terms(b,params,lb,ub,rb)
        if (ra /= 0 .or. rb /= 0) then
            value = tail_interval(a,b,params,gev_lcdf_c,gev_lsf_c)
            return
        end if
        lw = log_distance(b,a)-log(params(3))
        if (params(1) == 0.0_dp) then
            ld = la+log_exponential_mass(lw)
        else
            ! The reverse hazard decreases: retain log(T(a)-T(b)) from physical width.
            if (params(1) > 0.0_dp) then
                lr = log(params(1))+lw-ub
            else
                lr = log(-params(1))+lw-ua
            end if
            ! log(delta), delta=abs(log(u(a)/u(b))/c).
            lr = burr_logsoftplus(lr)-log(abs(params(1)))
            if (lr < log(epsilon(1.0_dp))-2.0_dp) then
                ld = la+lr
            else if (lr > log(huge(1.0_dp))) then
                ld = la
            else
                delta = exp(lr)
                ld = la+log(-expm1_safe(-delta))
            end if
        end if
        value = negative_hazard(lb)+log_exponential_mass(ld)
    end function gev_interval_c

    pure function cb_logratio(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! physical point at or below the tail transition
        real(dp), intent(in) :: params(:) !! beta, loc, scale, m, log denominator, log tail weight
        real(dp) :: value, z, ld
        value = 0.0_dp
        z = location_scale_z(x,params(2:3))
        if (z == -params(1)) return
        if (ieee_is_finite(z)) then
            ld = log_distance(-params(1),z)
        else
            ld = log_distance(params(2),x)-log(params(3))
        end if
        value = log1pexp(log(params(1))-log(params(4))+ld)
    end function cb_logratio

    pure function cb_lpdf_c(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! beta, loc, scale, m, log denominator, log tail weight
        real(dp) :: value, z
        z = location_scale_z(x,params(2:3))
        value = negative_infinity(0.0_dp)
        if (z > -params(1)) then
            if (abs(z) > sqrt(huge(1.0_dp))*sqrt(2.0_dp)) return
            value = (-0.5_dp*z)*z-params(5)-log(params(3))
        else
            value = -params(5)-0.5_dp*params(1)**2 &
                -params(4)*cb_logratio(x,params)-log(params(3))
        end if
    end function cb_lpdf_c

    pure function cb_lcdf_c(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite censoring limit
        real(dp), intent(in) :: params(:) !! beta, loc, scale, m, log denominator, log tail weight
        real(dp) :: value, z, q, core
        z = location_scale_z(x,params(2:3))
        if (z <= -params(1)) then
            value = params(6)-(params(4)-1.0_dp)*cb_logratio(x,params)
        else
            q = cb_lsf_c(x,params)
            if (q < -scifort_log_two) then
                value = ig_complement(q)
            else
                core = normal_interval(-params(1),z,[0.0_dp,1.0_dp]) &
                    +scifort_log_sqrt_two_pi-params(5)
                value = min(0.0_dp,bounded_logadd(params(6),core))
            end if
        end if
    end function cb_lcdf_c

    pure function cb_lsf_c(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite censoring limit
        real(dp), intent(in) :: params(:) !! beta, loc, scale, m, log denominator, log tail weight
        real(dp) :: value, z, p, core
        z = location_scale_z(x,params(2:3))
        if (z > -params(1)) then
            value = negative_infinity(0.0_dp)
            if (z > sqrt(huge(1.0_dp))*sqrt(2.0_dp)) return
            value = min(0.0_dp,scifort_log_sqrt_two_pi+normal_logsf(z)-params(5))
        else
            p = -(params(4)-1.0_dp)*cb_logratio(x,params)
            if (params(6)+p < -scifort_log_two) then
                value = ig_complement(params(6)+p)
                return
            end if
            core = scifort_log_sqrt_two_pi+normal_logcdf(params(1))-params(5)
            value = min(0.0_dp,bounded_logadd(core,params(6)+ig_complement(p)))
        end if
    end function cb_lsf_c

    pure function cb_interval_c(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper endpoint, greater than a
        real(dp), intent(in) :: params(:) !! beta, loc, scale, m, log denominator, log tail weight
        real(dp) :: value, za, zb, ld, lr, lm, left, right, lw
        za = location_scale_z(a,params(2:3))
        zb = location_scale_z(b,params(2:3))
        lw = log_distance(b,a)-log(params(3))
        if (za < -params(1) .and. zb > -params(1)) then
            if (lw+log(max(1.0_dp,params(1))) <= log(1.0e-6_dp)) then
                ! The joined density and slope agree. The Gaussian extension's
                ! relative error is at most O((1+beta**2/m)*width**2), < 1e-12.
                ! Retain the whole physical width instead of two rounded pieces.
                value = normal_small_interval(za,lw)+scifort_log_sqrt_two_pi-params(5)
                return
            end if
        end if
        if (za >= -params(1)) then
            if (za > 12.0_dp) then
                if (log(za)+log_distance(b,a)-log(params(3)) > log(1.0e8_dp)) then
                    value = cb_lsf_c(a,params)
                    return
                end if
            end if
            value = normal_small_interval(za,log_distance(b,a)-log(params(3)))
            if (ieee_is_nan(value)) &
                value = tail_interval(za,zb,[0.0_dp,1.0_dp],normal_lcdf,normal_lsf)
            value = value+scifort_log_sqrt_two_pi-params(5)
        else if (zb <= -params(1)) then
            ! log D(b), D=m/beta-beta-z, retaining the interval independently.
            ld = log(params(4))-log(params(1))+cb_logratio(b,params)
            lr = log_distance(b,a)-log(params(3))-ld
            lm = log(params(4)-1.0_dp)+burr_logsoftplus(lr)
            value = cb_lcdf_c(b,params)+log_exponential_mass(lm)
        else
            lr = cb_logratio(a,params)
            left = params(6)+log_exponential_mass(log(params(4)-1.0_dp)+log(lr))
            right = normal_interval(-params(1),zb,[0.0_dp,1.0_dp]) &
                +scifort_log_sqrt_two_pi-params(5)
            value = bounded_logadd(left,right)
        end if
    end function cb_interval_c

    pure function alpha_censored_loglikelihood(data, a, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: a !! positive alpha shape in [1e-6,100]
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(3)
        value = quiet_nan(0.0_dp)
        params = [a,0.0_dp,1.0_dp]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        if (.not. all(ieee_is_finite(params))) return
        if (params(3) <= 0.0_dp) return
        if (a < 1.0e-6_dp .or. a > 100.0_dp) return
        value = censored_loglikelihood(data,params,alpha_lpdf_c,alpha_lcdf_c, &
            alpha_lsf_c,alpha_interval_c)
    end function alpha_censored_loglikelihood

    pure function alpha_censored_nnlf(data, a, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: a !! positive alpha shape in [1e-6,100]
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -alpha_censored_loglikelihood(data,a,loc,scale)
    end function alpha_censored_nnlf

    pure function genextreme_censored_loglikelihood(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! GEV shape in [-100,100], zero or absolute value at least 1e-6
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(3)
        value = quiet_nan(0.0_dp)
        params = [c,0.0_dp,1.0_dp]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        if (.not. all(ieee_is_finite(params))) return
        if (params(3) <= 0.0_dp) return
        if (abs(c) > 100.0_dp) return
        if (c /= 0.0_dp .and. abs(c) < 1.0e-6_dp) return
        value = censored_loglikelihood(data,params,gev_lpdf_c,gev_lcdf_c, &
            gev_lsf_c,gev_interval_c)
    end function genextreme_censored_loglikelihood

    pure function genextreme_censored_nnlf(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! GEV shape in [-100,100], zero or absolute value at least 1e-6
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -genextreme_censored_loglikelihood(data,c,loc,scale)
    end function genextreme_censored_nnlf

    pure function crystalball_censored_loglikelihood(data, beta, m, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: beta !! positive transition magnitude in [1e-6,100]
        real(dp), intent(in) :: m !! power-law exponent in (1,1e6]
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(6)
        value = quiet_nan(0.0_dp)
        params = [beta,0.0_dp,1.0_dp,m,0.0_dp,0.0_dp]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        if (.not. all(ieee_is_finite(params))) return
        if (params(3) <= 0.0_dp) return
        if (beta < 1.0e-6_dp .or. beta > 100.0_dp) return
        if (m <= 1.0_dp .or. m > 1.0e6_dp) return
        params(6) = log(m)-log(beta)-log(m-1.0_dp)-0.5_dp*beta**2
        params(5) = bounded_logadd(params(6), &
            scifort_log_sqrt_two_pi+normal_logcdf(beta))
        if (params(6) >= scifort_log_sqrt_two_pi+normal_logcdf(beta)) then
            params(6) = -log1p_safe(exp(scifort_log_sqrt_two_pi+normal_logcdf(beta)-params(6)))
        else
            params(6) = params(6)-params(5)
        end if
        value = censored_loglikelihood(data,params,cb_lpdf_c,cb_lcdf_c, &
            cb_lsf_c,cb_interval_c)
    end function crystalball_censored_loglikelihood

    pure function crystalball_censored_nnlf(data, beta, m, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: beta !! positive transition magnitude in [1e-6,100]
        real(dp), intent(in) :: m !! power-law exponent in (1,1e6]
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -crystalball_censored_loglikelihood(data,beta, m,loc,scale)
    end function crystalball_censored_nnlf


    ! Independent transformed/conditional densities; see CODE_PROVENANCE.md.
    pure function gg_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! gamma shape, loc, scale, zero log rate, signed power
        real(dp) :: value
        value = negative_infinity(0.0_dp)
        if (x < params(2)) return
        if (x == params(2)) then
            if (params(5) < 0.0_dp) return
            if (params(1) < 1.0_dp/params(5)) value = positive_infinity(0.0_dp)
            if (params(1) == 1.0_dp/params(5)) value = log(params(5)) &
                -log_gamma(params(1))-log(params(3))
            return
        end if
        value = loggamma_lpdf(gamma_power_logy(x,params),[params(1),0.0_dp,1.0_dp]) &
            +log(abs(params(5)))-log_distance(x,params(2))
    end function gg_lpdf

    pure subroutine gg_tails(x, params, lp, lq)
        implicit none
        real(dp), intent(in) :: x !! finite censoring limit
        real(dp), intent(in) :: params(:) !! gamma shape, loc, scale, zero log rate, signed power
        real(dp), intent(out) :: lp !! log lower probability
        real(dp), intent(out) :: lq !! log upper probability
        real(dp) :: swap
        lp = negative_infinity(0.0_dp)
        lq = 0.0_dp
        if (x <= params(2)) return
        lq = gamma_power_lsf(x,params)
        if (params(1) <= epsilon(1.0_dp)/4096.0_dp) then
            lp = ig_complement(lq)
        else
            lp = loggamma_lcdf(gamma_power_logy(x,params),[params(1),0.0_dp,1.0_dp])
            if (lq < -scifort_log_two) lp = ig_complement(lq)
        end if
        if (params(5) < 0.0_dp) then
            swap = lp
            lp = lq
            lq = swap
        end if
    end subroutine gg_tails

    pure function gg_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! gamma shape, loc, scale, zero log rate, signed power
        real(dp) :: value, unused
        call gg_tails(x,params,value,unused)
    end function gg_lcdf

    pure function gg_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! gamma shape, loc, scale, zero log rate, signed power
        real(dp) :: value, unused
        call gg_tails(x,params,unused,value)
    end function gg_lsf

    pure function gg_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper endpoint greater than a
        real(dp), intent(in) :: params(:) !! gamma shape, loc, scale, zero log rate, signed power
        real(dp) :: value, width, logq, ly, lyb, lc, q, t(3), d(3), first(3), second(3), r(3)
        real(dp) :: correction, endpoint, rate, omitted
        value = negative_infinity(0.0_dp)
        if (b <= params(2)) return
        if (a <= params(2)) then
            value = gg_lcdf(b,params)
            return
        end if
        width = log_distance(b,a)
        logq = width-log_distance(a,params(2))
        ly = gamma_power_logy(a,params)
        lyb = gamma_power_logy(b,params)
        lc = log(abs(params(5)))
        if (logq+max(0.0_dp,lc+log(params(1)),lc+max(ly,lyb)) < log(0.001_dp)) then
            q = exp(logq)
            t = [0.5_dp,0.5_dp-sqrt(3.0_dp/5.0_dp)/2.0_dp, &
                0.5_dp+sqrt(3.0_dp/5.0_dp)/2.0_dp]
            d = log1p_safe(t*q)
            first = t
            second = params(5)*t
            if (q > 0.0_dp) then
                first = d/q
                second = expm1_safe(params(5)*d)/q
            end if
            r = sign(1.0_dp,params(5))*exp(log(params(1))+lc+logq)*first-d &
                -exp(ly+logq)*second
            correction = (4.0_dp/9.0_dp)*expm1_safe(r(1)) &
                +(5.0_dp/18.0_dp)*(expm1_safe(r(2))+expm1_safe(r(3)))
            value = gg_lpdf(a,params)+width+log1p_safe(correction)
        else
            ! Endpoint expansion when Y dominates both other density-slope terms.
            endpoint = a
            if (params(5) < 0.0_dp) endpoint = b
            ly = gamma_power_logy(endpoint,params)
            rate = lc+ly-log_distance(endpoint,params(2))
            if (logq < log(1.0e-6_dp) .and. ly > log(params(1))+log(1.0e16_dp) &
                    .and. lc+ly > log(1.0e16_dp) .and. rate+width > log(1.0e8_dp)) then
                value = gg_lpdf(endpoint,params)-rate+log_exponential_mass(rate+width)
            else
                value = tail_interval(a,b,params,gg_lcdf,gg_lsf)
                if (value > -scifort_log_two) then
                    omitted = bounded_logadd(gg_lcdf(a,params),gg_lsf(b,params))
                    value = log1p_safe(-exp(omitted))
                end if
            end if
        end if
    end function gg_interval

    pure subroutine trunc_setup(params, normal, valid)
        implicit none
        real(dp), intent(inout) :: params(:) !! power, a, b, loc, scale, physical bounds, normalizer
        logical, intent(in) :: normal !! true for normal; false for Weibull
        logical, intent(out) :: valid !! true when parameters and derived support are valid
        real(dp) :: product, width
        integer :: i
        valid = .false.
        if (.not. all(ieee_is_finite(params))) return
        if (params(3) <= params(2) .or. params(5) <= 0.0_dp) return
        if (normal) then
            if (maxval(abs(params(2:3))) > 100.0_dp) return
        else
            if (params(2) < 0.0_dp .or. params(1) < 1.0e-6_dp &
                    .or. params(1) > 1.0e6_dp) return
        end if
        if (abs(params(4)) > huge(1.0_dp)/4.0_dp) return
        do i = 2, 3
            if (abs(params(i)) > 1.0_dp) then
                if (params(5) > (huge(1.0_dp)/4.0_dp)/abs(params(i))) return
            end if
            product = params(5)*params(i)
            if (params(i) /= 0.0_dp .and. product == 0.0_dp) return
            if (product /= 0.0_dp .and. abs(product) < tiny(1.0_dp)) return
            if (abs(product) > huge(1.0_dp)/4.0_dp) return
            params(i+4) = params(4)+product
        end do
        if (params(7) <= params(6)) return
        if (params(6) /= 0.0_dp .and. abs(params(6)) < tiny(1.0_dp)) return
        if (params(7) /= 0.0_dp .and. abs(params(7)) < tiny(1.0_dp)) return
        width = params(3)-params(2)
        if (params(4) /= 0.0_dp) then
            if (log(abs(params(4)))-log(params(5))-log(width) > log(1.0e4_dp)) return
        end if
        if (width < 1.0e-4_dp*max(1.0_dp,abs(params(2)),abs(params(3)))) then
            ! Preserve narrow standardized windows under an exact binary scaling.
            if (params(4) /= 0.0_dp .or. fraction(params(5)) /= 0.5_dp) return
            if (params(6) /= 0.0_dp .and. abs(params(6)) < tiny(1.0_dp)) return
            if (abs(params(7)) < tiny(1.0_dp)) return
        end if
        if (normal) then
            params(8) = normal_interval(params(2),params(3),[0.0_dp,1.0_dp])
            if (.not. ieee_is_finite(params(8))) return
        else if (params(2) == 0.0_dp) then
            params(8) = params(1)*log(params(3))
        else
            params(8) = weibull_log_increment(params(2),params(3),[params(1),0.0_dp,1.0_dp])
        end if
        valid = .true.
    end subroutine trunc_setup

    pure function tn_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! unused power, a, b, loc, scale, physical bounds, log Z
        real(dp) :: value
        value = negative_infinity(0.0_dp)
        if (x < params(6) .or. x > params(7)) return
        value = normal_logpdf(x,params(4),params(5))-params(8)
    end function tn_lpdf

    pure function tn_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper interval endpoint greater than a
        real(dp), intent(in) :: params(:) !! unused power, a, b, loc, scale, physical bounds, log Z
        real(dp) :: value, lower, upper
        value = negative_infinity(0.0_dp)
        lower = max(a,params(6))
        upper = min(b,params(7))
        if (upper <= lower) return
        if (lower == params(6) .and. upper == params(7)) then
            value = 0.0_dp
        else
            value = min(0.0_dp,normal_interval(lower,upper,params(4:5))-params(8))
        end if
    end function tn_interval

    pure function tn_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! unused power, a, b, loc, scale, physical bounds, log Z
        real(dp) :: value, other
        value = negative_infinity(0.0_dp)
        if (x <= params(6)) return
        value = 0.0_dp
        if (x >= params(7)) return
        value = tn_interval(params(6),x,params)
        if (value > -scifort_log_two) then
            other = tn_interval(x,params(7),params)
            value = ig_complement(other)
        end if
    end function tn_lcdf

    pure function tn_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! unused power, a, b, loc, scale, physical bounds, log Z
        real(dp) :: value, other
        value = 0.0_dp
        if (x <= params(6)) return
        value = negative_infinity(0.0_dp)
        if (x >= params(7)) return
        value = tn_interval(x,params(7),params)
        if (value > -scifort_log_two) then
            other = tn_interval(params(6),x,params)
            value = ig_complement(other)
        end if
    end function tn_lsf

    pure function tw_increment(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! lower endpoint at or above parent loc
        real(dp), intent(in) :: b !! upper endpoint greater than a
        real(dp), intent(in) :: params(:) !! power, a, b, loc, scale, physical bounds, log hazard span
        real(dp) :: value
        if (a == params(4)) then
            value = weibull_log_hazard(b,[params(1),params(4),params(5)])
        else
            value = weibull_log_increment(a,b,[params(1),params(4),params(5)])
        end if
    end function tw_increment

    pure function tw_log_fraction(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! lower endpoint at or above parent location
        real(dp), intent(in) :: b !! upper endpoint greater than a
        real(dp), intent(in) :: params(:) !! power, a, b, loc, scale, physical bounds, log hazard span
        real(dp) :: value
        value = 0.0_dp
        if (a == params(4)) return
        value = log_exponential_mass(log(params(1))+burr_logsoftplus( &
            log_distance(b,a)-log_distance(a,params(4))))
    end function tw_log_fraction

    pure function tw_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! power, a, b, loc, scale, physical bounds, log hazard span
        real(dp) :: value, lz, shift, correction
        value = negative_infinity(0.0_dp)
        if (x < params(6) .or. x > params(7)) return
        if (x == params(4)) then
            if (params(1) < 1.0_dp) value = positive_infinity(0.0_dp)
            if (params(1) == 1.0_dp) value = -log(params(5))-log_exponential_mass(params(8))
            return
        end if
        lz = lognormal_logz(x,[params(1),params(4),params(5)])
        shift = 0.0_dp
        if (x > params(6)) shift = negative_hazard(tw_increment(params(6),x,params))
        if (params(8) < log(0.001_dp)) then
            lz = weibull_log_hazard(x,[params(1),params(4),params(7)-params(4)])
            correction = 0.0_dp
            if (params(8) > log(epsilon(1.0_dp))) &
                correction = log_exponential_mass(params(8))-params(8)
            value = log(params(1))-log_distance(x,params(4))+lz &
                -tw_log_fraction(params(6),params(7),params)+shift-correction
            return
        end if
        value = log(params(1))+(params(1)-1.0_dp)*lz-log(params(5)) &
            +shift-log_exponential_mass(params(8))
    end function tw_lpdf

    pure function tw_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper interval endpoint greater than a
        real(dp), intent(in) :: params(:) !! power, a, b, loc, scale, physical bounds, log hazard span
        real(dp) :: value, lower, upper, shift, increment, correction, normalization
        value = negative_infinity(0.0_dp)
        lower = max(a,params(6))
        upper = min(b,params(7))
        if (upper <= lower) return
        if (lower == params(6) .and. upper == params(7)) then
            value = 0.0_dp
            return
        end if
        shift = 0.0_dp
        if (lower > params(6)) shift = negative_hazard(tw_increment(params(6),lower,params))
        if (params(8) < log(0.001_dp)) then
            increment = tw_increment(lower,upper,params)
            correction = 0.0_dp
            normalization = 0.0_dp
            if (increment > log(epsilon(1.0_dp))) &
                correction = log_exponential_mass(increment)-increment
            if (params(8) > log(epsilon(1.0_dp))) &
                normalization = log_exponential_mass(params(8))-params(8)
            value = shift+weibull_log_hazard(upper,[params(1),params(4),params(7)-params(4)]) &
                +tw_log_fraction(lower,upper,params)-tw_log_fraction(params(6),params(7),params) &
                +correction-normalization
            value = min(0.0_dp,value)
            return
        end if
        value = min(0.0_dp,shift+log_exponential_mass(tw_increment(lower,upper,params)) &
            -log_exponential_mass(params(8)))
    end function tw_interval

    pure function tw_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! power, a, b, loc, scale, physical bounds, log hazard span
        real(dp) :: value, other
        value = negative_infinity(0.0_dp)
        if (x <= params(6)) return
        value = 0.0_dp
        if (x >= params(7)) return
        value = tw_interval(params(6),x,params)
        if (value > -scifort_log_two) then
            other = tw_interval(x,params(7),params)
            value = ig_complement(other)
        end if
    end function tw_lcdf

    pure function tw_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! power, a, b, loc, scale, physical bounds, log hazard span
        real(dp) :: value, other
        value = 0.0_dp
        if (x <= params(6)) return
        value = negative_infinity(0.0_dp)
        if (x >= params(7)) return
        value = tw_interval(x,params(7),params)
        if (value > -scifort_log_two) then
            other = tw_interval(params(6),x,params)
            value = ig_complement(other)
        end if
    end function tw_lsf

    pure function gengamma_censored_loglikelihood(data, a, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: a !! family shape or standardized bound; see API limits
        real(dp), intent(in) :: c !! family shape or standardized bound; see API limits
        real(dp), intent(in), optional :: loc !! finite parent location, default zero
        real(dp), intent(in), optional :: scale !! positive finite parent scale, default one
        real(dp) :: value, params(5)
        value = quiet_nan(0.0_dp)
        params = [a,0.0_dp,1.0_dp,0.0_dp,c]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        if (.not. all(ieee_is_finite(params))) return
        if (a < tiny(1.0_dp) .or. a > 1.0e6_dp) return
        if (abs(c) < 1.0e-6_dp .or. abs(c) > 1.0e6_dp .or. params(3) <= 0.0_dp) return
        value = censored_loglikelihood(data,params,gg_lpdf,gg_lcdf, &
            gg_lsf,gg_interval)
    end function gengamma_censored_loglikelihood

    pure function gengamma_censored_nnlf(data, a, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: a !! family shape or standardized bound; see API limits
        real(dp), intent(in) :: c !! family shape or standardized bound; see API limits
        real(dp), intent(in), optional :: loc !! finite parent location, default zero
        real(dp), intent(in), optional :: scale !! positive finite parent scale, default one
        real(dp) :: value
        value = -gengamma_censored_loglikelihood(data,a, c,loc,scale)
    end function gengamma_censored_nnlf

    pure function truncnorm_censored_loglikelihood(data, a, b, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: a !! family shape or standardized bound; see API limits
        real(dp), intent(in) :: b !! family shape or standardized bound; see API limits
        real(dp), intent(in), optional :: loc !! finite parent location, default zero
        real(dp), intent(in), optional :: scale !! positive finite parent scale, default one
        real(dp) :: value, params(8)
        logical :: valid
        value = quiet_nan(0.0_dp)
        params = [1.0_dp,a,b,0.0_dp,1.0_dp,0.0_dp,0.0_dp,0.0_dp]
        if (present(loc)) params(4) = loc
        if (present(scale)) params(5) = scale
        call trunc_setup(params,.true.,valid)
        if (.not. valid) return
        value = censored_loglikelihood(data,params,tn_lpdf,tn_lcdf, &
            tn_lsf,tn_interval)
    end function truncnorm_censored_loglikelihood

    pure function truncnorm_censored_nnlf(data, a, b, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: a !! family shape or standardized bound; see API limits
        real(dp), intent(in) :: b !! family shape or standardized bound; see API limits
        real(dp), intent(in), optional :: loc !! finite parent location, default zero
        real(dp), intent(in), optional :: scale !! positive finite parent scale, default one
        real(dp) :: value
        value = -truncnorm_censored_loglikelihood(data,a, b,loc,scale)
    end function truncnorm_censored_nnlf

    pure function truncweibull_min_censored_loglikelihood(data, c, a, b, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! family shape or standardized bound; see API limits
        real(dp), intent(in) :: a !! family shape or standardized bound; see API limits
        real(dp), intent(in) :: b !! family shape or standardized bound; see API limits
        real(dp), intent(in), optional :: loc !! finite parent location, default zero
        real(dp), intent(in), optional :: scale !! positive finite parent scale, default one
        real(dp) :: value, params(8)
        logical :: valid
        value = quiet_nan(0.0_dp)
        params = [c,a,b,0.0_dp,1.0_dp,0.0_dp,0.0_dp,0.0_dp]
        if (present(loc)) params(4) = loc
        if (present(scale)) params(5) = scale
        call trunc_setup(params,.false.,valid)
        if (.not. valid) return
        value = censored_loglikelihood(data,params,tw_lpdf,tw_lcdf, &
            tw_lsf,tw_interval)
    end function truncweibull_min_censored_loglikelihood

    pure function truncweibull_min_censored_nnlf(data, c, a, b, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! family shape or standardized bound; see API limits
        real(dp), intent(in) :: a !! family shape or standardized bound; see API limits
        real(dp), intent(in) :: b !! family shape or standardized bound; see API limits
        real(dp), intent(in), optional :: loc !! finite parent location, default zero
        real(dp), intent(in), optional :: scale !! positive finite parent scale, default one
        real(dp) :: value
        value = -truncweibull_min_censored_loglikelihood(data,c, a, b,loc,scale)
    end function truncweibull_min_censored_nnlf


    ! Independent density and normal-tail identities; see CODE_PROVENANCE.md.
    pure subroutine ig_coordinates(x, params, lz, la, lb, sign_a, logphi)
        implicit none
        real(dp), intent(in) :: x !! finite point strictly above loc
        real(dp), intent(in) :: params(:) !! positive shape, loc, scale, reciprocal direction
        real(dp), intent(out) :: lz !! log standardized inverse-Gaussian coordinate
        real(dp), intent(out) :: la !! log absolute normal coordinate A
        real(dp), intent(out) :: lb !! log positive normal coordinate B
        real(dp), intent(out) :: sign_a !! sign of A, zero at the shape
        real(dp), intent(out) :: logphi !! log standard normal density at A
        real(dp) :: lv, lw, delta
        lz = params(4)*(log_distance(x,params(2))-log(params(3)))
        lv = 0.5_dp*lz-log(params(1))
        lw = -0.5_dp*lz
        lb = bounded_logadd(lv,lw)
        delta = lz-log(params(1))
        sign_a = sign(1.0_dp,delta)
        if (delta == 0.0_dp) then
            la = negative_infinity(0.0_dp)
            sign_a = 0.0_dp
        else
            la = max(lv,lw)+log(-expm1_safe(-abs(delta)))
        end if
        if (abs(delta) < 0.01_dp) call ig_center_coordinate(x,params,la,sign_a)
        logphi = negative_hazard(2.0_dp*la-scifort_log_two)-scifort_log_sqrt_two_pi
    end subroutine ig_coordinates

    pure subroutine ig_two_sum(a, b, high, low)
        implicit none
        real(dp), intent(in) :: a, b !! finite addends without overflow
        real(dp), intent(out) :: high !! rounded sum
        real(dp), intent(out) :: low !! exact rounding remainder
        real(dp) :: v
        high = a+b
        v = high-a
        low = (a-(high-v))+(b-v)
    end subroutine ig_two_sum

    pure subroutine ig_center_coordinate(x, params, la, sign_a)
        implicit none
        real(dp), intent(in) :: x !! finite point above loc, near the shape coordinate
        real(dp), intent(in) :: params(:) !! positive shape, loc, scale, reciprocal direction
        real(dp), intent(out) :: la !! log absolute normal coordinate A
        real(dp), intent(out) :: sign_a !! sign of A, zero for exact cancellation
        real(dp) :: dh, dl, sh, sl, fm, ph, pl, c, ah, al, bh, bl, h, l, residual
        real(dp), parameter :: splitter = real(2**((digits(1.0_dp)+1)/2),dp)+1.0_dp
        integer :: ed, es, em, common
        ! Form d=x-loc as a scaled two-component sum. Near the mode,
        ! subtracting logarithms would erase a small but significant A.
        ed = exponent(x)
        if (x == 0.0_dp) ed = minexponent(x)-digits(x)
        if (params(2) /= 0.0_dp) ed = max(ed,exponent(params(2)))
        call ig_two_sum(scale(x,-ed),-scale(params(2),-ed),dh,dl)
        sh = fraction(params(3))
        sl = 0.0_dp
        es = exponent(params(3))
        if (params(4) < 0.0_dp) then
            c = dh
            dh = sh
            sh = c
            sl = dl
            dl = 0.0_dp
            common = ed
            ed = es
            es = common
        end if
        fm = fraction(params(1))
        em = exponent(params(1))
        common = max(ed,em+es)
        ph = fm*sh
        ! Error-free product of normalized factors by a binary split.
        c = splitter*fm
        ah = c-(c-fm)
        al = fm-ah
        c = splitter*sh
        bh = c-(c-sh)
        bl = sh-bh
        pl = ((ah*bh-ph)+ah*bl+al*bh)+al*bl+fm*sl
        call ig_two_sum(scale(dh,ed-common),-scale(ph,em+es-common),h,l)
        residual = h+(l+scale(dl,ed-common)-scale(pl,em+es-common))
        sign_a = sign(1.0_dp,residual)
        la = negative_infinity(0.0_dp)
        if (residual == 0.0_dp) then
            sign_a = 0.0_dp
        else
            la = log(abs(residual))+real(common,dp)*scifort_log_two-log(params(1)) &
                -0.5_dp*(log_distance(x,params(2))+log(params(3)))
        end if
    end subroutine ig_center_coordinate

    pure function ig_logmills(lx) result(value)
        implicit none
        real(dp), intent(in) :: lx !! log positive normal argument
        real(dp) :: value, term, total, inverse_square, next_term
        integer :: n
        if (lx < log(12.0_dp)) then
            value = log(erfcx(exp(lx)/sqrt(2.0_dp))) &
                +0.5_dp*log(scifort_pi/2.0_dp)
        else
            inverse_square = exp(-2.0_dp*lx)
            term = 1.0_dp
            total = term
            do n = 1, 100
                next_term = -term*real(2*n-1,dp)*inverse_square
                if (abs(next_term) >= abs(term)) exit
                total = total+next_term
                term = next_term
                if (abs(term) <= epsilon(total)*abs(total)/4.0_dp) exit
            end do
            value = -lx+log(total)
        end if
    end function ig_logmills

    pure function ig_mills_difference(la, lb, ldelta) result(value)
        implicit none
        real(dp), intent(in) :: la !! log positive A >= 12
        real(dp), intent(in) :: lb !! log B > log A, possibly rounded equal
        real(dp), intent(in) :: ldelta !! log exact B-A, retained independently
        real(dp) :: value, ratio, power, polynomial, term, total, next_term, inverse_square
        integer :: n
        ratio = exp(min(0.0_dp,la-lb))
        inverse_square = exp(-2.0_dp*la)
        term = 1.0_dp
        total = term
        polynomial = 1.0_dp
        power = 1.0_dp
        ! Divide each inverse-power difference by (1-A/B) before summing.
        do n = 1, 100
            power = power*ratio
            polynomial = polynomial+power
            power = power*ratio
            polynomial = polynomial+power
            next_term = -term*real(2*n-1,dp)*inverse_square
            if (abs(next_term*polynomial) >= abs(term)*(2.0_dp*n+1.0_dp)) exit
            term = next_term
            total = total+term*polynomial
            if (abs(term*polynomial) <= epsilon(total)*abs(total)/4.0_dp) exit
        end do
        value = -la+ldelta-lb+log(total)
    end function ig_mills_difference

    pure function ig_complement(logp) result(value)
        implicit none
        real(dp), intent(in) :: logp !! log probability in [-infinity,0]
        real(dp) :: value
        if (logp < -scifort_log_two) then
            value = log1p_safe(-exp(logp))
        else
            value = log_difference(0.0_dp,logp)
        end if
    end function ig_complement

    pure subroutine ig_tails(x, params, lp, lq)
        implicit none
        real(dp), intent(in) :: x !! finite censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, loc, scale, reciprocal direction
        real(dp), intent(out) :: lp !! log lower probability
        real(dp), intent(out) :: lq !! log upper probability
        real(dp) :: lz, la, lb, sign_a, logphi, a, delta, terms(3), t(3), v, other, swap
        integer :: i
        lp = negative_infinity(0.0_dp)
        lq = 0.0_dp
        if (x <= params(2)) return
        call ig_coordinates(x,params,lz,la,lb,sign_a,logphi)
        if (sign_a > 0.0_dp .and. la >= log(12.0_dp)) then
            lq = logphi+ig_mills_difference(la,lb,scifort_log_two-0.5_dp*lz)
            lp = ig_complement(lq)
        else if (sign_a < 0.0_dp .and. la >= log(12.0_dp)) then
            lp = logphi+bounded_logadd(ig_logmills(la),ig_logmills(lb))
            lq = ig_complement(lp)
        else
            a = sign_a*exp(la)
            if (scifort_log_two-0.5_dp*lz < log(0.001_dp/(1.0_dp+abs(a)))) then
                ! R(A)-R(B) = integral_A^B (1-t*R(t)) dt.
                ! Three-point Gauss-Legendre avoids cancellation as B-A -> 0.
                delta = exp(scifort_log_two-0.5_dp*lz)
                t = [0.5_dp,0.5_dp-sqrt(3.0_dp/5.0_dp)/2.0_dp, &
                    0.5_dp+sqrt(3.0_dp/5.0_dp)/2.0_dp]
                do i = 1, 3
                    v = a+t(i)*delta
                    terms(i) = 1.0_dp-v*sqrt(scifort_pi/2.0_dp)*erfcx(v/sqrt(2.0_dp))
                end do
                lq = logphi+scifort_log_two-0.5_dp*lz &
                    +log((4.0_dp/9.0_dp)*terms(1)+(5.0_dp/18.0_dp)*(terms(2)+terms(3)))
                lp = ig_complement(lq)
            else
                other = logphi+ig_logmills(lb)
                if (a <= 0.0_dp) then
                    lp = bounded_logadd(normal_logcdf(a),other)
                    lq = ig_complement(lp)
                else
                    lq = log_difference(normal_logsf(a),other)
                    lp = ig_complement(lq)
                end if
            end if
        end if
        if (params(4) < 0.0_dp) then
            swap = lp
            lp = lq
            lq = swap
        end if
    end subroutine ig_tails

    pure function ig_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! positive shape, loc, scale, reciprocal direction
        real(dp) :: value, lz, la, lb, sign_a, logphi
        value = negative_infinity(0.0_dp)
        if (x <= params(2)) return
        call ig_coordinates(x,params,lz,la,lb,sign_a,logphi)
        value = logphi-(1.0_dp+0.5_dp*params(4))*params(4)*lz-log(params(3))
    end function ig_lpdf

    pure function ig_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, loc, scale, reciprocal direction
        real(dp) :: value, unused
        call ig_tails(x,params,value,unused)
    end function ig_lcdf

    pure function ig_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, loc, scale, reciprocal direction
        real(dp) :: value, unused
        call ig_tails(x,params,unused,value)
    end function ig_lsf

    pure function ig_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper endpoint greater than a
        real(dp), intent(in) :: params(:) !! positive shape, loc, scale, reciprocal direction
        real(dp) :: value, lza, lzb, slope, width, t(3), points(3), densities(3), ref, omitted
        real(dp) :: positive, negative, endpoint, rate
        integer :: i
        if (a > params(2)) then
            lza = log_distance(a,params(2))-log(params(3))
            lzb = log_distance(b,params(2))-log(params(3))
            if (params(4) > 0.0_dp) then
                slope = bounded_logadd(log(1.5_dp), &
                    bounded_logadd(lzb-2.0_dp*log(params(1)),-lza)-scifort_log_two)
            else
                slope = bounded_logadd(-scifort_log_two, &
                    bounded_logadd(lzb,-lza-2.0_dp*log(params(1)))-scifort_log_two)
            end if
            width = log_distance(b,a)
            if (width-log_distance(a,params(2)) < log(1.0e-6_dp)) then
                if (params(4) > 0.0_dp) then
                    positive = -lza-scifort_log_two
                    negative = bounded_logadd(log(1.5_dp), &
                        lza-2.0_dp*log(params(1))-scifort_log_two)
                else
                    positive = -lza-2.0_dp*log(params(1))-scifort_log_two
                    negative = bounded_logadd(-scifort_log_two,lza-scifort_log_two)
                end if
                ! In an extremely steep monotone interval, rounded log tails
                ! can coincide. An endpoint Laplace integral has relative
                ! curvature error O(1/abs(d log(f)/d log(z))).
                if (abs(positive-negative) > log(1.0e16_dp)) then
                    endpoint = a
                    if (positive > negative) endpoint = b
                    lzb = log_distance(endpoint,params(2))-log(params(3))
                    if (params(4) > 0.0_dp) then
                        positive = -lzb-scifort_log_two
                        negative = bounded_logadd(log(1.5_dp), &
                            lzb-2.0_dp*log(params(1))-scifort_log_two)
                    else
                        positive = -lzb-2.0_dp*log(params(1))-scifort_log_two
                        negative = bounded_logadd(-scifort_log_two,lzb-scifort_log_two)
                    end if
                    slope = max(positive,negative)
                    rate = slope-log_distance(endpoint,params(2))
                    if (slope > log(1.0e16_dp) .and. rate+width > log(1.0e8_dp)) then
                        value = ig_lpdf(endpoint,params)-rate &
                            +log_exponential_mass(rate+width)
                        return
                    end if
                end if
                ! Restore the variation bound used by the smooth integral.
                lzb = log_distance(b,params(2))-log(params(3))
                if (params(4) > 0.0_dp) then
                    slope = bounded_logadd(log(1.5_dp), &
                        bounded_logadd(lzb-2.0_dp*log(params(1)),-lza)-scifort_log_two)
                else
                    slope = bounded_logadd(-scifort_log_two, &
                        bounded_logadd(lzb,-lza-2.0_dp*log(params(1)))-scifort_log_two)
                end if
            end if
            if (width-log_distance(a,params(2))+slope < log(0.001_dp)) then
                t = [0.5_dp,0.5_dp-sqrt(3.0_dp/5.0_dp)/2.0_dp, &
                    0.5_dp+sqrt(3.0_dp/5.0_dp)/2.0_dp]
                points = (1.0_dp-t)*a+t*b
                do i = 1, 3
                    densities(i) = ig_lpdf(points(i),params)
                end do
                ref = maxval(densities)
                value = ref
                if (.not. ieee_is_finite(ref)) return
                value = ref+width+log((4.0_dp/9.0_dp)*exp(densities(1)-ref) &
                    +(5.0_dp/18.0_dp)*(exp(densities(2)-ref)+exp(densities(3)-ref)))
                return
            end if
        end if
        value = tail_interval(a,b,params,ig_lcdf,ig_lsf)
        if (value > -scifort_log_two) then
            omitted = bounded_logadd(ig_lcdf(a,params),ig_lsf(b,params))
            value = log1p_safe(-exp(omitted))
        end if
    end function ig_interval

    pure function invgauss_censored_loglikelihood(data, mu_shape, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: mu_shape !! positive finite inverse-Gaussian shape
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(4)
        value = quiet_nan(0.0_dp)
        params = [mu_shape,0.0_dp,1.0_dp,1.0_dp]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        if (.not. all(ieee_is_finite(params))) return
        if (params(1) <= 0.0_dp .or. params(3) <= 0.0_dp) return
        value = censored_loglikelihood(data,params,ig_lpdf,ig_lcdf,ig_lsf,ig_interval)
    end function invgauss_censored_loglikelihood

    pure function invgauss_censored_nnlf(data, mu_shape, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: mu_shape !! positive finite inverse-Gaussian shape
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -invgauss_censored_loglikelihood(data,mu_shape, loc,scale)
    end function invgauss_censored_nnlf

    pure function wald_censored_loglikelihood(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(4)
        value = quiet_nan(0.0_dp)
        params = [1.0_dp,0.0_dp,1.0_dp,1.0_dp]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        if (.not. all(ieee_is_finite(params))) return
        if (params(1) <= 0.0_dp .or. params(3) <= 0.0_dp) return
        value = censored_loglikelihood(data,params,ig_lpdf,ig_lcdf,ig_lsf,ig_interval)
    end function wald_censored_loglikelihood

    pure function wald_censored_nnlf(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -wald_censored_loglikelihood(data,loc,scale)
    end function wald_censored_nnlf

    pure function recipinvgauss_censored_loglikelihood(data, mu_shape, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: mu_shape !! positive finite inverse-Gaussian shape
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(4)
        value = quiet_nan(0.0_dp)
        params = [mu_shape,0.0_dp,1.0_dp,-1.0_dp]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        if (.not. all(ieee_is_finite(params))) return
        if (params(1) <= 0.0_dp .or. params(3) <= 0.0_dp) return
        value = censored_loglikelihood(data,params,ig_lpdf,ig_lcdf,ig_lsf,ig_interval)
    end function recipinvgauss_censored_loglikelihood

    pure function recipinvgauss_censored_nnlf(data, mu_shape, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: mu_shape !! positive finite inverse-Gaussian shape
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -recipinvgauss_censored_loglikelihood(data,mu_shape, loc,scale)
    end function recipinvgauss_censored_nnlf


    ! Beta coordinates for Jones-Faddy and its symmetric Student-t special case.
    ! Independent mathematical derivation: see CODE_PROVENANCE.md.
    pure subroutine skew_t_coordinates(x, params, u, v, lu, lv, lh)
        implicit none
        real(dp), intent(in) :: x !! finite observation or limit
        real(dp), intent(in) :: params(:) !! positive a, b, finite loc, positive scale
        real(dp), intent(out) :: u !! lower beta coordinate
        real(dp), intent(out) :: v !! independently formed upper beta coordinate
        real(dp), intent(out) :: lu !! log lower coordinate, retained through underflow
        real(dp), intent(out) :: lv !! log upper coordinate, retained through underflow
        real(dp), intent(out) :: lh !! log sqrt(a+b+z*z), with z=(x-loc)/scale
        real(dp) :: lz, ls, ratio, small, lsmall
        ls = log(params(1)+params(2))
        lh = 0.5_dp*ls
        u = 0.5_dp
        v = u
        lu = -scifort_log_two
        lv = lu
        if (x == params(3)) return
        if (x > params(3)) then
            lz = log_distance(x,params(3))-log(params(4))
        else
            lz = log_distance(params(3),x)-log(params(4))
        end if
        lh = 0.5_dp*ls+0.5_dp*log1pexp(2.0_dp*lz-ls)
        ratio = exp(min(0.0_dp,lz-lh))
        lsmall = ls-scifort_log_two-2.0_dp*lh-log1p_safe(ratio)
        small = exp(lsmall)
        if (x > params(3)) then
            lv = lsmall
            v = small
            u = 1.0_dp-small
            lu = log1p_safe(-small)
        else
            lu = lsmall
            u = small
            v = 1.0_dp-small
            lv = log1p_safe(-small)
        end if
    end subroutine skew_t_coordinates

    pure function skew_t_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! normal positive a, b, finite loc, positive scale
        real(dp) :: value, u, v, lu, lv, lh, kernel
        call skew_t_coordinates(x,params,u,v,lu,lv,lh)
        if (min(u,v) > 0.0_dp) then
            kernel = log_beta_kernel(params(1),params(2),u,v)
        else
            kernel = params(1)*lu+params(2)*lv-log_beta(params(1),params(2))
        end if
        value = kernel+0.5_dp*(lu+lv)+2.0_dp*scifort_log_two &
            -0.5_dp*log(params(1)+params(2))-log(params(4))
    end function skew_t_lpdf

    pure subroutine skew_t_tails(x, params, lp, lq)
        implicit none
        real(dp), intent(in) :: x !! finite censoring limit
        real(dp), intent(in) :: params(:) !! normal positive a, b, finite loc, positive scale
        real(dp), intent(out) :: lp !! log lower probability
        real(dp), intent(out) :: lq !! log upper probability
        real(dp) :: u, v, lu, lv, lh, p, q
        call skew_t_coordinates(x,params,u,v,lu,lv,lh)
        if (max(params(1),params(2)) <= epsilon(1.0_dp)/16384.0_dp) then
            lp = -log1pexp(log(params(1))-log(params(2)))
            lq = -log1pexp(log(params(2))-log(params(1)))
        else if (params(1) == params(2) .and. x == params(3)) then
            lp = -scifort_log_two
            lq = lp
        else if (u == 0.0_dp) then
            lp = betaprime_endpoint_tail(params(1),params(2),lu)
            lq = log_difference(0.0_dp,lp)
        else if (v == 0.0_dp) then
            lq = betaprime_endpoint_tail(params(2),params(1),lv)
            lp = log_difference(0.0_dp,lq)
        else
            call incomplete_beta_xy(params(1),params(2),u,v,p,q,lp,lq)
        end if
    end subroutine skew_t_tails

    pure function skew_t_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! normal positive a, b, finite loc, positive scale
        real(dp) :: value, unused
        call skew_t_tails(x,params,value,unused)
    end function skew_t_lcdf

    pure function skew_t_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! normal positive a, b, finite loc, positive scale
        real(dp) :: value, unused
        call skew_t_tails(x,params,unused,value)
    end function skew_t_lsf

    pure function skew_t_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper endpoint greater than a
        real(dp), intent(in) :: params(:) !! normal positive a, b, finite loc, positive scale
        real(dp) :: value, lh, u, v, lu, lv, near, logwidth, shape_sum
        real(dp) :: t(3), points(3), densities(3), reference, omitted
        shape_sum = params(1)+params(2)
        near = max(a,min(b,params(3)))
        call skew_t_coordinates(near,params,u,v,lu,lv,lh)
        logwidth = log_distance(b,a)
        ! |d log(f)/dz| <= (2*(a_shape+b_shape)+1)/sqrt(a_shape+b_shape+z*z).
        if (logwidth-log(params(4))-lh+log(2.0_dp*shape_sum+1.0_dp) <= log(0.001_dp)) then
            t = [0.5_dp,0.5_dp-sqrt(3.0_dp/5.0_dp)/2.0_dp, &
                0.5_dp+sqrt(3.0_dp/5.0_dp)/2.0_dp]
            if (a >= 0.0_dp .or. b <= 0.0_dp) then
                points = a+t*(b-a)
            else
                points = (1.0_dp-t)*a+t*b
            end if
            densities(1) = skew_t_lpdf(points(1),params)
            densities(2) = skew_t_lpdf(points(2),params)
            densities(3) = skew_t_lpdf(points(3),params)
            reference = maxval(densities)
            value = reference+logwidth+log((4.0_dp/9.0_dp)*exp(densities(1)-reference) &
                +(5.0_dp/18.0_dp)*(exp(densities(2)-reference)+exp(densities(3)-reference)))
        else if (max(params(1),params(2)) <= epsilon(1.0_dp)/16384.0_dp) then
            value = skew_t_small_interval(a,b,params)
        else
            value = tail_interval(a,b,params,skew_t_lcdf,skew_t_lsf)
            if (value > -scifort_log_two) then
                omitted = bounded_logadd(skew_t_lcdf(a,params),skew_t_lsf(b,params))
                value = log1p_safe(-exp(omitted))
            end if
        end if
    end function skew_t_interval

    pure function skew_t_logasinh(lr) result(value)
        implicit none
        real(dp), intent(in) :: lr !! log positive standardized distance divided by sqrt(a+b)
        real(dp) :: value
        if (lr < log(tiny(1.0_dp))) then
            value = lr
        else if (lr <= 0.0_dp) then
            value = log(asinh(exp(lr)))
        else
            value = log(lr+log1p_safe(sqrt(1.0_dp+exp(-2.0_dp*lr))))
        end if
    end function skew_t_logasinh

    pure function skew_t_small_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper endpoint greater than a
        real(dp), intent(in) :: params(:) !! tiny positive a, b, finite loc, positive scale
        real(dp) :: value, first, last, lnear, lfar, hnear, hfar, ls, logs, delta, ratio
        real(dp) :: u, v, lu, lv
        ls = log(params(1)+params(2))
        logs = log(params(4))
        if (a < params(3) .and. b > params(3)) then
            first = skew_t_logasinh(log_distance(params(3),a)-logs-0.5_dp*ls)
            last = skew_t_logasinh(log_distance(b,params(3))-logs-0.5_dp*ls)
            delta = bounded_logadd(first,last)
        else
            first = a
            last = b
            if (b <= params(3)) then
                first = b
                last = a
            end if
            call skew_t_coordinates(first,params,u,v,lu,lv,hnear)
            call skew_t_coordinates(last,params,u,v,lu,lv,hfar)
            if (first == params(3)) then
                lnear = negative_infinity(0.0_dp)
            else
                lnear = log_distance(max(first,params(3)),min(first,params(3)))-logs
            end if
            lfar = log_distance(max(last,params(3)),min(last,params(3)))-logs
            ratio = exp(min(0.0_dp,bounded_logadd(lnear,lfar)-bounded_logadd(hnear,hfar)))
            delta = burr_logsoftplus(log_distance(b,a)-logs-bounded_logadd(lnear,hnear) &
                +log1p_safe(ratio))
        end if
        value = scifort_log_two+log(params(1))+log(params(2))-ls+delta
    end function skew_t_small_interval

    pure function f_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper endpoint greater than a
        real(dp), intent(in) :: params(:) !! beta half-shapes, loc, scale, log dfn/dfd
        real(dp) :: value
        if (max(params(1),params(2)) > epsilon(1.0_dp)/16384.0_dp .or. a <= params(3)) then
            value = betaprime_interval(a,b,params)
        else
            value = log(params(1))+log(params(2))-log(params(1)+params(2)) &
                +burr_logsoftplus(log_distance(b,a)-log_distance(a,params(3)))
        end if
    end function f_interval

    pure function t_censored_loglikelihood(data, df, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: df !! degrees of freedom in [2*tiny,1e6], or +infinity
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(4)
        value = quiet_nan(0.0_dp)
        if (ieee_is_nan(df)) return
        if (df == positive_infinity(0.0_dp)) then
            value = normal_censored_loglikelihood(data,loc,scale)
            return
        end if
        if (.not. ieee_is_finite(df)) return
        if (df <= 0.0_dp) return
        params = [0.5_dp*df,0.5_dp*df,0.0_dp,1.0_dp]
        if (present(loc)) params(3) = loc
        if (present(scale)) params(4) = scale
        if (.not. all(ieee_is_finite(params))) return
        if (minval(params(1:2)) < tiny(1.0_dp)) return
        if (maxval(params(1:2)) > 1.0e6_dp) return
        if (params(1)+params(2) > 1.0e6_dp .or. params(4) <= 0.0_dp) return
        value = censored_loglikelihood(data,params,skew_t_lpdf,skew_t_lcdf, &
            skew_t_lsf,skew_t_interval)
    end function t_censored_loglikelihood

    pure function t_censored_nnlf(data, df, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: df !! degrees of freedom in [2*tiny,1e6], or +infinity
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -t_censored_loglikelihood(data,df,loc,scale)
    end function t_censored_nnlf

    pure function f_censored_loglikelihood(data, dfn, dfd, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: dfn !! numerator degrees >=2*tiny; dfn+dfd <=2e6
        real(dp), intent(in) :: dfd !! denominator degrees >=2*tiny; dfn+dfd <=2e6
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(5)
        value = quiet_nan(0.0_dp)
        if (.not. ieee_is_finite(dfn) .or. .not. ieee_is_finite(dfd)) return
        if (min(dfn,dfd) <= 0.0_dp) return
        params = [0.5_dp*dfn,0.5_dp*dfd,0.0_dp,1.0_dp,log(dfn)-log(dfd)]
        if (present(loc)) params(3) = loc
        if (present(scale)) params(4) = scale
        if (.not. all(ieee_is_finite(params))) return
        if (minval(params(1:2)) < tiny(1.0_dp)) return
        if (maxval(params(1:2)) > 1.0e6_dp) return
        if (params(1)+params(2) > 1.0e6_dp .or. params(4) <= 0.0_dp) return
        value = censored_loglikelihood(data,params,betaprime_lpdf,betaprime_lcdf, &
            betaprime_lsf,f_interval)
    end function f_censored_loglikelihood

    pure function f_censored_nnlf(data, dfn, dfd, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: dfn !! numerator degrees >=2*tiny; dfn+dfd <=2e6
        real(dp), intent(in) :: dfd !! denominator degrees >=2*tiny; dfn+dfd <=2e6
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -f_censored_loglikelihood(data,dfn, dfd,loc,scale)
    end function f_censored_nnlf

    pure function jf_skew_t_censored_loglikelihood(data, a, b, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: a !! first beta shape >=tiny; a+b <=1e6
        real(dp), intent(in) :: b !! second beta shape >=tiny; a+b <=1e6
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(4)
        value = quiet_nan(0.0_dp)
        params = [a,b,0.0_dp,1.0_dp]
        if (present(loc)) params(3) = loc
        if (present(scale)) params(4) = scale
        if (.not. all(ieee_is_finite(params))) return
        if (minval(params(1:2)) < tiny(1.0_dp)) return
        if (maxval(params(1:2)) > 1.0e6_dp) return
        if (params(1)+params(2) > 1.0e6_dp .or. params(4) <= 0.0_dp) return
        value = censored_loglikelihood(data,params,skew_t_lpdf,skew_t_lcdf, &
            skew_t_lsf,skew_t_interval)
    end function jf_skew_t_censored_loglikelihood

    pure function jf_skew_t_censored_nnlf(data, a, b, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: a !! first beta shape >=tiny; a+b <=1e6
        real(dp), intent(in) :: b !! second beta shape >=tiny; a+b <=1e6
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -jf_skew_t_censored_loglikelihood(data,a, b,loc,scale)
    end function jf_skew_t_censored_nnlf


    ! Gamma-power laws: Y = rate*((x-loc)/scale)**power ~ Gamma(a,1).
    ! Derived from the public density specifications; see CODE_PROVENANCE.md.
    pure function gamma_power_logy(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite point strictly above location
        real(dp), intent(in) :: params(:) !! gamma shape, loc, scale, log rate, power
        real(dp) :: value
        value = params(5)*lognormal_logz(x,params(1:3))+params(4)
    end function gamma_power_logy

    pure function gamma_power_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! gamma shape, loc, scale, log rate, power
        real(dp) :: value, logy
        value = negative_infinity(0.0_dp)
        if (x < params(2)) return
        if (x == params(2)) then
            if (params(1) < 1.0_dp/params(5)) value = positive_infinity(0.0_dp)
            if (params(1) == 1.0_dp/params(5)) value = log(params(5)) &
                +params(1)*params(4)-log_gamma(params(1))-log(params(3))
            return
        end if
        logy = gamma_power_logy(x,params)
        value = loggamma_lpdf(logy,[params(1),0.0_dp,1.0_dp]) &
            +log(params(5))-log_distance(x,params(2))
    end function gamma_power_lpdf

    pure function gamma_power_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! gamma shape, loc, scale, log rate, power
        real(dp) :: value
        value = negative_infinity(0.0_dp)
        if (x <= params(2)) return
        if (params(1) <= epsilon(1.0_dp)/4096.0_dp) then
            value = log_difference(0.0_dp,gamma_power_lsf(x,params))
            return
        end if
        value = loggamma_lcdf(gamma_power_logy(x,params),[params(1),0.0_dp,1.0_dp])
    end function gamma_power_lcdf

    pure function gamma_power_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! gamma shape, loc, scale, log rate, power
        real(dp) :: value, logy
        real(dp), parameter :: euler_gamma = 0.577215664901532860606512090082402431_dp
        value = 0.0_dp
        if (x <= params(2)) return
        logy = gamma_power_logy(x,params)
        if (params(1) <= epsilon(1.0_dp)/4096.0_dp .and. logy <= log(1.5_dp)) then
            ! Q(a,Y)=a*E1(Y)*(1+O(a*(1+abs(log(Y))))). Keep log(a).
            if (logy < log(tiny(1.0_dp))) then
                value = log(params(1))+log(-logy-euler_gamma)
            else
                value = log(params(1))+log(exp1(exp(logy)))
            end if
        else
            value = loggamma_lsf(logy,[params(1),0.0_dp,1.0_dp])
        end if
    end function gamma_power_lsf

    pure function gamma_power_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper endpoint greater than a
        real(dp), intent(in) :: params(:) !! gamma shape, loc, scale, log rate, power
        real(dp) :: value, logwidth, logq, logy, q, t(3), d(3), r(3)
        real(dp) :: first(3), second(3), correction, omitted
        value = negative_infinity(0.0_dp)
        if (b <= params(2)) return
        if (a <= params(2)) then
            value = gamma_power_lcdf(b,params)
            return
        end if
        logwidth = log_distance(b,a)
        logq = logwidth-log_distance(a,params(2))
        logy = gamma_power_logy(a,params)
        ! Bound relative width and both terms of the log-density variation.
        if (logq+max(0.0_dp,log(params(5))+log(params(1)), &
                log(params(5))+logy) <= log(0.001_dp)) then
            q = exp(logq)
            t = [0.5_dp,0.5_dp-sqrt(3.0_dp/5.0_dp)/2.0_dp, &
                0.5_dp+sqrt(3.0_dp/5.0_dp)/2.0_dp]
            d = log1p_safe(t*q)
            first = t
            second = params(5)*t
            if (q > 0.0_dp) then
                first = d/q
                second = expm1_safe(params(5)*d)/q
            end if
            r = exp(log(params(1))+log(params(5))+logq)*first-d &
                -exp(logy+logq)*second
            correction = (4.0_dp/9.0_dp)*expm1_safe(r(1)) &
                +(5.0_dp/18.0_dp)*(expm1_safe(r(2))+expm1_safe(r(3)))
            value = gamma_power_lpdf(a,params)+logwidth+log1p_safe(correction)
        else
            value = tail_interval(a,b,params,gamma_power_lcdf,gamma_power_lsf)
            if (value > -scifort_log_two) then
                omitted = bounded_logadd(gamma_power_lcdf(a,params),gamma_power_lsf(b,params))
                value = log1p_safe(-exp(omitted))
            end if
        end if
    end function gamma_power_interval

    pure function chi2_censored_loglikelihood(data, df, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: df !! degrees of freedom in (0,2e6], with positive fl(df/2)
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(5)
        value = quiet_nan(0.0_dp)
        if (.not. ieee_is_finite(df)) return
        if (df <= 0.0_dp) return
        params = [0.5_dp*df,0.0_dp,1.0_dp,-scifort_log_two,1.0_dp]
        if (params(1) > 1.0e6_dp) return
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        if (.not. all(ieee_is_finite(params))) return
        if (params(1) <= 0.0_dp .or. params(3) <= 0.0_dp) return
        value = censored_loglikelihood(data,params,gamma_power_lpdf,gamma_power_lcdf, &
            gamma_power_lsf,gamma_power_interval)
    end function chi2_censored_loglikelihood

    pure function chi2_censored_nnlf(data, df, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: df !! degrees of freedom in (0,2e6], with positive fl(df/2)
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -chi2_censored_loglikelihood(data,df,loc,scale)
    end function chi2_censored_nnlf

    pure function chi_censored_loglikelihood(data, df, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: df !! degrees of freedom in (0,2e6], with positive fl(df/2)
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(5)
        value = quiet_nan(0.0_dp)
        if (.not. ieee_is_finite(df)) return
        if (df <= 0.0_dp) return
        params = [0.5_dp*df,0.0_dp,1.0_dp,-scifort_log_two,2.0_dp]
        if (params(1) > 1.0e6_dp) return
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        if (.not. all(ieee_is_finite(params))) return
        if (params(1) <= 0.0_dp .or. params(3) <= 0.0_dp) return
        value = censored_loglikelihood(data,params,gamma_power_lpdf,gamma_power_lcdf, &
            gamma_power_lsf,gamma_power_interval)
    end function chi_censored_loglikelihood

    pure function chi_censored_nnlf(data, df, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: df !! degrees of freedom in (0,2e6], with positive fl(df/2)
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -chi_censored_loglikelihood(data,df,loc,scale)
    end function chi_censored_nnlf

    pure function nakagami_censored_loglikelihood(data, nu, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: nu !! Nakagami shape in (0,1e6]
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(5)
        value = quiet_nan(0.0_dp)
        if (.not. ieee_is_finite(nu)) return
        if (nu <= 0.0_dp) return
        params = [nu,0.0_dp,1.0_dp,log(nu),2.0_dp]
        if (params(1) > 1.0e6_dp) return
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        if (.not. all(ieee_is_finite(params))) return
        if (params(1) <= 0.0_dp .or. params(3) <= 0.0_dp) return
        value = censored_loglikelihood(data,params,gamma_power_lpdf,gamma_power_lcdf, &
            gamma_power_lsf,gamma_power_interval)
    end function nakagami_censored_loglikelihood

    pure function nakagami_censored_nnlf(data, nu, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: nu !! Nakagami shape in (0,1e6]
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -nakagami_censored_loglikelihood(data,nu,loc,scale)
    end function nakagami_censored_nnlf


    pure function rdist_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! equal half-shapes, rounded lower edge, full width
        real(dp) :: value, left, right, z, w, lz, lw
        if (params(1) > epsilon(1.0_dp)/4096.0_dp) then
            value = beta_lpdf(x,params)
            return
        end if
        value = negative_infinity(0.0_dp)
        call beta_distances(x,params,left,right)
        if (left < 0.0_dp .or. right < 0.0_dp) return
        value = positive_infinity(0.0_dp)
        if (left == 0.0_dp .or. right == 0.0_dp) return
        call beta_coordinates(x,params,left,right,z,w,lz,lw)
        value = log(params(1))-scifort_log_two-log(params(4))-lz-lw
    end function rdist_lpdf

    pure subroutine rdist_tails(x, params, lp, lq)
        implicit none
        real(dp), intent(in) :: x !! finite censoring limit
        real(dp), intent(in) :: params(:) !! equal half-shapes, rounded lower edge, full width
        real(dp), intent(out) :: lp !! log lower probability
        real(dp), intent(out) :: lq !! log upper probability
        real(dp) :: left, right
        call beta_distances(x,params,left,right)
        lp = negative_infinity(0.0_dp)
        lq = 0.0_dp
        if (left <= 0.0_dp) return
        lp = 0.0_dp
        lq = negative_infinity(0.0_dp)
        if (right <= 0.0_dp) return
        if (params(1) <= epsilon(1.0_dp)/4096.0_dp .or. left == right) then
            lp = -scifort_log_two
            lq = lp
        else
            call beta_logtails(x,params,lp,lq)
        end if
    end subroutine rdist_tails

    pure function rdist_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! equal half-shapes, rounded lower edge, full width
        real(dp) :: value, unused
        call rdist_tails(x,params,value,unused)
    end function rdist_lcdf

    pure function rdist_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! equal half-shapes, rounded lower edge, full width
        real(dp) :: value, unused
        call rdist_tails(x,params,unused,value)
    end function rdist_lsf

    pure function rdist_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper endpoint greater than a
        real(dp), intent(in) :: params(:) !! equal half-shapes, rounded lower edge, full width
        real(dp) :: value, left, right, lower, upper, span
        if (params(1) > epsilon(1.0_dp)/4096.0_dp) then
            value = beta_interval(a,b,params)
            return
        end if
        value = negative_infinity(0.0_dp)
        call beta_distances(b,params,left,right)
        if (left <= 0.0_dp) return
        if (right <= 0.0_dp) then
            value = rdist_lsf(a,params)
            return
        end if
        upper = log(right)
        call beta_distances(a,params,left,right)
        if (right <= 0.0_dp) return
        if (left <= 0.0_dp) then
            value = rdist_lcdf(b,params)
            return
        end if
        lower = log_distance(a,params(3))
        span = log_distance(b,a)
        value = log(params(1))-scifort_log_two &
            +bounded_logadd(burr_logsoftplus(span-lower),burr_logsoftplus(span-upper))
    end function rdist_interval

    pure function generalized_normal_norm(a) result(value)
        implicit none
        real(dp), intent(in) :: a !! positive finite reciprocal power
        real(dp) :: value, product
        if (a < stirling_threshold) then
            value = -log_gamma_one_plus(a)
        else
            product = bounded_cutoff_product(a,log(a)-1.0_dp)
            value = negative_infinity(0.0_dp)
            if (.not. ieee_is_finite(product)) return
            value = -product-0.5_dp*log(a)-scifort_log_sqrt_two_pi-stirling_remainder(a)
        end if
    end function generalized_normal_norm

    pure function generalized_normal_power(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite point at or above location
        real(dp), intent(in) :: params(:) !! power, location, scale, reciprocal power, log normalizer, overflow flag
        real(dp) :: value, lz
        value = negative_infinity(0.0_dp)
        if (x == params(2)) return
        call kappa3_power(x,params,lz,value)
    end function generalized_normal_power

    pure subroutine generalized_normal_gamma(ly, params, lp, lq)
        implicit none
        real(dp), intent(in) :: ly !! log powered coordinate, possibly infinite
        real(dp), intent(in) :: params(:) !! power, location, scale, reciprocal power, log normalizer, overflow flag
        real(dp), intent(out) :: lp !! log probability between location and positive observation
        real(dp), intent(out) :: lq !! log positive half-law survival probability
        real(dp) :: y
        lp = negative_infinity(0.0_dp)
        lq = 0.0_dp
        if (ly == negative_infinity(0.0_dp)) return
        if (params(6) /= 0.0_dp) return
        if (ly > log(huge(1.0_dp))) then
            lp = 0.0_dp
            lq = negative_infinity(0.0_dp)
        else if (ly < log(tiny(1.0_dp))) then
            lp = gamma_small_kernel(params(4),ly)
            lq = log_difference(0.0_dp,lp)
        else
            y = exp(ly)
            if (params(4) > 1.0_dp) then
                if (params(4) > huge(1.0_dp)/max(1.0_dp,log(params(4)))) then
                    ! Here y is near one: the lower-gamma series correction is O(1/a).
                    lp = params(5)+params(4)*ly-y
                    lq = log_difference(0.0_dp,lp)
                    return
                end if
            end if
            lp = log_gammainc(params(4),y)
            lq = log_gammaincc(params(4),y)
        end if
    end subroutine generalized_normal_gamma

    pure subroutine generalized_normal_point(x, params, lp, lq)
        implicit none
        real(dp), intent(in) :: x !! finite point at or above the half-law location
        real(dp), intent(in) :: params(:) !! power, location, scale, reciprocal power, normalizer, flag
        real(dp), intent(out) :: lp !! log positive half-law CDF
        real(dp), intent(out) :: lq !! log positive half-law SF
        real(dp) :: lz, ly
        lp = negative_infinity(0.0_dp)
        lq = 0.0_dp
        if (x == params(2) .or. params(6) /= 0.0_dp) return
        call kappa3_power(x,params,lz,ly)
        if (ly < log(tiny(1.0_dp))) then
            ! (z**beta)**(1/beta) = z, even if beta*log(z) overflows.
            lp = lz+params(5)
            lq = log_difference(0.0_dp,lp)
        else
            call generalized_normal_gamma(ly,params,lp,lq)
        end if
    end subroutine generalized_normal_point

    pure function halfgennorm_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! power, location, scale, reciprocal power, log normalizer, overflow flag
        real(dp) :: value
        value = negative_infinity(0.0_dp)
        if (x < params(2)) return
        if (params(6) /= 0.0_dp) return
        value = params(5)-log(params(3))+negative_hazard(generalized_normal_power(x,params))
    end function halfgennorm_lpdf

    pure function gennorm_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! power, location, scale, reciprocal power, log normalizer, overflow flag
        real(dp) :: value, reflected(6)
        reflected = params
        reflected(2) = -params(2)
        if (x < params(2)) then
            value = halfgennorm_lpdf(-x,reflected)-scifort_log_two
        else
            value = halfgennorm_lpdf(x,params)-scifort_log_two
        end if
    end function gennorm_lpdf

    pure subroutine generalized_normal_tails(x, params, half, lp, lq)
        implicit none
        real(dp), intent(in) :: x !! finite censoring limit
        real(dp), intent(in) :: params(:) !! power, location, scale, reciprocal power, log normalizer, overflow flag
        logical, intent(in) :: half !! restrict support to the positive half
        real(dp), intent(out) :: lp !! log lower-tail probability
        real(dp), intent(out) :: lq !! log upper-tail probability
        real(dp) :: first, second, reflected(6)
        if (half) then
            lp = negative_infinity(0.0_dp)
            lq = 0.0_dp
            if (x <= params(2)) return
            call generalized_normal_point(x,params,lp,lq)
        else
            reflected = params
            reflected(2) = -params(2)
            if (x < params(2)) then
                call generalized_normal_point(-x,reflected,first,second)
                lp = second-scifort_log_two
                lq = log1p_safe(-exp(lp))
            else
                call generalized_normal_point(x,params,first,second)
                lq = second-scifort_log_two
                lp = log1p_safe(-exp(lq))
            end if
        end if
    end subroutine generalized_normal_tails

    pure function halfgennorm_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper endpoint greater than a
        real(dp), intent(in) :: params(:) !! power, location, scale, reciprocal power, log normalizer, overflow flag
        real(dp) :: value, lo, la, lb, lr, ld, lp, lq, rp, rq, t(3), r(3), correction
        integer :: j
        value = negative_infinity(0.0_dp)
        if (b <= params(2)) return
        lo = max(a,params(2))
        lb = generalized_normal_power(b,params)
        if (lb < log(epsilon(1.0_dp))-2.0_dp) then
            value = params(5)-log(params(3))+log_distance(b,lo)
            if (params(6) /= 0.0_dp) value = negative_infinity(0.0_dp)
            return
        end if
        call generalized_normal_point(b,params,rp,rq)
        if (lo == params(2)) then
            value = rp
            return
        end if
        la = generalized_normal_power(lo,params)
        call generalized_normal_point(lo,params,lp,lq)
        if (la > log(huge(1.0_dp))) return
        lr = log_distance(b,lo)-log_distance(lo,params(2))
        ld = la+log_expm1_from_log(log(params(1))+burr_logsoftplus(lr))
        if (lr <= log(0.001_dp) .and. ld <= log(0.001_dp) &
                .and. log(params(1))+burr_logsoftplus(lr) <= log(0.001_dp)) then
            t = [0.5_dp,0.5_dp-sqrt(3.0_dp/5.0_dp)/2.0_dp, &
                0.5_dp+sqrt(3.0_dp/5.0_dp)/2.0_dp]
            do j = 1, 3
                ld = la+log_expm1_from_log(log(params(1))+burr_logsoftplus(log(t(j))+lr))
                r(j) = negative_hazard(ld)
            end do
            correction = (4.0_dp/9.0_dp)*expm1_safe(r(1)) &
                +(5.0_dp/18.0_dp)*(expm1_safe(r(2))+expm1_safe(r(3)))
            value = halfgennorm_lpdf(lo,params)+log_distance(b,lo)+log1p_safe(correction)
        else if (rp <= lq) then
            value = log_difference(rp,lp)
        else
            value = log_difference(lq,rq)
        end if
        if (value > -scifort_log_two) value = log1p_safe(-exp(bounded_logadd(lp,rq)))
    end function halfgennorm_interval

    pure function gennorm_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper endpoint greater than a
        real(dp), intent(in) :: params(:) !! power, location, scale, reciprocal power, log normalizer, overflow flag
        real(dp) :: value, reflected(6), left, right, unused
        reflected = params
        reflected(2) = -params(2)
        if (b <= params(2)) then
            value = halfgennorm_interval(-b,-a,reflected)-scifort_log_two
        else if (a >= params(2)) then
            value = halfgennorm_interval(a,b,params)-scifort_log_two
        else
            call generalized_normal_point(-a,reflected,left,unused)
            call generalized_normal_point(b,params,right,unused)
            value = bounded_logadd(left,right)-scifort_log_two
        end if
        if (value > -scifort_log_two) then
            call generalized_normal_tails(a,params,.false.,left,unused)
            call generalized_normal_tails(b,params,.false.,unused,right)
            value = log1p_safe(-exp(bounded_logadd(left,right)))
        end if
    end function gennorm_interval

    pure function gennorm_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite censoring limit
        real(dp), intent(in) :: params(:) !! power, location, scale, reciprocal power, log normalizer, overflow flag
        real(dp) :: value, unused
        call generalized_normal_tails(x,params,.false.,value,unused)
    end function gennorm_lcdf

    pure function gennorm_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite censoring limit
        real(dp), intent(in) :: params(:) !! power, location, scale, reciprocal power, log normalizer, overflow flag
        real(dp) :: value, unused
        call generalized_normal_tails(x,params,.false.,unused,value)
    end function gennorm_lsf

    pure function halfgennorm_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite censoring limit
        real(dp), intent(in) :: params(:) !! power, location, scale, reciprocal power, log normalizer, overflow flag
        real(dp) :: value, unused
        call generalized_normal_tails(x,params,.true.,value,unused)
    end function halfgennorm_lcdf

    pure function halfgennorm_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite censoring limit
        real(dp), intent(in) :: params(:) !! power, location, scale, reciprocal power, log normalizer, overflow flag
        real(dp) :: value, unused
        call generalized_normal_tails(x,params,.true.,unused,value)
    end function halfgennorm_lsf


    pure function gennorm_censored_loglikelihood(data, beta, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: beta !! positive finite generalized-normal power
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(6)
        params = [beta,0.0_dp,1.0_dp,0.0_dp,0.0_dp,0.0_dp]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (beta <= 0.0_dp) return
        if (params(3) <= 0.0_dp) return
        params(4) = location_scale_z(1.0_dp,[0.0_dp,beta])
        if (.not. ieee_is_finite(params(4))) return
        params(5) = generalized_normal_norm(params(4))
        if (.not. ieee_is_finite(params(5))) then
            params(5) = 0.0_dp
            params(6) = 1.0_dp
        end if
        value = censored_loglikelihood(data,params,gennorm_lpdf,gennorm_lcdf, &
            gennorm_lsf,gennorm_interval)
    end function gennorm_censored_loglikelihood

    pure function gennorm_censored_nnlf(data, beta, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: beta !! positive finite generalized-normal power
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -gennorm_censored_loglikelihood(data, beta, loc, scale)
    end function gennorm_censored_nnlf

    pure function halfgennorm_censored_loglikelihood(data, beta, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: beta !! positive finite generalized-normal power
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(6)
        params = [beta,0.0_dp,1.0_dp,0.0_dp,0.0_dp,0.0_dp]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (beta <= 0.0_dp) return
        if (params(3) <= 0.0_dp) return
        params(4) = location_scale_z(1.0_dp,[0.0_dp,beta])
        if (.not. ieee_is_finite(params(4))) return
        params(5) = generalized_normal_norm(params(4))
        if (.not. ieee_is_finite(params(5))) then
            params(5) = 0.0_dp
            params(6) = 1.0_dp
        end if
        value = censored_loglikelihood(data,params,halfgennorm_lpdf,halfgennorm_lcdf, &
            halfgennorm_lsf,halfgennorm_interval)
    end function halfgennorm_censored_loglikelihood

    pure function halfgennorm_censored_nnlf(data, beta, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: beta !! positive finite generalized-normal power
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -halfgennorm_censored_loglikelihood(data, beta, loc, scale)
    end function halfgennorm_censored_nnlf

    pure function rdist_censored_loglikelihood(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! positive finite symmetric-beta shape, at most 1e6
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(4), center, width
        params = [c,0.0_dp,1.0_dp,0.0_dp]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (c <= 0.0_dp) return
        if (c > 1.0e6_dp) return
        if (params(3) <= 0.0_dp) return
        if (0.5_dp*c <= 0.0_dp) return
        if (params(3) > 0.5_dp*huge(1.0_dp)) return
        center = params(2)
        width = 2.0_dp*params(3)
        if (center < 0.0_dp) then
            if (params(3) > huge(1.0_dp)+center) return
        end if
        params = [0.5_dp*c,0.5_dp*c,center-params(3),width]
        if (.not. beta_parameters_valid(params)) return
        if (params(3) >= center .or. params(3)+width <= center) return
        value = censored_loglikelihood(data,params,rdist_lpdf,rdist_lcdf, &
            rdist_lsf,rdist_interval)
    end function rdist_censored_loglikelihood

    pure function rdist_censored_nnlf(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! positive finite symmetric-beta shape, at most 1e6
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -rdist_censored_loglikelihood(data, c, loc, scale)
    end function rdist_censored_nnlf


    pure subroutine folded_coordinates(x, params, minus, plus, lm, lp, lz)
        implicit none
        real(dp), intent(in) :: x !! finite point at or above location
        real(dp), intent(in) :: params(:) !! nonnegative offset, location, positive scale
        real(dp), intent(out) :: minus !! standardized distance minus offset
        real(dp), intent(out) :: plus !! standardized distance plus offset
        real(dp), intent(out) :: lm !! log absolute minus coordinate
        real(dp), intent(out) :: lp !! log plus coordinate
        real(dp), intent(out) :: lz !! log standardized distance
        real(dp) :: z, half, lc
        lz = negative_infinity(0.0_dp)
        if (x > params(2)) lz = log_distance(x,params(2))-log(params(3))
        lc = bounded_log(params(1))
        lp = bounded_logadd(lz,lc)
        z = location_scale_z(x,params(2:3))
        plus = finite_sum(z,params(1))
        minus = z-params(1)
        if (.not. ieee_is_finite(z)) then
            half = location_scale_z(0.5_dp*x,[0.5_dp*params(2),params(3)])
            minus = half-0.5_dp*params(1)
            minus = finite_sum(minus,minus)
        end if
        if (ieee_is_finite(minus)) then
            lm = bounded_log(abs(minus))
        else
            lm = lz+log1p_safe(-exp(lc-lz))
        end if
    end subroutine folded_coordinates

    pure function folded_normal_mass(a, b, lw) result(value)
        implicit none
        real(dp), intent(in) :: a !! standardized lower endpoint, possibly infinite
        real(dp), intent(in) :: b !! standardized upper endpoint, possibly coincident after rounding
        real(dp), intent(in) :: lw !! log exact standardized positive interval width
        real(dp) :: value, omitted
        value = normal_small_interval(a,lw)
        if (.not. ieee_is_nan(value)) return
        if (a >= 0.0_dp) then
            value = log_difference(normal_logsf(a),normal_logsf(b))
        else
            value = log_difference(normal_logcdf(b),normal_logcdf(a))
            if (b > 0.0_dp) then
                omitted = exp(normal_logcdf(a))+exp(normal_logsf(b))
                if (omitted <= 0.5_dp) value = log1p_safe(-omitted)
            end if
        end if
    end function folded_normal_mass

    pure function folded_cauchy_mass(a, b, la, lb, lw) result(value)
        implicit none
        real(dp), intent(in) :: a !! standardized lower endpoint
        real(dp), intent(in) :: b !! standardized upper endpoint
        real(dp), intent(in) :: la !! log absolute exact lower coordinate
        real(dp), intent(in) :: lb !! log absolute exact upper coordinate
        real(dp), intent(in) :: lw !! log exact positive standardized interval width
        real(dp) :: value, omitted
        if (a >= 0.0_dp .or. b <= 0.0_dp) then
            value = cauchy_logatan_exp(lw-log1pexp(la+lb))-log(scifort_pi)
        else if (min(la,lb) >= 0.0_dp) then
            omitted = exp(cauchy_log_small_tail(la))+exp(cauchy_log_small_tail(lb))
            value = log1p_safe(-omitted)
        else
            value = bounded_logadd(cauchy_logatan_exp(la),cauchy_logatan_exp(lb)) &
                -log(scifort_pi)
        end if
    end function folded_cauchy_mass

    pure function foldnorm_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! nonnegative offset, location, positive scale
        real(dp) :: value, a, b, la, lb, lz
        value = negative_infinity(0.0_dp)
        if (x < params(2)) return
        call folded_coordinates(x,params,a,b,la,lb,lz)
        value = bounded_logadd(finite_normal_lpdf(a),finite_normal_lpdf(b))-log(params(3))
    end function foldnorm_lpdf

    pure function foldcauchy_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! nonnegative offset, location, positive scale
        real(dp) :: value, a, b, la, lb, lz
        value = negative_infinity(0.0_dp)
        if (x < params(2)) return
        call folded_coordinates(x,params,a,b,la,lb,lz)
        value = bounded_logadd(-log1pexp(2.0_dp*la),-log1pexp(2.0_dp*lb)) &
            -log(scifort_pi)-log(params(3))
    end function foldcauchy_lpdf

    pure subroutine folded_tails(x, params, normal, lcdf, lsf)
        implicit none
        real(dp), intent(in) :: x !! finite censoring limit
        real(dp), intent(in) :: params(:) !! nonnegative offset, location, positive scale
        logical, intent(in) :: normal !! normal rather than Cauchy base law
        real(dp), intent(out) :: lcdf !! log lower-tail probability
        real(dp), intent(out) :: lsf !! log upper-tail probability
        real(dp) :: a, b, la, lb, lz, first
        lcdf = negative_infinity(0.0_dp)
        lsf = 0.0_dp
        if (x <= params(2)) return
        call folded_coordinates(x,params,a,b,la,lb,lz)
        if (normal) then
            lcdf = folded_normal_mass(-b,a,scifort_log_two+lz)
            lsf = bounded_logadd(normal_logsf(a),normal_logsf(b))
        else
            lcdf = folded_cauchy_mass(-b,a,lb,la,scifort_log_two+lz)
            first = cauchy_log_small_tail(la)
            if (a < 0.0_dp) first = log1p_safe(-exp(first))
            lsf = bounded_logadd(first,cauchy_log_small_tail(lb))
        end if
        if (lcdf <= lsf) then
            lcdf = min(lcdf,0.0_dp)
            lsf = log1p_safe(-exp(lcdf))
        else
            lsf = min(lsf,0.0_dp)
            lcdf = log1p_safe(-exp(lsf))
        end if
    end subroutine folded_tails

    pure function folded_interval(a, b, params, normal) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper interval endpoint greater than a
        real(dp), intent(in) :: params(:) !! nonnegative offset, location, positive scale
        logical, intent(in) :: normal !! normal rather than Cauchy base law
        real(dp) :: value, lo, am, ap, bm, bp, lam, lap, lbm, lbp, lz, lw, omitted, unused
        value = negative_infinity(0.0_dp)
        if (b <= params(2)) return
        lo = max(a,params(2))
        call folded_coordinates(lo,params,am,ap,lam,lap,lz)
        call folded_coordinates(b,params,bm,bp,lbm,lbp,lz)
        lw = log_distance(b,lo)-log(params(3))
        if (normal) then
            value = bounded_logadd(folded_normal_mass(am,bm,lw),folded_normal_mass(ap,bp,lw))
        else
            value = bounded_logadd(folded_cauchy_mass(am,bm,lam,lbm,lw), &
                folded_cauchy_mass(ap,bp,lap,lbp,lw))
        end if
        if (value > -scifort_log_two) then
            call folded_tails(lo,params,normal,omitted,unused)
            call folded_tails(b,params,normal,unused,lz)
            value = log1p_safe(-(exp(omitted)+exp(lz)))
        end if
    end function folded_interval

    pure function genpareto_logh(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! interior point strictly above location
        real(dp), intent(in) :: params(:) !! shape, location, scale, bounded width if negative
        real(dp) :: value, lu, lw, nh
        if (params(1) < 0.0_dp) then
            call genhalflogistic_coordinates(x,[-params(1),params(2),params(4)],lu,lw,nh,value)
        else
            value = log_distance(x,params(2))-log(params(3))
            if (params(1) > 0.0_dp) value = &
                burr_logsoftplus(log(params(1))+value)-log(params(1))
        end if
    end function genpareto_logh

    pure function genpareto_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! shape, location, scale, bounded width if negative
        real(dp) :: value, left, right, lu, lw, nh, lh
        value = negative_infinity(0.0_dp)
        if (x < params(2)) return
        if (params(1) < 0.0_dp) then
            call beta_distances(x,[1.0_dp,1.0_dp,params(2),params(4)],left,right)
            if (right < 0.0_dp) return
            if (right == 0.0_dp) then
                if (params(1) == -1.0_dp) value = -log(params(4))
                if (params(1) < -1.0_dp) value = positive_infinity(0.0_dp)
                return
            end if
            value = -log(-params(1))-log(params(4))
            if (left == 0.0_dp) return
            call genhalflogistic_coordinates(x,[-params(1),params(2),params(4)],lu,lw,nh,lh)
            if (params(1) <= -0.5_dp) then
                value = value+((1.0_dp+params(1))/(-params(1)))*lw
            else
                value = value+nh-lw
            end if
        else
            value = -log(params(3))
            if (x == params(2)) return
            value = value+negative_hazard(genpareto_logh(x,params))
            if (params(1) > 0.0_dp) value = value &
                -log1pexp(log(params(1))+log_distance(x,params(2))-log(params(3)))
        end if
    end function genpareto_lpdf

    pure function genpareto_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! shape, location, scale, bounded width if negative
        real(dp) :: value, left, right
        value = negative_infinity(0.0_dp)
        if (x <= params(2)) return
        if (params(1) < 0.0_dp) then
            call beta_distances(x,[1.0_dp,1.0_dp,params(2),params(4)],left,right)
            value = 0.0_dp
            if (right <= 0.0_dp) return
        end if
        value = log_hazard_mass(genpareto_logh(x,params))
    end function genpareto_lcdf

    pure function genpareto_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! shape, location, scale, bounded width if negative
        real(dp) :: value, left, right
        value = 0.0_dp
        if (x <= params(2)) return
        if (params(1) < 0.0_dp) then
            call beta_distances(x,[1.0_dp,1.0_dp,params(2),params(4)],left,right)
            value = negative_infinity(0.0_dp)
            if (right <= 0.0_dp) return
        end if
        value = negative_hazard(genpareto_logh(x,params))
    end function genpareto_lsf

    pure function genpareto_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper endpoint greater than a
        real(dp), intent(in) :: params(:) !! shape, location, scale, bounded width if negative
        real(dp) :: value, lo, left, right, delta, denom
        value = negative_infinity(0.0_dp)
        if (b <= params(2)) return
        lo = max(a,params(2))
        if (params(1) < 0.0_dp) then
            call beta_distances(b,[1.0_dp,1.0_dp,params(2),params(4)],left,right)
            if (right <= 0.0_dp) then
                value = genpareto_lsf(lo,params)
                return
            end if
            delta = burr_logsoftplus(log_distance(b,lo)-log(right))-log(-params(1))
        else if (params(1) == 0.0_dp) then
            delta = log_distance(b,lo)-log(params(3))
        else
            denom = log(params(3))
            if (lo > params(2)) denom = bounded_logadd(denom, &
                log(params(1))+log_distance(lo,params(2)))
            delta = burr_logsoftplus(log(params(1))+log_distance(b,lo)-denom)-log(params(1))
        end if
        value = genpareto_lsf(lo,params)+log_hazard_mass(delta)
    end function genpareto_interval

    pure function foldnorm_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite censoring limit
        real(dp), intent(in) :: params(:) !! nonnegative offset, location, positive scale
        real(dp) :: value, unused
        call folded_tails(x,params,.true.,value,unused)
    end function foldnorm_lcdf

    pure function foldnorm_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite censoring limit
        real(dp), intent(in) :: params(:) !! nonnegative offset, location, positive scale
        real(dp) :: value, unused
        call folded_tails(x,params,.true.,unused,value)
    end function foldnorm_lsf

    pure function foldnorm_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper endpoint greater than a
        real(dp), intent(in) :: params(:) !! nonnegative offset, location, positive scale
        real(dp) :: value
        value = folded_interval(a,b,params,.true.)
    end function foldnorm_interval

    pure function foldcauchy_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite censoring limit
        real(dp), intent(in) :: params(:) !! nonnegative offset, location, positive scale
        real(dp) :: value, unused
        call folded_tails(x,params,.false.,value,unused)
    end function foldcauchy_lcdf

    pure function foldcauchy_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite censoring limit
        real(dp), intent(in) :: params(:) !! nonnegative offset, location, positive scale
        real(dp) :: value, unused
        call folded_tails(x,params,.false.,unused,value)
    end function foldcauchy_lsf

    pure function foldcauchy_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper endpoint greater than a
        real(dp), intent(in) :: params(:) !! nonnegative offset, location, positive scale
        real(dp) :: value
        value = folded_interval(a,b,params,.false.)
    end function foldcauchy_interval


    pure function foldnorm_censored_loglikelihood(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! nonnegative finite folding offset
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(3)
        params = [c,0.0_dp,1.0_dp]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (c < 0.0_dp) return
        if (params(3) <= 0.0_dp) return
        value = censored_loglikelihood(data,params,foldnorm_lpdf,foldnorm_lcdf, &
            foldnorm_lsf,foldnorm_interval)
    end function foldnorm_censored_loglikelihood

    pure function foldnorm_censored_nnlf(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! nonnegative finite folding offset
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -foldnorm_censored_loglikelihood(data, c, loc, scale)
    end function foldnorm_censored_nnlf

    pure function foldcauchy_censored_loglikelihood(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! nonnegative finite folding offset
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(3)
        params = [c,0.0_dp,1.0_dp]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (c < 0.0_dp) return
        if (params(3) <= 0.0_dp) return
        value = censored_loglikelihood(data,params,foldcauchy_lpdf,foldcauchy_lcdf, &
            foldcauchy_lsf,foldcauchy_interval)
    end function foldcauchy_censored_loglikelihood

    pure function foldcauchy_censored_nnlf(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! nonnegative finite folding offset
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -foldcauchy_censored_loglikelihood(data, c, loc, scale)
    end function foldcauchy_censored_nnlf

    pure function genpareto_censored_loglikelihood(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! finite generalized Pareto shape, including zero
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(4)
        params = [c,0.0_dp,1.0_dp,0.0_dp]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (params(3) <= 0.0_dp) return
        if (c < 0.0_dp) then
            if (-c < 1.0_dp) then
                if (params(3) > huge(1.0_dp)*(-c)) return
            end if
            params(4) = params(3)/(-c)
            if (.not. beta_parameters_valid([1.0_dp,1.0_dp,params(2),params(4)])) return
            if (params(2)+params(4) <= params(2)) return
        end if
        value = censored_loglikelihood(data,params,genpareto_lpdf,genpareto_lcdf, &
            genpareto_lsf,genpareto_interval)
    end function genpareto_censored_loglikelihood

    pure function genpareto_censored_nnlf(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! finite generalized Pareto shape, including zero
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -genpareto_censored_loglikelihood(data, c, loc, scale)
    end function genpareto_censored_nnlf


    pure function bounded_cutoff_product(a, b) result(value)
        implicit none
        real(dp), intent(in) :: a !! positive finite cutoff multiplier
        real(dp), intent(in) :: b !! positive finite scale
        real(dp) :: value, half
        value = quiet_nan(0.0_dp)
        if (a > 1.0_dp) then
            if (b > huge(1.0_dp)/a) return
            if (b <= (0.5_dp*huge(1.0_dp))/a) then
                value = a*b
                return
            end if
            ! The quotient guard can round upward. Test the half product before doubling.
            half = (0.5_dp*a)*b
            if (half > 0.5_dp*huge(1.0_dp)) return
            value = 2.0_dp*half
        else
            value = a*b
        end if
    end function bounded_cutoff_product

    pure function bounded_power_valid(params) result(valid)
        implicit none
        real(dp), intent(in) :: params(:) !! power, location, positive physical cutoffs
        logical :: valid
        valid = .false.
        if (params(3) <= 0.0_dp .or. params(4) <= params(3)) return
        if (.not. beta_parameters_valid([1.0_dp,1.0_dp,params(2),params(4)])) return
        if (params(2)+params(3) <= params(2)) return
        if (params(2)+params(4) <= params(2)+params(3)) return
        valid = .true.
    end function bounded_power_valid

    pure subroutine bounded_power_gaps(x, params, left, right)
        implicit none
        real(dp), intent(in) :: x !! observation
        real(dp), intent(in) :: params(:) !! power, location, positive physical cutoffs
        real(dp), intent(out) :: left !! signed gap from lower cutoff
        real(dp), intent(out) :: right !! signed gap to upper cutoff
        real(dp) :: unused
        call beta_distances(x,[1.0_dp,1.0_dp,params(2),params(4)],unused,right)
        call beta_distances(x,[1.0_dp,1.0_dp,params(2),params(3)],unused,left)
        left = -left
    end subroutine bounded_power_gaps

    pure function bounded_power_range(params) result(value)
        implicit none
        real(dp), intent(in) :: params(:) !! power, location, positive physical cutoffs
        real(dp) :: value
        value = burr_logsoftplus(log(params(4)-params(3))-log(params(3)))
    end function bounded_power_range

    pure function bounded_power_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! observation
        real(dp), intent(in) :: params(:) !! power, location, positive physical cutoffs
        real(dp) :: value, left, right, lr, lv, anchor, lk, lb, beta, v2, term, total
        integer :: j
        value = negative_infinity(0.0_dp)
        call bounded_power_gaps(x,params,left,right)
        if (left < 0.0_dp .or. right < 0.0_dp) return
        lr = bounded_power_range(params)
        value = -log_distance(x,params(2))-lr
        if (params(1) == 0.0_dp) return
        lk = log(abs(params(1)))
        if (params(1) > 0.0_dp) then
            lv = burr_logsoftplus(bounded_log(left)-log(params(3)))
        else
            lv = burr_logsoftplus(bounded_log(right)-log_distance(x,params(2)))
        end if
        anchor = negative_hazard(lk+lv)
        lb = lk+lr
        if (lb <= 0.0_dp) then
            beta = exp(lb)
            v2 = (0.5_dp*beta)**2
            term = 1.0_dp
            total = 0.0_dp
            do j = 1, 10
                term = term*v2/real((2*j)*(2*j+1),dp)
                total = total+term
            end do
            value = value+0.5_dp*beta-log1p_safe(total)+anchor
        else
            value = lk-log_distance(x,params(2))+anchor-log_hazard_mass(lb)
        end if
    end function bounded_power_lpdf

    pure subroutine bounded_power_tails(x, params, lp, lq)
        implicit none
        real(dp), intent(in) :: x !! observation
        real(dp), intent(in) :: params(:) !! power, location, positive physical cutoffs
        real(dp), intent(out) :: lp !! log lower probability
        real(dp), intent(out) :: lq !! log upper probability
        real(dp) :: left, right, lv, lw, lr, lk, norm
        call bounded_power_gaps(x,params,left,right)
        lp = negative_infinity(0.0_dp)
        lq = 0.0_dp
        if (left <= 0.0_dp) return
        lp = 0.0_dp
        lq = negative_infinity(0.0_dp)
        if (right <= 0.0_dp) return
        lr = bounded_power_range(params)
        lv = burr_logsoftplus(log(left)-log(params(3)))
        lw = burr_logsoftplus(log(right)-log_distance(x,params(2)))
        if (params(1) == 0.0_dp) then
            lp = lv-lr
            lq = lw-lr
        else
            lk = log(abs(params(1)))
            norm = log_hazard_mass(lk+lr)
            lp = log_hazard_mass(lk+lv)-norm
            lq = log_hazard_mass(lk+lw)-norm
            if (params(1) > 0.0_dp) then
                lq = negative_hazard(lk+lv)+lq
            else
                lp = negative_hazard(lk+lw)+lp
            end if
        end if
        if (lp <= lq) then
            lp = min(lp,0.0_dp)
            lq = log1p_safe(-exp(lp))
        else
            lq = min(lq,0.0_dp)
            lp = log1p_safe(-exp(lq))
        end if
    end subroutine bounded_power_tails

    pure function bounded_power_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! observation
        real(dp), intent(in) :: params(:) !! power, location, positive physical cutoffs
        real(dp) :: value, unused
        call bounded_power_tails(x,params,value,unused)
    end function bounded_power_lcdf

    pure function bounded_power_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! observation
        real(dp), intent(in) :: params(:) !! power, location, positive physical cutoffs
        real(dp) :: value, unused
        call bounded_power_tails(x,params,unused,value)
    end function bounded_power_lsf

    pure function bounded_power_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper endpoint greater than a
        real(dp), intent(in) :: params(:) !! power, location, positive physical cutoffs
        real(dp) :: value, left, right, omitted, span, lk, anchor
        value = negative_infinity(0.0_dp)
        call bounded_power_gaps(b,params,left,right)
        if (left <= 0.0_dp) return
        if (right <= 0.0_dp) then
            value = bounded_power_lsf(a,params)
            return
        end if
        call bounded_power_gaps(a,params,left,right)
        if (right <= 0.0_dp) return
        if (left <= 0.0_dp) then
            value = bounded_power_lcdf(b,params)
            return
        end if
        omitted = exp(bounded_power_lcdf(a,params))+exp(bounded_power_lsf(b,params))
        if (omitted <= 0.5_dp) then
            value = log1p_safe(-omitted)
            return
        end if
        span = burr_logsoftplus(log_distance(b,a)-log_distance(a,params(2)))
        if (params(1) == 0.0_dp) then
            value = span-bounded_power_range(params)
        else
            lk = log(abs(params(1)))
            if (params(1) > 0.0_dp) then
                anchor = burr_logsoftplus(log(left)-log(params(3)))
            else
                call bounded_power_gaps(b,params,left,right)
                anchor = burr_logsoftplus(log(right)-log_distance(b,params(2)))
            end if
            value = negative_hazard(lk+anchor)+log_hazard_mass(lk+span) &
                -log_hazard_mass(lk+bounded_power_range(params))
        end if
    end function bounded_power_interval

    pure subroutine genhalflogistic_coordinates(x, params, lu, lw, nh, lh)
        implicit none
        real(dp), intent(in) :: x !! interior observation
        real(dp), intent(in) :: params(:) !! positive shape, location, rounded physical width
        real(dp), intent(out) :: lu !! log lower coordinate
        real(dp), intent(out) :: lw !! log upper coordinate
        real(dp), intent(out) :: nh !! negative transformed hazard
        real(dp), intent(out) :: lh !! log transformed hazard
        real(dp) :: left, right, u, w, bp(4)
        bp = [1.0_dp,1.0_dp,params(2),params(3)]
        call beta_distances(x,bp,left,right)
        call beta_coordinates(x,bp,left,right,u,w,lu,lw)
        if (lu < log(epsilon(1.0_dp))-2.0_dp) then
            lh = lu-log(params(1))
        else
            lh = log(-lw)-log(params(1))
        end if
        nh = negative_hazard(lh)
    end subroutine genhalflogistic_coordinates

    pure function genhalflogistic_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! observation
        real(dp), intent(in) :: params(:) !! positive shape, location, rounded physical width
        real(dp) :: value, left, right, lu, lw, nh, lh
        call beta_distances(x,[1.0_dp,1.0_dp,params(2),params(3)],left,right)
        value = negative_infinity(0.0_dp)
        if (left < 0.0_dp .or. right < 0.0_dp) return
        if (left == 0.0_dp) then
            value = -scifort_log_two-log(params(1))-log(params(3))
        else if (right == 0.0_dp) then
            if (params(1) == 1.0_dp) value = scifort_log_two-log(params(3))
            if (params(1) > 1.0_dp) value = positive_infinity(0.0_dp)
        else
            call genhalflogistic_coordinates(x,params,lu,lw,nh,lh)
            value = scifort_log_two-log(params(1))-log(params(3))-lw+nh &
                -2.0_dp*log1pexp(nh)
        end if
    end function genhalflogistic_lpdf

    pure subroutine genhalflogistic_tails(x, params, lp, lq)
        implicit none
        real(dp), intent(in) :: x !! observation
        real(dp), intent(in) :: params(:) !! positive shape, location, rounded physical width
        real(dp), intent(out) :: lp !! log lower probability
        real(dp), intent(out) :: lq !! log upper probability
        real(dp) :: left, right, lu, lw, nh, lh
        call beta_distances(x,[1.0_dp,1.0_dp,params(2),params(3)],left,right)
        lp = negative_infinity(0.0_dp)
        lq = 0.0_dp
        if (left <= 0.0_dp) return
        lp = 0.0_dp
        lq = negative_infinity(0.0_dp)
        if (right <= 0.0_dp) return
        call genhalflogistic_coordinates(x,params,lu,lw,nh,lh)
        lp = log_hazard_mass(lh)-log1pexp(nh)
        lq = scifort_log_two+nh-log1pexp(nh)
        if (lp <= lq) then
            lq = log1p_safe(-exp(lp))
        else
            lp = log1p_safe(-exp(lq))
        end if
    end subroutine genhalflogistic_tails

    pure function genhalflogistic_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! observation
        real(dp), intent(in) :: params(:) !! positive shape, location, rounded physical width
        real(dp) :: value, unused
        call genhalflogistic_tails(x,params,value,unused)
    end function genhalflogistic_lcdf

    pure function genhalflogistic_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! observation
        real(dp), intent(in) :: params(:) !! positive shape, location, rounded physical width
        real(dp) :: value, unused
        call genhalflogistic_tails(x,params,unused,value)
    end function genhalflogistic_lsf

    pure function genhalflogistic_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper endpoint greater than a
        real(dp), intent(in) :: params(:) !! positive shape, location, rounded physical width
        real(dp) :: value, left, right, omitted, delta, lu, lw, nha, nhb, lh
        value = negative_infinity(0.0_dp)
        call beta_distances(b,[1.0_dp,1.0_dp,params(2),params(3)],left,right)
        if (left <= 0.0_dp) return
        if (right <= 0.0_dp) then
            value = genhalflogistic_lsf(a,params)
            return
        end if
        delta = burr_logsoftplus(log_distance(b,a)-log(right))-log(params(1))
        call beta_distances(a,[1.0_dp,1.0_dp,params(2),params(3)],left,right)
        if (right <= 0.0_dp) return
        if (left <= 0.0_dp) then
            value = genhalflogistic_lcdf(b,params)
            return
        end if
        omitted = exp(genhalflogistic_lcdf(a,params))+exp(genhalflogistic_lsf(b,params))
        if (omitted <= 0.5_dp) then
            value = log1p_safe(-omitted)
            return
        end if
        call genhalflogistic_coordinates(a,params,lu,lw,nha,lh)
        call genhalflogistic_coordinates(b,params,lu,lw,nhb,lh)
        value = scifort_log_two+nha+log_hazard_mass(delta) &
            -log1pexp(nha)-log1pexp(nhb)
    end function genhalflogistic_interval


    pure function loguniform_censored_loglikelihood(data, a, b, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: a !! positive finite lower shape, less than b
        real(dp), intent(in) :: b !! finite upper shape, greater than a
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(4)
        params = [a,b,0.0_dp,1.0_dp]
        if (present(loc)) params(3) = loc
        if (present(scale)) params(4) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (min(a,b) <= 0.0_dp) return
        if (params(4) <= 0.0_dp) return
        if (a >= b) return
        if (b > 1.0_dp) then
            if (params(4) > huge(1.0_dp)/b) return
        end if
        params = [0.0_dp,params(3),bounded_cutoff_product(a,params(4)), &
            bounded_cutoff_product(b,params(4))]
        if (.not. all(ieee_is_finite(params))) return
        if (.not. bounded_power_valid(params)) return
        value = censored_loglikelihood(data,params,bounded_power_lpdf,bounded_power_lcdf, &
            bounded_power_lsf,bounded_power_interval)
    end function loguniform_censored_loglikelihood

    pure function loguniform_censored_nnlf(data, a, b, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: a !! positive finite lower shape, less than b
        real(dp), intent(in) :: b !! finite upper shape, greater than a
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -loguniform_censored_loglikelihood(data, a, b, loc, scale)
    end function loguniform_censored_nnlf

    pure function truncpareto_censored_loglikelihood(data, b, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: b !! finite nonzero power, either sign
        real(dp), intent(in) :: c !! finite upper cutoff, greater than one
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(4)
        params = [b,c,0.0_dp,1.0_dp]
        if (present(loc)) params(3) = loc
        if (present(scale)) params(4) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (b == 0.0_dp .or. c <= 1.0_dp) return
        if (params(4) <= 0.0_dp) return
        if (c > 1.0_dp) then
            if (params(4) > huge(1.0_dp)/c) return
        end if
        params = [b,params(3),params(4),bounded_cutoff_product(c,params(4))]
        if (.not. all(ieee_is_finite(params))) return
        if (.not. bounded_power_valid(params)) return
        value = censored_loglikelihood(data,params,bounded_power_lpdf,bounded_power_lcdf, &
            bounded_power_lsf,bounded_power_interval)
    end function truncpareto_censored_loglikelihood

    pure function truncpareto_censored_nnlf(data, b, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: b !! finite nonzero power, either sign
        real(dp), intent(in) :: c !! finite upper cutoff, greater than one
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -truncpareto_censored_loglikelihood(data, b, c, loc, scale)
    end function truncpareto_censored_nnlf

    pure function genhalflogistic_censored_loglikelihood(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! positive finite shape
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(3)
        params = [c,0.0_dp,1.0_dp]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (c <= 0.0_dp) return
        if (params(3) <= 0.0_dp) return
        if (c < 1.0_dp) then
            if (params(3) > huge(1.0_dp)*c) return
        end if
        params(3) = params(3)/c
        if (.not. beta_parameters_valid([1.0_dp,1.0_dp,params(2),params(3)])) return
        if (params(2)+params(3) <= params(2)) return
        value = censored_loglikelihood(data,params,genhalflogistic_lpdf,genhalflogistic_lcdf, &
            genhalflogistic_lsf,genhalflogistic_interval)
    end function genhalflogistic_censored_loglikelihood

    pure function genhalflogistic_censored_nnlf(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! positive finite shape
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -genhalflogistic_censored_loglikelihood(data, c, loc, scale)
    end function genhalflogistic_censored_nnlf


    ! Independently integrated public density specifications; see CODE_PROVENANCE.md.
    pure subroutine truncexpon_gaps(x, params, left, right)
        implicit none
        real(dp), intent(in) :: x !! finite observation or censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, location, rate scale, rounded width
        real(dp), intent(out) :: left !! signed distance above lower support edge
        real(dp), intent(out) :: right !! compensated distance below upper support edge
        call beta_distances(x,[1.0_dp,1.0_dp,params(2),params(4)],left,right)
    end subroutine truncexpon_gaps

    pure function truncexpon_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! positive shape, location, rate scale, rounded width
        real(dp) :: value, left, right, cutoff, logcutoff, term, series, v2, gap, unused
        integer :: k
        call truncexpon_gaps(x,params,left,right)
        value = negative_infinity(0.0_dp)
        if (left < 0.0_dp .or. right < 0.0_dp) return
        logcutoff = log(params(4))-log(params(3))
        if (logcutoff <= 0.0_dp) then
            ! exp(-z)/(1-exp(-cutoff)) = exp(cutoff/2-z)
            ! / (cutoff*sinh(cutoff/2)/(cutoff/2)). Retain the nearly uniform log density.
            cutoff = exp(logcutoff)
            v2 = (0.5_dp*cutoff)**2
            term = 1.0_dp
            series = 0.0_dp
            do k = 1, 10
                term = term*v2/real((2*k)*(2*k+1),dp)
                series = series+term
            end do
            if (params(4) >= 2.0_dp*tiny(1.0_dp)) then
                call beta_distances(x,[1.0_dp,1.0_dp,params(2),0.5_dp*params(4)],unused,gap)
                value = gap/params(3)
            else
                value = 0.5_dp*cutoff-location_scale_z(x,params(2:3))
            end if
            value = value-log(params(4))-log1p_safe(series)
            return
        end if
        value = -location_scale_z(x,params(2:3))-log(params(3)) &
            -log_hazard_mass(logcutoff)
    end function truncexpon_lpdf

    pure subroutine truncexpon_tails(x, params, logp, logq)
        implicit none
        real(dp), intent(in) :: x !! finite censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, location, rate scale, rounded width
        real(dp), intent(out) :: logp !! log lower-tail probability
        real(dp), intent(out) :: logq !! log upper-tail probability
        real(dp) :: left, right, normalization
        logp = negative_infinity(0.0_dp)
        logq = 0.0_dp
        call truncexpon_gaps(x,params,left,right)
        if (left <= 0.0_dp) return
        if (right <= 0.0_dp) then
            logp = 0.0_dp
            logq = negative_infinity(0.0_dp)
            return
        end if
        normalization = log_hazard_mass(log(params(4))-log(params(3)))
        logp = log_hazard_mass(log_distance(x,params(2))-log(params(3)))-normalization
        logq = -location_scale_z(x,params(2:3)) &
            +log_hazard_mass(log(right)-log(params(3)))-normalization
        if (logp <= logq) then
            logp = min(0.0_dp,logp)
            logq = log1p_safe(-exp(logp))
        else
            logq = min(0.0_dp,logq)
            logp = log1p_safe(-exp(logq))
        end if
    end subroutine truncexpon_tails

    pure function truncexpon_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, location, rate scale, rounded width
        real(dp) :: value, unused
        call truncexpon_tails(x,params,value,unused)
    end function truncexpon_lcdf

    pure function truncexpon_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, location, rate scale, rounded width
        real(dp) :: value, unused
        call truncexpon_tails(x,params,unused,value)
    end function truncexpon_lsf

    pure function truncexpon_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper endpoint greater than a
        real(dp), intent(in) :: params(:) !! positive shape, location, rate scale, rounded width
        real(dp) :: value, left, right, omitted
        call truncexpon_gaps(b,params,left,right)
        value = negative_infinity(0.0_dp)
        if (left <= 0.0_dp) return
        if (right <= 0.0_dp) then
            value = truncexpon_lsf(a,params)
            return
        end if
        call truncexpon_gaps(a,params,left,right)
        if (right <= 0.0_dp) return
        if (left <= 0.0_dp) then
            value = truncexpon_lcdf(b,params)
            return
        end if
        omitted = exp(truncexpon_lcdf(a,params))+exp(truncexpon_lsf(b,params))
        if (omitted <= 0.5_dp) then
            value = log1p_safe(-omitted)
        else
            value = -location_scale_z(a,params(2:3)) &
                +log_hazard_mass(log_distance(b,a)-log(params(3))) &
                -log_hazard_mass(log(params(4))-log(params(3)))
        end if
    end function truncexpon_interval

    pure function johnsonsb_affine(a, b, t) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite first Johnson shape
        real(dp), intent(in) :: b !! positive finite second Johnson shape
        real(dp), intent(in) :: t !! finite logit coordinate
        real(dp) :: value, halfsum
        if (abs(t) > 1.0_dp) then
            if (b > huge(1.0_dp)/abs(t)) then
                value = sign(positive_infinity(0.0_dp),t)
                if (0.5_dp*b > huge(1.0_dp)/abs(t)) return
                halfsum = finite_sum(0.5_dp*a,(0.5_dp*b)*t)
                if (abs(halfsum) > 0.5_dp*huge(1.0_dp)) return
                value = 2.0_dp*halfsum
                return
            end if
        end if
        value = finite_sum(a,b*t)
    end function johnsonsb_affine

    pure subroutine johnsonsb_coordinates(x, params, q, lz, lw)
        implicit none
        real(dp), intent(in) :: x !! finite point strictly inside support
        real(dp), intent(in) :: params(:) !! finite a, positive b, lower edge, positive width
        real(dp), intent(out) :: q !! transformed normal coordinate, possibly infinite
        real(dp), intent(out) :: lz !! log standardized distance above lower edge
        real(dp), intent(out) :: lw !! log standardized distance below upper edge
        real(dp) :: left, right, z, w, t, gap, unused, v, product
        call beta_distances(x,params,left,right)
        call beta_coordinates(x,params,left,right,z,w,lz,lw)
        t = lz-lw
        if (abs(t) < 0.1_dp .and. params(4) < 2.0_dp*tiny(1.0_dp)) then
            ! Subnormal widths can have an unrepresentable half-step midpoint.
            v = (left-right)/params(4)
            t = log1p_safe(v)-log1p_safe(-v)
            q = johnsonsb_affine(params(1),params(2),t)
            return
        end if
        if (abs(t) < 0.1_dp .and. 0.5_dp*params(4) > 0.0_dp) then
            call beta_distances(x,[1.0_dp,1.0_dp,params(3),0.5_dp*params(4)],unused,gap)
            if (gap == 0.0_dp) then
                q = params(1)
                return
            end if
            v = (-2.0_dp*gap)/params(4)
            if (abs(v) < sqrt(epsilon(1.0_dp))) then
                product = signed_exp(log(params(2))+log(abs(gap))-log(params(4)) &
                    +2.0_dp*scifort_log_two,-sign(1.0_dp,gap))
                q = finite_sum(params(1),product)
                return
            end if
            t = log1p_safe(v)-log1p_safe(-v)
        end if
        q = johnsonsb_affine(params(1),params(2),t)
    end subroutine johnsonsb_coordinates

    pure function johnsonsb_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! finite a, positive b, lower edge, positive width
        real(dp) :: value, left, right, q, lz, lw
        value = negative_infinity(0.0_dp)
        call beta_distances(x,params,left,right)
        if (min(left,right) <= 0.0_dp) return
        call johnsonsb_coordinates(x,params,q,lz,lw)
        value = finite_normal_lpdf(q)+log(params(2))-log(params(4))-lz-lw
    end function johnsonsb_lpdf

    pure subroutine johnsonsb_tails(x, params, logp, logq)
        implicit none
        real(dp), intent(in) :: x !! finite censoring limit
        real(dp), intent(in) :: params(:) !! finite a, positive b, lower edge, positive width
        real(dp), intent(out) :: logp !! log lower-tail probability
        real(dp), intent(out) :: logq !! log upper-tail probability
        real(dp) :: left, right, q, lz, lw
        logp = negative_infinity(0.0_dp)
        logq = 0.0_dp
        call beta_distances(x,params,left,right)
        if (left <= 0.0_dp) return
        if (right <= 0.0_dp) then
            logp = 0.0_dp
            logq = negative_infinity(0.0_dp)
            return
        end if
        call johnsonsb_coordinates(x,params,q,lz,lw)
        logp = normal_logcdf(q)
        logq = normal_logsf(q)
    end subroutine johnsonsb_tails

    pure function johnsonsb_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! finite a, positive b, lower edge, positive width
        real(dp) :: value, unused
        call johnsonsb_tails(x,params,value,unused)
    end function johnsonsb_lcdf

    pure function johnsonsb_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! finite a, positive b, lower edge, positive width
        real(dp) :: value, unused
        call johnsonsb_tails(x,params,unused,value)
    end function johnsonsb_lsf

    pure function johnsonsb_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper endpoint greater than a
        real(dp), intent(in) :: params(:) !! finite a, positive b, lower edge, positive width
        real(dp) :: value, left, right, q, lz, lw, span, first, second, omitted
        call beta_distances(b,params,left,right)
        value = negative_infinity(0.0_dp)
        if (left <= 0.0_dp) return
        if (right <= 0.0_dp) then
            value = johnsonsb_lsf(a,params)
            return
        end if
        second = log(right)
        call beta_distances(a,params,left,right)
        if (right <= 0.0_dp) return
        if (left <= 0.0_dp) then
            value = johnsonsb_lcdf(b,params)
            return
        end if
        span = log_distance(b,a)
        first = burr_logsoftplus(span-log_distance(a,params(3)))
        second = burr_logsoftplus(span-second)
        call johnsonsb_coordinates(a,params,q,lz,lw)
        value = normal_small_interval(q,log(params(2))+bounded_logadd(first,second))
        if (ieee_is_nan(value)) then
            value = tail_interval(a,b,params,johnsonsb_lcdf,johnsonsb_lsf)
            if (value > -scifort_log_two) then
                omitted = bounded_logadd(johnsonsb_lcdf(a,params),johnsonsb_lsf(b,params))
                value = log1p_safe(-exp(omitted))
            end if
        end if
    end function johnsonsb_interval

    pure subroutine kappa3_power(x, params, lz, power)
        implicit none
        real(dp), intent(in) :: x !! finite point strictly above lower support edge
        real(dp), intent(in) :: params(:) !! positive shape, location, positive scale
        real(dp), intent(out) :: lz !! log standardized distance, possibly rounded to zero
        real(dp), intent(out) :: power !! shape times log standardized distance
        real(dp) :: z, gap, unused, logr
        z = location_scale_z(x,params(2:3))
        if (z >= 0.5_dp .and. z <= 2.0_dp) then
            if (params(2) > 0.0_dp) then
                if (params(3) > huge(1.0_dp)-params(2)) then
                    ! Reorder scale-(x-loc) when loc+scale is unrepresentable.
                    call beta_distances(-params(2),[1.0_dp,1.0_dp,-x,params(3)],unused,gap)
                else
                    call beta_distances(x,[1.0_dp,1.0_dp,params(2),params(3)],unused,gap)
                end if
            else
                call beta_distances(x,[1.0_dp,1.0_dp,params(2),params(3)],unused,gap)
            end if
            lz = 0.0_dp
            power = 0.0_dp
            if (gap == 0.0_dp) return
            logr = log(abs(gap))-log(params(3))
            if (logr < log(epsilon(1.0_dp))-2.0_dp) then
                lz = -gap/params(3)
                power = signed_exp(log(params(1))+logr,-sign(1.0_dp,gap))
                return
            end if
            lz = log1p_safe(-gap/params(3))
        else
            lz = log_distance(x,params(2))-log(params(3))
        end if
        power = sign(positive_infinity(0.0_dp),lz)
        if (params(1) > 1.0_dp) then
            if (abs(lz) > huge(1.0_dp)/params(1)) return
        end if
        power = params(1)*lz
    end subroutine kappa3_power

    pure subroutine kappa3_hazard(x, params, q, logh, nh)
        implicit none
        real(dp), intent(in) :: x !! finite point strictly above lower support edge
        real(dp), intent(in) :: params(:) !! positive shape, location, positive scale
        real(dp), intent(out) :: q !! log(a)-a*log(z), possibly infinite
        real(dp), intent(out) :: logh !! log reversed cumulative hazard
        real(dp), intent(out) :: nh !! negative reversed cumulative hazard
        real(dp) :: power, lz
        call kappa3_power(x,params,lz,power)
        q = -power
        if (ieee_is_finite(power)) q = finite_sum(log(params(1)),-power)
        if (q >= 40.0_dp) then
            nh = lz-log(params(1))/params(1)-log1pexp(-q)/params(1)
            logh = log(-nh)
        else
            logh = burr_logsoftplus(q)-log(params(1))
            nh = negative_hazard(logh)
        end if
    end subroutine kappa3_hazard

    pure function kappa3_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! positive shape, location, positive scale
        real(dp) :: value, q, logh, nh, la
        value = negative_infinity(0.0_dp)
        if (x < params(2)) return
        la = log(params(1))
        if (x == params(2)) then
            value = -log(params(3))
            if (la < 0.0_dp) value = signed_exp(log(-la)-la,1.0_dp)-log(params(3))
            if (la > 0.0_dp) value = -la/params(1)-log(params(3))
            return
        end if
        call kappa3_hazard(x,params,q,logh,nh)
        if (q >= 40.0_dp) then
            value = -la/params(1)-log(params(3)) &
                -(1.0_dp+1.0_dp/params(1))*log1pexp(-q)
        else
            value = nh-log_distance(x,params(2))-log1pexp(-q)
        end if
    end function kappa3_lpdf

    pure function kappa3_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, location, positive scale
        real(dp) :: value, q, logh
        value = negative_infinity(0.0_dp)
        if (x <= params(2)) return
        call kappa3_hazard(x,params,q,logh,value)
    end function kappa3_lcdf

    pure function kappa3_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, location, positive scale
        real(dp) :: value, q, logh, nh
        value = 0.0_dp
        if (x <= params(2)) return
        call kappa3_hazard(x,params,q,logh,nh)
        value = hazard_complement(logh,nh)
    end function kappa3_lsf

    pure function kappa3_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper endpoint greater than a
        real(dp), intent(in) :: params(:) !! positive shape, location, positive scale
        real(dp) :: value, qa, qb, logh, nh, span, delta, diff, u, omitted, d
        value = kappa3_lcdf(b,params)
        if (a <= params(2)) return
        omitted = exp(kappa3_lcdf(a,params))+exp(kappa3_lsf(b,params))
        if (omitted <= 0.5_dp) then
            value = log1p_safe(-omitted)
            return
        end if
        call kappa3_hazard(a,params,qa,logh,nh)
        call kappa3_hazard(b,params,qb,logh,value)
        span = burr_logsoftplus(log_distance(b,a)-log_distance(a,params(2)))
        if (qb >= 40.0_dp) then
            diff = span
        else
            delta = log(params(1))+span
            if (delta < log(0.5_dp)) then
                u = log_expm1_from_log(delta)
                diff = burr_logsoftplus(-log1pexp(-qb)+u)-log(params(1))
            else if (qb >= 0.0_dp) then
                diff = delta
                if (delta < log(huge(1.0_dp))) then
                    d = exp(delta)
                    diff = log(d+log1pexp(-qa)-log1pexp(-qb))
                end if
                diff = diff-log(params(1))
            else if (qa <= 0.0_dp) then
                u = -log1pexp(-qa)+log_hazard_mass(delta)
                diff = u
                if (u >= log(epsilon(1.0_dp))-2.0_dp) &
                    diff = log(-log1p_safe(-exp(u)))
                diff = diff-log(params(1))
            else if (ieee_is_finite(qa)) then
                diff = log(log1pexp(qa)-log1pexp(qb))-log(params(1))
            else
                diff = log(-lognormal_logz(a,params))
            end if
        end if
        value = value+log_hazard_mass(diff)
    end function kappa3_interval


    pure function truncexpon_censored_loglikelihood(data, b, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: b !! positive finite standardized cutoff
        real(dp), intent(in), optional :: loc !! finite lower support edge, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(4)
        params = [b,0.0_dp,1.0_dp,0.0_dp]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (b <= 0.0_dp) return
        if (params(3) <= 0.0_dp) return
        if (b > 1.0_dp) then
            if (params(3) > huge(1.0_dp)/b) return
        end if
        params(4) = b*params(3)
        if (params(4) <= 0.0_dp) return
        if (params(2) > 0.0_dp) then
            if (params(4) > huge(1.0_dp)-params(2)) return
        end if
        if (params(2)+params(4) <= params(2)) return
        value = censored_loglikelihood(data,params,truncexpon_lpdf,truncexpon_lcdf, &
            truncexpon_lsf,truncexpon_interval)
    end function truncexpon_censored_loglikelihood

    pure function truncexpon_censored_nnlf(data, b, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: b !! positive finite standardized cutoff
        real(dp), intent(in), optional :: loc !! finite lower support edge, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -truncexpon_censored_loglikelihood(data, b, loc, scale)
    end function truncexpon_censored_nnlf

    pure function johnsonsb_censored_loglikelihood(data, a, b, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: a !! finite normal-logit shift
        real(dp), intent(in) :: b !! positive finite normal-logit scale
        real(dp), intent(in), optional :: loc !! finite lower support edge, default zero
        real(dp), intent(in), optional :: scale !! positive finite support width, default one
        real(dp) :: value, params(4)
        params = [a,b,0.0_dp,1.0_dp]
        if (present(loc)) params(3) = loc
        if (present(scale)) params(4) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (b <= 0.0_dp) return
        if (params(4) <= 0.0_dp) return
        if (.not. beta_parameters_valid([1.0_dp,1.0_dp,params(3),params(4)])) return
        if (params(3)+params(4) <= params(3)) return
        value = censored_loglikelihood(data,params,johnsonsb_lpdf,johnsonsb_lcdf, &
            johnsonsb_lsf,johnsonsb_interval)
    end function johnsonsb_censored_loglikelihood

    pure function johnsonsb_censored_nnlf(data, a, b, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: a !! finite normal-logit shift
        real(dp), intent(in) :: b !! positive finite normal-logit scale
        real(dp), intent(in), optional :: loc !! finite lower support edge, default zero
        real(dp), intent(in), optional :: scale !! positive finite support width, default one
        real(dp) :: value
        value = -johnsonsb_censored_loglikelihood(data, a, b, loc, scale)
    end function johnsonsb_censored_nnlf

    pure function kappa3_censored_loglikelihood(data, a, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: a !! positive finite kappa exponent
        real(dp), intent(in), optional :: loc !! finite lower support edge, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(3)
        params = [a,0.0_dp,1.0_dp]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (a <= 0.0_dp) return
        if (params(3) <= 0.0_dp) return
        value = censored_loglikelihood(data,params,kappa3_lpdf,kappa3_lcdf, &
            kappa3_lsf,kappa3_interval)
    end function kappa3_censored_loglikelihood

    pure function kappa3_censored_nnlf(data, a, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: a !! positive finite kappa exponent
        real(dp), intent(in), optional :: loc !! finite lower support edge, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -kappa3_censored_loglikelihood(data, a, loc, scale)
    end function kappa3_censored_nnlf


    pure function cosine_censored_loglikelihood(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite support center, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(3), sigma
        params = [1.0_dp,0.0_dp,1.0_dp]
        sigma = 1.0_dp
        if (present(loc)) params(2) = loc
        if (present(scale)) sigma = scale
        value = quiet_nan(0.0_dp)
        if (.not. valid_loc_scale(params(2),sigma)) return
        if (sigma > huge(1.0_dp)/3.1415926535897931_dp) return
        params(3) = 3.1415926535897931_dp*sigma
        if (params(3) > huge(1.0_dp)-abs(params(2))) return
        if (params(2)+params(3) <= params(2)) return
        if (params(2)-params(3) >= params(2)) return
        value = censored_loglikelihood(data,params,compact_lpdf,compact_lcdf, &
            compact_lsf,compact_interval)
    end function cosine_censored_loglikelihood

    pure function cosine_censored_nnlf(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite support center, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -cosine_censored_loglikelihood(data,loc,scale)
    end function cosine_censored_nnlf

    pure function anglit_censored_loglikelihood(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite support center, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(3), sigma
        params = [2.0_dp,0.0_dp,1.0_dp]
        sigma = 1.0_dp
        if (present(loc)) params(2) = loc
        if (present(scale)) sigma = scale
        value = quiet_nan(0.0_dp)
        if (.not. valid_loc_scale(params(2),sigma)) return
        params(3) = 0.78539816339744828_dp*sigma
        if (params(3) > huge(1.0_dp)-abs(params(2))) return
        if (params(2)+params(3) <= params(2)) return
        if (params(2)-params(3) >= params(2)) return
        value = censored_loglikelihood(data,params,compact_lpdf,compact_lcdf, &
            compact_lsf,compact_interval)
    end function anglit_censored_loglikelihood

    pure function anglit_censored_nnlf(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite support center, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -anglit_censored_loglikelihood(data,loc,scale)
    end function anglit_censored_nnlf

    pure function semicircular_censored_loglikelihood(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite support center, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(3), sigma
        params = [3.0_dp,0.0_dp,1.0_dp]
        sigma = 1.0_dp
        if (present(loc)) params(2) = loc
        if (present(scale)) sigma = scale
        value = quiet_nan(0.0_dp)
        if (.not. valid_loc_scale(params(2),sigma)) return
        params(3) = 1_dp*sigma
        if (params(3) > huge(1.0_dp)-abs(params(2))) return
        if (params(2)+params(3) <= params(2)) return
        if (params(2)-params(3) >= params(2)) return
        value = censored_loglikelihood(data,params,compact_lpdf,compact_lcdf, &
            compact_lsf,compact_interval)
    end function semicircular_censored_loglikelihood

    pure function semicircular_censored_nnlf(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite support center, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -semicircular_censored_loglikelihood(data,loc,scale)
    end function semicircular_censored_nnlf

    ! Independent integration of the SciPy public density specifications.
    ! The rounded half-width is the numerical support contract; see the guide.
    pure subroutine compact_gap(x, params, gap, logq)
        implicit none
        real(dp), intent(in) :: x !! finite observation or censoring limit
        real(dp), intent(in) :: params(:) !! family number, center, rounded half-width
        real(dp), intent(out) :: gap !! signed distance below the nearer support edge
        real(dp), intent(out) :: logq !! log distance to the nearer edge divided by half-width
        real(dp) :: unused
        if (x >= params(2)) then
            call beta_distances(x,[1.0_dp,1.0_dp,params(2),params(3)],unused,gap)
        else
            call beta_distances(-x,[1.0_dp,1.0_dp,-params(2),params(3)],unused,gap)
        end if
        logq = negative_infinity(0.0_dp)
        if (gap > 0.0_dp) logq = min(0.0_dp,log(gap)-log(params(3)))
    end subroutine compact_gap

    pure function compact_logsin(logx) result(value)
        implicit none
        real(dp), intent(in) :: logx !! log angle in (0, pi/2]
        real(dp) :: value
        value = logx
        ! The omitted relative correction is x**2/6 < epsilon/6.
        if (logx > 0.5_dp*log(epsilon(1.0_dp))) value = log(sin(exp(logx)))
    end function compact_logsin

    pure function compact_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! family number, center, rounded half-width
        real(dp) :: value, gap, lq
        call compact_gap(x,params,gap,lq)
        value = negative_infinity(0.0_dp)
        if (gap <= 0.0_dp) return
        value = compact_density(lq,params)
    end function compact_lpdf

    pure function compact_density(lq, params) result(value)
        implicit none
        real(dp), intent(in) :: lq !! log nearer endpoint gap divided by half-width
        real(dp), intent(in) :: params(:) !! family number, center, rounded half-width
        real(dp) :: value
        select case (int(params(1)))
        case (1)
            value = 2.0_dp*compact_logsin(lq+log(scifort_pi/2.0_dp))
        case (2)
            value = log(scifort_pi/4.0_dp)+compact_logsin(lq+log(scifort_pi/2.0_dp))
        case (3)
            value = log(2.0_dp/scifort_pi)+0.5_dp*(lq+log(2.0_dp-exp(lq)))
        end select
        value = value-log(params(3))
    end function compact_density

    pure function compact_small_tail(logq, family) result(value)
        implicit none
        real(dp), intent(in) :: logq !! log distance to the nearer endpoint in half-width units
        integer, intent(in) :: family !! 1 cosine, 2 anglit, 3 semicircular
        real(dp) :: value, q, term, total, v
        integer :: k
        q = exp(logq)
        select case (family)
        case (1)
            if (q <= 0.5_dp) then
                ! Integrate sin(pi*q/2)**2: q-sin(pi*q)/pi, without subtraction.
                v = (scifort_pi*q)**2
                term = 1.0_dp
                total = term
                do k = 2, 16
                    term = -term*v/real((2*k)*(2*k+1),dp)
                    total = total+term
                end do
                value = 3.0_dp*logq+2.0_dp*log(scifort_pi)-log(12.0_dp)+log(total)
            else
                value = log((q-sin(scifort_pi*q)/scifort_pi)/2.0_dp)
            end if
        case (2)
            value = 2.0_dp*compact_logsin(logq+log(scifort_pi/4.0_dp))
        case (3)
            if (q <= 0.25_dp) then
                ! Integrate sqrt(q*(2-q)) using the binomial series for sqrt(1-q/2).
                term = 1.0_dp
                total = term
                do k = 1, 24
                    term = term*(real(k,dp)-1.5_dp)*q/(2.0_dp*real(k,dp))
                    total = total+term*3.0_dp/real(2*k+3,dp)
                end do
                value = 1.5_dp*logq+log(4.0_dp*sqrt(2.0_dp)/(3.0_dp*scifort_pi)) &
                    +log(total)
            else
                value = log((acos(1.0_dp-q)-(1.0_dp-q)*sqrt(q*(2.0_dp-q)))/scifort_pi)
            end if
        end select
    end function compact_small_tail

    pure subroutine compact_tails(x, params, logp, logq)
        implicit none
        real(dp), intent(in) :: x !! finite censoring limit
        real(dp), intent(in) :: params(:) !! family number, center, rounded half-width
        real(dp), intent(out) :: logp !! log lower-tail probability
        real(dp), intent(out) :: logq !! log upper-tail probability
        real(dp) :: gap, lq, small, large
        call compact_gap(x,params,gap,lq)
        small = negative_infinity(0.0_dp)
        large = 0.0_dp
        if (gap > 0.0_dp) then
            small = min(-scifort_log_two,compact_small_tail(lq,int(params(1))))
            large = log1p_safe(-exp(small))
        end if
        logp = small
        logq = large
        if (x >= params(2)) then
            logp = large
            logq = small
        end if
    end subroutine compact_tails

    pure function compact_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! family number, center, rounded half-width
        real(dp) :: value, unused
        call compact_tails(x,params,value,unused)
    end function compact_lcdf

    pure function compact_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! family number, center, rounded half-width
        real(dp) :: value, unused
        call compact_tails(x,params,unused,value)
    end function compact_lsf

    pure function compact_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper endpoint greater than a
        real(dp), intent(in) :: params(:) !! family number, center, rounded half-width
        real(dp) :: value, ga, gb, la, lb, lw, anchor, correction, rho, lq, t(3), r(3)
        integer :: k
        call compact_gap(a,params,ga,la)
        call compact_gap(b,params,gb,lb)
        lw = log_distance(b,a)
        if (min(ga,gb) > 0.0_dp .and. lw < log(huge(1.0_dp))-1.0_dp) then
            if (lw-min(la,lb)-log(params(3)) <= log(0.001_dp)) then
                anchor = compact_lpdf(a,params)
                t = [0.5_dp,0.5_dp-sqrt(3.0_dp/5.0_dp)/2.0_dp, &
                    0.5_dp+sqrt(3.0_dp/5.0_dp)/2.0_dp]
                do k = 1, 3
                    ! Integrate in gap coordinates: physical quadrature nodes can
                    ! round to the original endpoints for adjacent observations.
                    if (a >= params(2)) then
                        rho = exp(lw-log(params(3))-la)
                        lq = la+log1p_safe(-t(k)*rho)
                    else if (b <= params(2)) then
                        rho = exp(lw-log(params(3))-la)
                        lq = min(0.0_dp,la+log1p_safe(t(k)*rho))
                    else
                        rho = exp(lw-log(params(3)))
                        lq = log1p_safe(-abs(expm1_safe(la)+t(k)*rho))
                    end if
                    r(k) = compact_density(lq,params)-anchor
                end do
                correction = (4.0_dp/9.0_dp)*expm1_safe(r(1)) &
                    +(5.0_dp/18.0_dp)*(expm1_safe(r(2))+expm1_safe(r(3)))
                value = anchor+lw+log1p_safe(correction)
                return
            end if
        end if
        value = tail_interval(a,b,params,compact_lcdf,compact_lsf)
        if (value > -scifort_log_two) then
            value = log1p_safe(-exp(bounded_logadd(compact_lcdf(a,params),compact_lsf(b,params))))
        end if
    end function compact_interval


    pure function finite_normal_lpdf(z) result(value)
        implicit none
        real(dp), intent(in) :: z !! normal coordinate, possibly infinite
        real(dp) :: value
        value = negative_infinity(0.0_dp)
        if (abs(z) > sqrt(huge(1.0_dp))*sqrt(2.0_dp)) return
        value = (-0.5_dp*z)*z-scifort_log_sqrt_two_pi
    end function finite_normal_lpdf

    pure function signed_exp(logabs, direction) result(value)
        implicit none
        real(dp), intent(in) :: logabs !! log nonnegative magnitude
        real(dp), intent(in) :: direction !! sign to apply, zero treated as positive
        real(dp) :: value
        value = sign(positive_infinity(0.0_dp),direction)
        if (logabs <= log(huge(1.0_dp))) value = sign(exp(logabs),direction)
    end function signed_exp

    pure subroutine fatiguelife_coordinates(x, params, w, v, logjac)
        implicit none
        real(dp), intent(in) :: x !! finite point strictly above location
        real(dp), intent(in) :: params(:) !! positive shape, location, scale
        real(dp), intent(out) :: w !! transformed normal coordinate
        real(dp), intent(out) :: v !! log standardized distance
        real(dp), intent(out) :: logjac !! log derivative of the normal coordinate
        real(dp) :: logw
        v = lognormal_logz(x,params)
        logjac = 0.5_dp*abs(v)+log1p_safe(exp(-abs(v)))-scifort_log_two &
            -log(params(1))-log_distance(x,params(2))
        w = 0.0_dp
        if (v == 0.0_dp) return
        if (abs(v) < 0.5_dp) then
            logw = log(abs(expm1_safe(v)))-0.5_dp*v-log(params(1))
        else
            logw = 0.5_dp*abs(v)+log1p_safe(-exp(-abs(v)))-log(params(1))
        end if
        w = signed_exp(logw,v)
    end subroutine fatiguelife_coordinates

    pure function fatiguelife_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! positive shape, location, scale
        real(dp) :: value, w, v, logjac
        value = negative_infinity(0.0_dp)
        if (x <= params(2)) return
        call fatiguelife_coordinates(x,params,w,v,logjac)
        value = finite_normal_lpdf(w)+logjac
    end function fatiguelife_lpdf

    pure function fatiguelife_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, location, scale
        real(dp) :: value, w, v, logjac
        value = negative_infinity(0.0_dp)
        if (x <= params(2)) return
        call fatiguelife_coordinates(x,params,w,v,logjac)
        value = normal_logcdf(w)
    end function fatiguelife_lcdf

    pure function fatiguelife_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, location, scale
        real(dp) :: value, w, v, logjac
        value = 0.0_dp
        if (x <= params(2)) return
        call fatiguelife_coordinates(x,params,w,v,logjac)
        value = normal_logsf(w)
    end function fatiguelife_lsf

    pure function fatiguelife_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper interval endpoint, > a
        real(dp), intent(in) :: params(:) !! positive shape, location, scale
        real(dp) :: value, wa, wb, va, vb, logjac, logspan, logwidth, omitted
        value = fatiguelife_lcdf(b,params)
        if (a <= params(2)) return
        call fatiguelife_coordinates(a,params,wa,va,logjac)
        call fatiguelife_coordinates(b,params,wb,vb,logjac)
        logspan = burr_logsoftplus(log_distance(b,a)-log_distance(a,params(2)))
        logwidth = log_expm1_from_log(logspan-scifort_log_two) &
            +bounded_logadd(0.5_dp*va,-0.5_dp*vb)-log(params(1))
        value = normal_small_interval(wa,logwidth)
        if (ieee_is_nan(value)) then
            value = tail_interval(a,b,params,fatiguelife_lcdf,fatiguelife_lsf)
            if (value > -scifort_log_two) then
                omitted = bounded_logadd(fatiguelife_lcdf(a,params),fatiguelife_lsf(b,params))
                value = log1p_safe(-exp(omitted))
            end if
        end if
    end function fatiguelife_interval

    pure function finite_sum(a, b) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite first summand
        real(dp), intent(in) :: b !! finite second summand
        real(dp) :: value
        if (a > 0.0_dp .and. b > 0.0_dp) then
            if (b > huge(1.0_dp)-a) then
                value = positive_infinity(0.0_dp)
                return
            end if
        else if (a < 0.0_dp .and. b < 0.0_dp) then
            if (b < -huge(1.0_dp)-a) then
                value = negative_infinity(0.0_dp)
                return
            end if
        end if
        value = a+b
    end function finite_sum

    pure subroutine johnsonsu_coordinates(x, params, q, t, loghypot)
        implicit none
        real(dp), intent(in) :: x !! finite evaluation point
        real(dp), intent(in) :: params(:) !! finite a, positive b, location, scale
        real(dp), intent(out) :: q !! transformed normal coordinate
        real(dp), intent(out) :: t !! inverse hyperbolic sine of standardized point
        real(dp), intent(out) :: loghypot !! log sqrt(1+z*z)
        real(dp) :: z, logz, product, halfsum
        logz = cauchy_log_distance(x,params(3:4))
        z = location_scale_z(x,params(3:4))
        loghypot = 0.5_dp*log1pexp(2.0_dp*logz)
        if (logz < log(epsilon(1.0_dp))-2.0_dp) then
            t = z
            product = signed_exp(log(params(2))+logz,sign(1.0_dp,z))
            if (x < params(3)) product = -abs(product)
            q = finite_sum(params(1),product)
            return
        else if (logz > log(1.0e8_dp)) then
            t = sign(logz+scifort_log_two,z)
        else
            t = asinh(z)
        end if
        if (abs(t) > 1.0_dp) then
            if (params(2) > huge(1.0_dp)/abs(t)) then
                ! Halving before the product allows an opposing finite a
                ! to rescue a product between huge() and 2*huge().
                q = sign(positive_infinity(0.0_dp),t)
                if (0.5_dp*params(2) > huge(1.0_dp)/abs(t)) return
                halfsum = finite_sum(0.5_dp*params(1),(0.5_dp*params(2))*t)
                if (abs(halfsum) > 0.5_dp*huge(1.0_dp)) return
                q = 2.0_dp*halfsum
                return
            end if
        end if
        q = finite_sum(params(1),params(2)*t)
    end subroutine johnsonsu_coordinates

    pure function johnsonsu_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! finite a, positive b, location, scale
        real(dp) :: value, q, t, loghypot
        call johnsonsu_coordinates(x,params,q,t,loghypot)
        value = finite_normal_lpdf(q)+log(params(2))-log(params(4))-loghypot
    end function johnsonsu_lpdf

    pure function johnsonsu_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! finite a, positive b, location, scale
        real(dp) :: value, q, t, loghypot
        call johnsonsu_coordinates(x,params,q,t,loghypot)
        value = normal_logcdf(q)
    end function johnsonsu_lcdf

    pure function johnsonsu_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! finite a, positive b, location, scale
        real(dp) :: value, q, t, loghypot
        call johnsonsu_coordinates(x,params,q,t,loghypot)
        value = normal_logsf(q)
    end function johnsonsu_lsf

    pure function johnsonsu_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper interval endpoint, > a
        real(dp), intent(in) :: params(:) !! finite a, positive b, location, scale
        real(dp) :: value, qa, qb, ta, tb, la, lb, logr, logwidth, r, omitted
        call johnsonsu_coordinates(a,params,qa,ta,la)
        call johnsonsu_coordinates(b,params,qb,tb,lb)
        ! tanh((asinh(zb)-asinh(za))/2)=(zb-za)/(hypot(1,za)+hypot(1,zb)).
        logr = log_distance(b,a)-log(params(4))-bounded_logadd(la,lb)
        if (logr < log(epsilon(1.0_dp))-2.0_dp) then
            logwidth = scifort_log_two+logr
        else if (logr <= -scifort_log_two) then
            r = exp(logr)
            logwidth = log(log1p_safe(2.0_dp*r/(1.0_dp-r)))
        else
            logwidth = log(tb-ta)
        end if
        logwidth = log(params(2))+logwidth
        value = normal_small_interval(qa,logwidth)
        if (ieee_is_nan(value)) then
            value = tail_interval(a,b,params,johnsonsu_lcdf,johnsonsu_lsf)
            if (value > -scifort_log_two) then
                omitted = bounded_logadd(johnsonsu_lcdf(a,params),johnsonsu_lsf(b,params))
                value = log1p_safe(-exp(omitted))
            end if
        end if
    end function johnsonsu_interval

    pure function genexpon_log_remainder(logv) result(value)
        implicit none
        real(dp), intent(in) :: logv !! log nonnegative v in v+expm1(-v)
        real(dp) :: value, v, term, series
        integer :: k
        if (logv < log(epsilon(1.0_dp))-2.0_dp) then
            value = 2.0_dp*logv-scifort_log_two
        else if (logv <= log(0.5_dp)) then
            v = exp(logv)
            term = 1.0_dp
            series = 1.0_dp
            do k = 3, 22
                term = -term*v/real(k,dp)
                series = series+term
            end do
            value = 2.0_dp*logv-scifort_log_two+log(series)
        else if (logv > log(40.0_dp)) then
            value = logv+log1p_safe(-exp(-logv))
        else
            v = exp(logv)
            value = log(v+expm1_safe(-v))
        end if
    end function genexpon_log_remainder

    pure function genexpon_logh(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite evaluation point
        real(dp), intent(in) :: params(:) !! positive a, b, c, location, scale
        real(dp) :: value, logz
        value = negative_infinity(0.0_dp)
        if (x <= params(4)) return
        logz = lognormal_logz(x,[1.0_dp,params(4),params(5)])
        value = bounded_logadd(log(params(1))+logz,log(params(2))-log(params(3)) &
            +genexpon_log_remainder(log(params(3))+logz))
    end function genexpon_logh

    pure function genexpon_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! positive a, b, c, location, scale
        real(dp) :: value, logv
        value = negative_infinity(0.0_dp)
        if (x < params(4)) return
        value = log(params(1))-log(params(5))
        if (x == params(4)) return
        logv = log(params(3))+lognormal_logz(x,[1.0_dp,params(4),params(5)])
        value = bounded_logadd(log(params(1)),log(params(2))+log_hazard_mass(logv)) &
            -log(params(5))+negative_hazard(genexpon_logh(x,params))
    end function genexpon_lpdf

    pure function genexpon_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! positive a, b, c, location, scale
        real(dp) :: value
        value = log_hazard_mass(genexpon_logh(x,params))
    end function genexpon_lcdf

    pure function genexpon_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! positive a, b, c, location, scale
        real(dp) :: value
        value = negative_hazard(genexpon_logh(x,params))
    end function genexpon_lsf

    pure function genexpon_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper interval endpoint, > a
        real(dp), intent(in) :: params(:) !! positive a, b, c, location, scale
        real(dp) :: value, lower, logdelta, logv, logdv, first, second, third, omitted
        value = negative_infinity(0.0_dp)
        if (b <= params(4)) return
        lower = max(a,params(4))
        omitted = exp(genexpon_lcdf(lower,params))+exp(genexpon_lsf(b,params))
        if (omitted <= 0.5_dp) then
            value = log1p_safe(-omitted)
            return
        end if
        logdelta = log_distance(b,lower)-log(params(5))
        logv = negative_infinity(0.0_dp)
        if (lower > params(4)) then
            logv = log(params(3))+lognormal_logz(lower,[1.0_dp,params(4),params(5)])
        end if
        logdv = log(params(3))+logdelta
        first = log(params(1))+logdelta
        second = log(params(2))+logdelta+log_hazard_mass(logv)
        third = log(params(2))-log(params(3))+negative_hazard(logv) &
            +genexpon_log_remainder(logdv)
        value = negative_hazard(genexpon_logh(lower,params)) &
            +log_hazard_mass(bounded_logadd(first,bounded_logadd(second,third)))
    end function genexpon_interval


    pure function fatiguelife_censored_loglikelihood(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! positive finite fatigue-life shape
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(3)
        params = [c,0.0_dp,1.0_dp]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (c <= 0.0_dp) return
        if (params(3) <= 0.0_dp) return
        value = censored_loglikelihood(data,params,fatiguelife_lpdf,fatiguelife_lcdf, &
            fatiguelife_lsf,fatiguelife_interval)
    end function fatiguelife_censored_loglikelihood

    pure function fatiguelife_censored_nnlf(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! positive finite fatigue-life shape
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -fatiguelife_censored_loglikelihood(data, c, loc, scale)
    end function fatiguelife_censored_nnlf

    pure function johnsonsu_censored_loglikelihood(data, a, b, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: a !! finite first Johnson shape
        real(dp), intent(in) :: b !! positive finite Johnson scale shape
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(4)
        params = [a,b,0.0_dp,1.0_dp]
        if (present(loc)) params(3) = loc
        if (present(scale)) params(4) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (b <= 0.0_dp) return
        if (params(4) <= 0.0_dp) return
        value = censored_loglikelihood(data,params,johnsonsu_lpdf,johnsonsu_lcdf, &
            johnsonsu_lsf,johnsonsu_interval)
    end function johnsonsu_censored_loglikelihood

    pure function johnsonsu_censored_nnlf(data, a, b, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: a !! finite first Johnson shape
        real(dp), intent(in) :: b !! positive finite Johnson scale shape
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -johnsonsu_censored_loglikelihood(data, a, b, loc, scale)
    end function johnsonsu_censored_nnlf

    pure function genexpon_censored_loglikelihood(data, a, b, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: a !! positive finite baseline rate
        real(dp), intent(in) :: b !! positive finite additional rate
        real(dp), intent(in) :: c !! positive finite rate-transition parameter
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(5)
        params = [a,b,c,0.0_dp,1.0_dp]
        if (present(loc)) params(4) = loc
        if (present(scale)) params(5) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (min(a,b,c) <= 0.0_dp) return
        if (params(5) <= 0.0_dp) return
        value = censored_loglikelihood(data,params,genexpon_lpdf,genexpon_lcdf, &
            genexpon_lsf,genexpon_interval)
    end function genexpon_censored_loglikelihood

    pure function genexpon_censored_nnlf(data, a, b, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: a !! positive finite baseline rate
        real(dp), intent(in) :: b !! positive finite additional rate
        real(dp), intent(in) :: c !! positive finite rate-transition parameter
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -genexpon_censored_loglikelihood(data, a, b, c, loc, scale)
    end function genexpon_censored_nnlf


    pure function weibull_log_increment(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower endpoint strictly above location
        real(dp), intent(in) :: b !! finite upper endpoint, > a
        real(dp), intent(in) :: params(:) !! positive Weibull shape, location, scale
        real(dp) :: value, logspan, logdelta
        ! exp(logt(b))-exp(logt(a)), preserving the original endpoint gap.
        logspan = log(params(1))+burr_logsoftplus( &
            log_distance(b,a)-log_distance(a,params(2)))
        if (logspan > log(50.0_dp)) then
            value = weibull_log_hazard(b,params)
        else
            logdelta = log_expm1_from_log(logspan)
            value = weibull_log_hazard(a,params)+logdelta
        end if
    end function weibull_log_increment

    pure function exponpow_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! positive power shape, location, scale
        real(dp) :: value, logt, t
        value = negative_infinity(0.0_dp)
        if (x < params(2)) return
        if (x == params(2)) then
            if (params(1) < 1.0_dp) value = positive_infinity(0.0_dp)
            if (params(1) == 1.0_dp) value = -log(params(3))
            return
        end if
        logt = weibull_log_hazard(x,params)
        if (logt > log(log(huge(1.0_dp)))) return
        t = exp(logt)
        value = log(params(1))-log_distance(x,params(2))+logt+t-expm1_safe(t)
    end function exponpow_lpdf

    pure function exponpow_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! positive power shape, location, scale
        real(dp) :: value
        value = negative_infinity(0.0_dp)
        if (x <= params(2)) return
        value = log_hazard_mass(log_expm1_from_log(weibull_log_hazard(x,params)))
    end function exponpow_lcdf

    pure function exponpow_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! positive power shape, location, scale
        real(dp) :: value, logt
        value = 0.0_dp
        if (x <= params(2)) return
        logt = weibull_log_hazard(x,params)
        value = negative_hazard(log_expm1_from_log(logt))
        if (logt >= log(tiny(1.0_dp)) .and. logt <= log(log(huge(1.0_dp)))) &
            value = -expm1_safe(exp(logt))
    end function exponpow_lsf

    pure function exponpow_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper interval endpoint, > a
        real(dp), intent(in) :: params(:) !! positive power shape, location, scale
        real(dp) :: value, logt, logdelta, omitted
        value = exponpow_lcdf(b,params)
        if (a <= params(2)) return
        omitted = exp(exponpow_lcdf(a,params))+exp(exponpow_lsf(b,params))
        if (omitted <= 0.5_dp) then
            value = log1p_safe(-omitted)
            return
        end if
        value = exponpow_lsf(a,params)
        if (.not. ieee_is_finite(value)) return
        logt = weibull_log_hazard(a,params)
        logdelta = weibull_log_increment(a,b,params)
        value = value+log_hazard_mass(exp(logt)+log_expm1_from_log(logdelta))
    end function exponpow_interval

    pure subroutine exponweib_hazard(x, params, logt, logh, nh)
        implicit none
        real(dp), intent(in) :: x !! finite point strictly above location
        real(dp), intent(in) :: params(:) !! positive a, c, location, scale
        real(dp), intent(out) :: logt !! c*log(z), possibly infinite
        real(dp), intent(out) :: logh !! log(-a*log(1-exp(-z**c)))
        real(dp), intent(out) :: nh !! negative reversed cumulative hazard
        real(dp) :: v, base, a
        a = params(1)
        logt = weibull_log_hazard(x,params(2:4))
        if (logt < log(epsilon(1.0_dp))-2.0_dp) then
            ! For tiny z**c, log(F_weibull) = c*log(z) to relative roundoff.
            ! Form a*c*abs(log(z)) without requiring c*log(z) to be finite.
            v = lognormal_logz(x,params(2:4))
            logh = log(a)+log(params(2))+log(-v)
            nh = negative_hazard(logh)
            if (ieee_is_finite(logt)) then
                if (a <= 1.0_dp) then
                    nh = a*logt
                else if (abs(logt) <= huge(1.0_dp)/a) then
                    nh = a*logt
                end if
            end if
        else if (logt > log(36.0_dp)) then
            ! -log(1-exp(-t)) / exp(-t) differs from one by < 2e-16.
            logh = log(a)+negative_hazard(logt)
            nh = negative_hazard(logh)
        else
            base = -log_hazard_mass(logt)
            logh = log(a)+log(base)
            nh = negative_hazard(logh)
            if (base <= 1.0_dp) then
                nh = -a*base
            else if (a <= huge(1.0_dp)/base) then
                nh = -a*base
            end if
        end if
    end subroutine exponweib_hazard

    pure function exponweib_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! positive a, c, location, scale
        real(dp) :: value, logt, logh, nh, product
        value = negative_infinity(0.0_dp)
        if (x < params(3)) return
        if (x == params(3)) then
            if (params(1) > 1.0_dp) then
                if (params(2) > huge(1.0_dp)/params(1)) return
            end if
            product = params(1)*params(2)
            if (product < 1.0_dp) value = positive_infinity(0.0_dp)
            if (product == 1.0_dp) value = -log(params(4))
            return
        end if
        call exponweib_hazard(x,params,logt,logh,nh)
        if (.not. ieee_is_finite(nh)) return
        if (logt > log(huge(1.0_dp))) return
        value = log(params(1))+log(params(2))-log_distance(x,params(3))+nh
        if (logt >= log(epsilon(1.0_dp))-2.0_dp) &
            value = value+logt-log_expm1_from_log(logt)
    end function exponweib_lpdf

    pure function exponweib_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! positive a, c, location, scale
        real(dp) :: value, logt, logh
        value = negative_infinity(0.0_dp)
        if (x <= params(3)) return
        call exponweib_hazard(x,params,logt,logh,value)
    end function exponweib_lcdf

    pure function exponweib_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! positive a, c, location, scale
        real(dp) :: value, logt, logh, nh
        value = 0.0_dp
        if (x <= params(3)) return
        call exponweib_hazard(x,params,logt,logh,nh)
        value = hazard_complement(logh,nh)
    end function exponweib_lsf

    pure function exponweib_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper interval endpoint, > a
        real(dp), intent(in) :: params(:) !! positive a, c, location, scale
        real(dp) :: value, logta, logtb, logha, loghb, nha, nhb, logdelta, logu, omitted
        value = exponweib_lcdf(b,params)
        if (a <= params(3)) return
        omitted = exp(exponweib_lcdf(a,params))+exp(exponweib_lsf(b,params))
        if (omitted <= 0.5_dp) then
            value = log1p_safe(-omitted)
            return
        end if
        call exponweib_hazard(a,params,logta,logha,nha)
        call exponweib_hazard(b,params,logtb,loghb,nhb)
        if (logtb < log(epsilon(1.0_dp))-2.0_dp) then
            logdelta = log(params(1))+log(params(2))+burr_logsoftplus( &
                log_distance(b,a)-log_distance(a,params(3)))
        else if (logta == negative_infinity(0.0_dp)) then
            logdelta = logha
        else
            value = negative_infinity(0.0_dp)
            if (logta > log(huge(1.0_dp))) return
            logdelta = weibull_log_increment(a,b,params(2:4))
            logu = -exp(logta)+log_hazard_mass(logdelta)-log_hazard_mass(logta)
            logdelta = log(params(1))+burr_logsoftplus(logu)
        end if
        value = nhb+log_hazard_mass(logdelta)
    end function exponweib_interval

    pure subroutine mielke_hazard(x, params, q, logh, nh)
        implicit none
        real(dp), intent(in) :: x !! finite point strictly above location
        real(dp), intent(in) :: params(:) !! positive k, s, location, scale
        real(dp), intent(out) :: q !! -s*log(z), possibly infinite
        real(dp), intent(out) :: logh !! log((k/s)*log1p(z**(-s)))
        real(dp), intent(out) :: nh !! negative reversed cumulative hazard
        real(dp) :: v, base, ratio
        q = -weibull_log_hazard(x,params(2:4))
        if (q >= 40.0_dp) then
            v = -lognormal_logz(x,params(2:4))
            logh = log(params(1))+log(v)
            nh = negative_hazard(logh)
            if (params(1) <= 1.0_dp) then
                nh = -params(1)*v
            else if (v <= huge(1.0_dp)/params(1)) then
                nh = -params(1)*v
            end if
        else
            logh = log(params(1))-log(params(2))+burr_logsoftplus(q)
            nh = negative_hazard(logh)
            ! Direct finite ratios/products retain ordinary-tail accuracy.
            if (params(2) < 1.0_dp) then
                if (params(1) > huge(1.0_dp)*params(2)) return
            end if
            ratio = params(1)/params(2)
            base = log1pexp(q)
            if (ratio < tiny(1.0_dp) .or. base < tiny(1.0_dp)) return
            if (ratio > 1.0_dp) then
                if (base > huge(1.0_dp)/ratio) return
            end if
            nh = -ratio*base
        end if
    end subroutine mielke_hazard

    pure function mielke_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! positive k, s, location, scale
        real(dp) :: value, q, logh, nh
        value = negative_infinity(0.0_dp)
        if (x < params(3)) return
        if (x == params(3)) then
            if (params(1) < 1.0_dp) value = positive_infinity(0.0_dp)
            if (params(1) == 1.0_dp) value = -log(params(4))
            return
        end if
        call mielke_hazard(x,params,q,logh,nh)
        value = log(params(1))-log_distance(x,params(3))-log1pexp(-q)+nh
    end function mielke_lpdf

    pure function mielke_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! positive k, s, location, scale
        real(dp) :: value, q, logh
        value = negative_infinity(0.0_dp)
        if (x <= params(3)) return
        call mielke_hazard(x,params,q,logh,value)
    end function mielke_lcdf

    pure function mielke_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! positive k, s, location, scale
        real(dp) :: value, q, logh, nh
        value = 0.0_dp
        if (x <= params(3)) return
        call mielke_hazard(x,params,q,logh,nh)
        value = hazard_complement(logh,nh)
    end function mielke_lsf

    pure function mielke_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper interval endpoint, > a
        real(dp), intent(in) :: params(:) !! positive k, s, location, scale
        real(dp) :: value, qa, qb, logh, nh, logspan, logdelta, logdiff, logu, delta, omitted
        value = mielke_lcdf(b,params)
        if (a <= params(3)) return
        omitted = exp(mielke_lcdf(a,params))+exp(mielke_lsf(b,params))
        if (omitted <= 0.5_dp) then
            value = log1p_safe(-omitted)
            return
        end if
        call mielke_hazard(a,params,qa,logh,nh)
        call mielke_hazard(b,params,qb,logh,value)
        logspan = burr_logsoftplus(log_distance(b,a)-log_distance(a,params(3)))
        if (qb >= 40.0_dp) then
            logdiff = log(params(1))+logspan
        else
            logdelta = log(params(2))+logspan
            if (logdelta < log(0.5_dp)) then
                logu = log_expm1_from_log(logdelta)
                logdiff = burr_logsoftplus(-log1pexp(-qb)+logu)
            else if (qb >= 0.0_dp) then
                logdiff = logdelta
                if (logdelta < log(huge(1.0_dp))) then
                    delta = exp(logdelta)
                    logdiff = log(delta+log1pexp(-qa)-log1pexp(-qb))
                end if
            else if (qa <= 0.0_dp) then
                logu = -log1pexp(-qa)+log_hazard_mass(logdelta)
                logdiff = logu
                if (logu >= log(epsilon(1.0_dp))-2.0_dp) &
                    logdiff = log(-log1p_safe(-exp(logu)))
            else if (ieee_is_finite(qa)) then
                logdiff = log(log1pexp(qa)-log1pexp(qb))
            else
                logdiff = log(params(2))+log(-lognormal_logz(a,params(2:4)))
            end if
            logdiff = log(params(1))-log(params(2))+logdiff
        end if
        value = value+log_hazard_mass(logdiff)
    end function mielke_interval


    pure function exponpow_censored_loglikelihood(data, b, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: b !! positive finite power shape
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(3)
        params = [b,0.0_dp,1.0_dp]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (b <= 0.0_dp) return
        if (params(3) <= 0.0_dp) return
        value = censored_loglikelihood(data,params,exponpow_lpdf,exponpow_lcdf, &
            exponpow_lsf,exponpow_interval)
    end function exponpow_censored_loglikelihood

    pure function exponpow_censored_nnlf(data, b, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: b !! positive finite power shape
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -exponpow_censored_loglikelihood(data, b, loc, scale)
    end function exponpow_censored_nnlf

    pure function exponweib_censored_loglikelihood(data, a, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: a !! positive finite exponentiation shape
        real(dp), intent(in) :: c !! positive finite Weibull shape
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(4)
        params = [a,c,0.0_dp,1.0_dp]
        if (present(loc)) params(3) = loc
        if (present(scale)) params(4) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (min(a,c) <= 0.0_dp) return
        if (params(4) <= 0.0_dp) return
        value = censored_loglikelihood(data,params,exponweib_lpdf,exponweib_lcdf, &
            exponweib_lsf,exponweib_interval)
    end function exponweib_censored_loglikelihood

    pure function exponweib_censored_nnlf(data, a, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: a !! positive finite exponentiation shape
        real(dp), intent(in) :: c !! positive finite Weibull shape
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -exponweib_censored_loglikelihood(data, a, c, loc, scale)
    end function exponweib_censored_nnlf

    pure function mielke_censored_loglikelihood(data, k, s, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: k !! positive finite first shape
        real(dp), intent(in) :: s !! positive finite second shape
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(4)
        params = [k,s,0.0_dp,1.0_dp]
        if (present(loc)) params(3) = loc
        if (present(scale)) params(4) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (min(k,s) <= 0.0_dp) return
        if (params(4) <= 0.0_dp) return
        value = censored_loglikelihood(data,params,mielke_lpdf,mielke_lcdf, &
            mielke_lsf,mielke_interval)
    end function mielke_censored_loglikelihood

    pure function mielke_censored_nnlf(data, k, s, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: k !! positive finite first shape
        real(dp), intent(in) :: s !! positive finite second shape
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -mielke_censored_loglikelihood(data, k, s, loc, scale)
    end function mielke_censored_nnlf


    pure function normal_log_mills(logz) result(value)
        implicit none
        real(dp), intent(in) :: logz !! log positive standardized point, z >= 20
        real(dp) :: value, inverse_square, term, series
        integer :: k
        ! Q(z)*z/phi(z), DLMF 7.12.1. At z >= 20 the first omitted
        ! term after twelve corrections is below 2e-20.
        inverse_square = exp(-2.0_dp*logz)
        term = 1.0_dp
        series = 1.0_dp
        do k = 1, 12
            term = -term*real(2*k-1,dp)*inverse_square
            series = series+term
        end do
        value = log(series)
    end function normal_log_mills

    pure function normal_log_hazard(z, logz) result(value)
        implicit none
        real(dp), intent(in) :: z !! standardized point, possibly infinite
        real(dp), intent(in) :: logz !! log absolute standardized point
        real(dp) :: value
        if (z >= 20.0_dp) then
            value = logz-normal_log_mills(logz)
        else if (z < -20.0_dp) then
            value = negative_hazard(2.0_dp*logz-scifort_log_two)-scifort_log_sqrt_two_pi
        else
            value = normal_logpdf(z)-normal_logsf(z)
        end if
    end function normal_log_hazard

    pure function power_normal_logh(z, logz, c) result(value)
        implicit none
        real(dp), intent(in) :: z !! standardized point, possibly infinite
        real(dp), intent(in) :: logz !! log absolute standardized point
        real(dp), intent(in) :: c !! positive finite power
        real(dp) :: value, correction
        if (z >= 20.0_dp) then
            correction = 2.0_dp*(logz+scifort_log_sqrt_two_pi-normal_log_mills(logz))
            correction = correction*exp(-2.0_dp*logz)
            value = log(c)+2.0_dp*logz-scifort_log_two+log1p_safe(correction)
        else if (z < -8.0_dp) then
            ! -log(1-Phi(z)) / Phi(z) differs from one by < 4e-16.
            value = log(c)+normal_logcdf(z)
        else
            value = log(c)+log(-normal_logsf(z))
        end if
    end function power_normal_logh

    pure subroutine power_coordinates(x, params, z, logz)
        implicit none
        real(dp), intent(in) :: x !! finite evaluation point
        real(dp), intent(in) :: params(:) !! power, location, positive scale
        real(dp), intent(out) :: z !! standardized point, possibly infinite
        real(dp), intent(out) :: logz !! log absolute standardized point
        z = location_scale_z(x,params(2:3))
        logz = cauchy_log_distance(x,params(2:3))
        if (z /= 0.0_dp .and. ieee_is_finite(z)) logz = log(abs(z))
    end subroutine power_coordinates

    pure function powernorm_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! power, location, positive scale
        real(dp) :: value, z, logz
        call power_coordinates(x,params,z,logz)
        value = log(params(1))+normal_log_hazard(z,logz)-log(params(3)) &
            +negative_hazard(power_normal_logh(z,logz,params(1)))
    end function powernorm_lpdf

    pure function powernorm_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! power, location, positive scale
        real(dp) :: value, z, logz
        call power_coordinates(x,params,z,logz)
        value = log_hazard_mass(power_normal_logh(z,logz,params(1)))
    end function powernorm_lcdf

    pure function powernorm_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! power, location, positive scale
        real(dp) :: value, z, logz
        call power_coordinates(x,params,z,logz)
        value = negative_hazard(power_normal_logh(z,logz,params(1)))
    end function powernorm_lsf

    pure function power_normal_small_interval(z, logz, logwidth, c) result(value)
        implicit none
        real(dp), intent(in) :: z !! standardized lower endpoint, possibly infinite
        real(dp), intent(in) :: logz !! log absolute standardized lower endpoint
        real(dp), intent(in) :: logwidth !! log width in normal coordinates
        real(dp), intent(in) :: c !! positive finite power
        real(dp) :: value, logdelta, q, nodes(3), ratios(3), w, hazard, correction
        integer :: i
        value = quiet_nan(0.0_dp)
        if (z >= 20.0_dp) then
            if (logwidth-logz > log(0.001_dp)) return
        else
            if (logwidth+max(0.0_dp,logz) > log(0.001_dp)) return
        end if
        hazard = normal_log_hazard(z,logz)
        if (.not. ieee_is_finite(hazard)) then
            value = negative_infinity(0.0_dp)
            return
        end if
        if (z >= 1.0e8_dp) then
            ! The Mills correction to the integrated hazard is <= 1/z**2.
            q = exp(logwidth-logz)
            correction = log1p_safe(0.5_dp*q)
        else
            ! Integrate the normal hazard, rather than subtract nearly equal
            ! log survivals. Three-point Gauss-Legendre, DLMF 3.5(v).
            nodes = [0.5_dp,0.5_dp-sqrt(3.0_dp/5.0_dp)/2.0_dp, &
                0.5_dp+sqrt(3.0_dp/5.0_dp)/2.0_dp]
            do i = 1, 3
                w = z+nodes(i)*exp(logwidth)
                if (w == 0.0_dp) then
                    ratios(i) = normal_log_hazard(w,negative_infinity(0.0_dp))-hazard
                else
                    ratios(i) = normal_log_hazard(w,log(abs(w)))-hazard
                end if
            end do
            correction = log1p_safe((4.0_dp/9.0_dp)*expm1_safe(ratios(1)) &
                +(5.0_dp/18.0_dp)*(expm1_safe(ratios(2))+expm1_safe(ratios(3))))
        end if
        logdelta = log(c)+logwidth+hazard+correction
        value = negative_hazard(power_normal_logh(z,logz,c))+log_hazard_mass(logdelta)
    end function power_normal_small_interval

    pure function powernorm_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper interval endpoint, > a
        real(dp), intent(in) :: params(:) !! power, location, positive scale
        real(dp) :: value, z, logz, omitted
        call power_coordinates(a,params,z,logz)
        value = power_normal_small_interval(z,logz, &
            log_distance(b,a)-log(params(3)),params(1))
        if (ieee_is_nan(value)) then
            value = tail_interval(a,b,params,powernorm_lcdf,powernorm_lsf)
            if (value > -scifort_log_two) then
                omitted = bounded_logadd(powernorm_lcdf(a,params),powernorm_lsf(b,params))
                value = log1p_safe(-exp(omitted))
            end if
        end if
    end function powernorm_interval

    pure subroutine powerlog_coordinates(x, params, z, logz)
        implicit none
        real(dp), intent(in) :: x !! finite point strictly above location
        real(dp), intent(in) :: params(:) !! power, lognormal shape, location, scale
        real(dp), intent(out) :: z !! normal coordinate, possibly infinite
        real(dp), intent(out) :: logz !! log absolute normal coordinate
        real(dp) :: lz
        lz = lognormal_logz(x,params(2:4))
        z = lognormal_standardized(x,params(2:4))
        logz = negative_infinity(0.0_dp)
        if (lz /= 0.0_dp) logz = log(abs(lz))-log(params(2))
        if (z /= 0.0_dp .and. ieee_is_finite(z)) logz = log(abs(z))
    end subroutine powerlog_coordinates

    pure function powerlognorm_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! power, lognormal shape, location, scale
        real(dp) :: value, z, logz
        value = negative_infinity(0.0_dp)
        if (x <= params(3)) return
        call powerlog_coordinates(x,params,z,logz)
        value = log(params(1))+normal_log_hazard(z,logz)-log(params(2)) &
            -log_distance(x,params(3))+negative_hazard(power_normal_logh(z,logz,params(1)))
    end function powerlognorm_lpdf

    pure function powerlognorm_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! power, lognormal shape, location, scale
        real(dp) :: value, z, logz
        value = negative_infinity(0.0_dp)
        if (x <= params(3)) return
        call powerlog_coordinates(x,params,z,logz)
        value = log_hazard_mass(power_normal_logh(z,logz,params(1)))
    end function powerlognorm_lcdf

    pure function powerlognorm_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! power, lognormal shape, location, scale
        real(dp) :: value, z, logz
        value = 0.0_dp
        if (x <= params(3)) return
        call powerlog_coordinates(x,params,z,logz)
        value = negative_hazard(power_normal_logh(z,logz,params(1)))
    end function powerlognorm_lsf

    pure function powerlognorm_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper interval endpoint, > a
        real(dp), intent(in) :: params(:) !! power, lognormal shape, location, scale
        real(dp) :: value, z, logz, logr, logspan, omitted
        value = negative_infinity(0.0_dp)
        if (b <= params(3)) return
        if (a <= params(3)) then
            value = powerlognorm_lcdf(b,params)
            return
        end if
        logr = log_distance(b,a)-log_distance(a,params(3))
        logspan = logr
        if (logr >= log(epsilon(1.0_dp))-2.0_dp) logspan = log(log1pexp(logr))
        call powerlog_coordinates(a,params,z,logz)
        value = power_normal_small_interval(z,logz,logspan-log(params(2)),params(1))
        if (ieee_is_nan(value)) then
            value = tail_interval(a,b,params,powerlognorm_lcdf,powerlognorm_lsf)
            if (value > -scifort_log_two) then
                omitted = bounded_logadd(powerlognorm_lcdf(a,params),powerlognorm_lsf(b,params))
                value = log1p_safe(-exp(omitted))
            end if
        end if
    end function powerlognorm_interval

    pure function skewcauchy_logside(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite evaluation point
        real(dp), intent(in) :: params(:) !! skewness, location, positive scale
        real(dp) :: value
        value = log1p_safe(params(1))
        if (x < params(2)) value = log1p_safe(-params(1))
    end function skewcauchy_logside

    pure function skewcauchy_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! skewness, location, positive scale
        real(dp) :: value, logz
        logz = cauchy_log_distance(x,params(2:3))-skewcauchy_logside(x,params)
        value = -log(scifort_pi)-log(params(3))-log1pexp(2.0_dp*logz)
    end function skewcauchy_lpdf

    pure function skewcauchy_outer_tail(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite evaluation point
        real(dp), intent(in) :: params(:) !! skewness, location, positive scale
        real(dp) :: value, logh, logz
        logh = skewcauchy_logside(x,params)
        logz = cauchy_log_distance(x,params(2:3))-logh
        value = logh-log(scifort_pi)+cauchy_logatan_exp(-logz)
    end function skewcauchy_outer_tail

    pure function skewcauchy_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! skewness, location, positive scale
        real(dp) :: value, center, logh, logz
        center = log1p_safe(-params(1))-scifort_log_two
        if (params(1) < 0.0_dp) center = log1p_safe(-0.5_dp*(1.0_dp+params(1)))
        value = center
        if (x == params(2)) return
        value = skewcauchy_outer_tail(x,params)
        if (x > params(2)) then
            if (value < -scifort_log_two) then
                value = log1p_safe(-exp(value))
            else
                logh = log1p_safe(params(1))
                logz = cauchy_log_distance(x,params(2:3))-logh
                value = bounded_logadd(center,logh-log(scifort_pi)+cauchy_logatan_exp(logz))
            end if
        end if
    end function skewcauchy_lcdf

    pure function skewcauchy_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! skewness, location, positive scale
        real(dp) :: value, center, logh, logz
        center = log1p_safe(params(1))-scifort_log_two
        if (params(1) > 0.0_dp) center = log1p_safe(-0.5_dp*(1.0_dp-params(1)))
        value = center
        if (x == params(2)) return
        value = skewcauchy_outer_tail(x,params)
        if (x < params(2)) then
            if (value < -scifort_log_two) then
                value = log1p_safe(-exp(value))
            else
                logh = log1p_safe(-params(1))
                logz = cauchy_log_distance(x,params(2:3))-logh
                value = bounded_logadd(center,logh-log(scifort_pi)+cauchy_logatan_exp(logz))
            end if
        end if
    end function skewcauchy_lsf

    pure function skewcauchy_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper interval endpoint, > a
        real(dp), intent(in) :: params(:) !! skewness, location, positive scale
        real(dp) :: value, logh, loga, logb, ratio, left, right, omitted
        if (b <= params(2) .or. a >= params(2)) then
            logh = log1p_safe(params(1))
            if (b <= params(2)) logh = log1p_safe(-params(1))
            loga = cauchy_log_distance(a,params(2:3))-logh
            logb = cauchy_log_distance(b,params(2:3))-logh
            ratio = log_distance(b,a)-log(params(3))-logh-log1pexp(loga+logb)
            value = logh-log(scifort_pi)+cauchy_logatan_exp(ratio)
        else
            logh = log1p_safe(-params(1))
            loga = cauchy_log_distance(a,params(2:3))-logh
            left = logh-log(scifort_pi)+cauchy_logatan_exp(loga)
            logh = log1p_safe(params(1))
            logb = cauchy_log_distance(b,params(2:3))-logh
            right = logh-log(scifort_pi)+cauchy_logatan_exp(logb)
            value = min(0.0_dp,bounded_logadd(left,right))
            if (value > -scifort_log_two) then
                omitted = bounded_logadd(skewcauchy_lcdf(a,params),skewcauchy_lsf(b,params))
                value = log1p_safe(-exp(omitted))
            end if
        end if
    end function skewcauchy_interval


    pure function powernorm_censored_loglikelihood(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! positive finite normal power
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(3)
        params = [c,0.0_dp,1.0_dp]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (c <= 0.0_dp) return
        if (params(3) <= 0.0_dp) return
        value = censored_loglikelihood(data,params,powernorm_lpdf,powernorm_lcdf, &
            powernorm_lsf,powernorm_interval)
    end function powernorm_censored_loglikelihood

    pure function powernorm_censored_nnlf(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! positive finite normal power
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -powernorm_censored_loglikelihood(data, c, loc, scale)
    end function powernorm_censored_nnlf

    pure function powerlognorm_censored_loglikelihood(data, c, s, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! positive finite normal power
        real(dp), intent(in) :: s !! positive finite lognormal shape
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(4)
        params = [c,s,0.0_dp,1.0_dp]
        if (present(loc)) params(3) = loc
        if (present(scale)) params(4) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (min(c,s) <= 0.0_dp) return
        if (params(4) <= 0.0_dp) return
        value = censored_loglikelihood(data,params,powerlognorm_lpdf,powerlognorm_lcdf, &
            powerlognorm_lsf,powerlognorm_interval)
    end function powerlognorm_censored_loglikelihood

    pure function powerlognorm_censored_nnlf(data, c, s, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! positive finite normal power
        real(dp), intent(in) :: s !! positive finite lognormal shape
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -powerlognorm_censored_loglikelihood(data, c, s, loc, scale)
    end function powerlognorm_censored_nnlf

    pure function skewcauchy_censored_loglikelihood(data, a, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: a !! skewness, strictly between -1 and 1
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(3)
        params = [a,0.0_dp,1.0_dp]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (abs(a) >= 1.0_dp) return
        if (params(3) <= 0.0_dp) return
        value = censored_loglikelihood(data,params,skewcauchy_lpdf,skewcauchy_lcdf, &
            skewcauchy_lsf,skewcauchy_interval)
    end function skewcauchy_censored_loglikelihood

    pure function skewcauchy_censored_nnlf(data, a, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: a !! skewness, strictly between -1 and 1
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -skewcauchy_censored_loglikelihood(data, a, loc, scale)
    end function skewcauchy_censored_nnlf


    pure function invgamma_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite observation or censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, location, scale
        real(dp) :: value, logy
        value = negative_infinity(0.0_dp)
        if (x <= params(2)) return
        logy = -lognormal_logz(x,params)
        value = loggamma_lpdf(logy,[params(1),0.0_dp,1.0_dp]) &
            -log_distance(x,params(2))
    end function invgamma_lpdf

    pure function invgamma_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite observation or censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, location, scale
        real(dp) :: value, logy
        value = negative_infinity(0.0_dp)
        if (x <= params(2)) return
        logy = -lognormal_logz(x,params)
        value = loggamma_lsf(logy,[params(1),0.0_dp,1.0_dp])
    end function invgamma_lcdf

    pure function invgamma_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite observation or censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, location, scale
        real(dp) :: value, logy
        value = 0.0_dp
        if (x <= params(2)) return
        logy = -lognormal_logz(x,params)
        value = loggamma_lcdf(logy,[params(1),0.0_dp,1.0_dp])
    end function invgamma_lsf

    pure function invgamma_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper interval endpoint, greater than a
        real(dp), intent(in) :: params(:) !! positive shape, location, scale
        real(dp) :: value, logwidth, logq, logy, q, t(3), r(3), correction, omitted
        value = negative_infinity(0.0_dp)
        if (b <= params(2)) return
        if (a <= params(2)) then
            value = invgamma_lcdf(b,params)
            return
        end if
        logwidth = log_distance(b,a)
        logq = logwidth-log_distance(a,params(2))
        logy = -lognormal_logz(a,params)
        if (logq+max(log(params(1)+1.0_dp),logy) <= log(0.001_dp)) then
            q = exp(logq)
            t = [0.5_dp,0.5_dp-sqrt(3.0_dp/5.0_dp)/2.0_dp, &
                0.5_dp+sqrt(3.0_dp/5.0_dp)/2.0_dp]
            r = -(params(1)+1.0_dp)*log1p_safe(t*q)+exp(logy+logq)*t/(1.0_dp+t*q)
            correction = (4.0_dp/9.0_dp)*expm1_safe(r(1)) &
                +(5.0_dp/18.0_dp)*(expm1_safe(r(2))+expm1_safe(r(3)))
            value = invgamma_lpdf(a,params)+logwidth+log1p_safe(correction)
        else
            value = tail_interval(a,b,params,invgamma_lcdf,invgamma_lsf)
            if (value > -scifort_log_two) then
                omitted = bounded_logadd(invgamma_lcdf(a,params),invgamma_lsf(b,params))
                value = log1p_safe(-exp(omitted))
            end if
        end if
    end function invgamma_interval

    pure function maxwell_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite observation or censoring limit
        real(dp), intent(in) :: params(:) !! finite location and positive scale
        real(dp) :: value, logy
        value = negative_infinity(0.0_dp)
        if (x <= params(1)) return
        logy = 2.0_dp*lognormal_logz(x,[1.0_dp,params(1),params(2)])-scifort_log_two
        value = loggamma_lpdf(logy,[1.5_dp,0.0_dp,1.0_dp]) &
            -log_distance(x,params(1))+scifort_log_two
    end function maxwell_lpdf

    pure function maxwell_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite observation or censoring limit
        real(dp), intent(in) :: params(:) !! finite location and positive scale
        real(dp) :: value, logy
        value = negative_infinity(0.0_dp)
        if (x <= params(1)) return
        logy = 2.0_dp*lognormal_logz(x,[1.0_dp,params(1),params(2)])-scifort_log_two
        value = loggamma_lcdf(logy,[1.5_dp,0.0_dp,1.0_dp])
    end function maxwell_lcdf

    pure function maxwell_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite observation or censoring limit
        real(dp), intent(in) :: params(:) !! finite location and positive scale
        real(dp) :: value, logy
        value = 0.0_dp
        if (x <= params(1)) return
        logy = 2.0_dp*lognormal_logz(x,[1.0_dp,params(1),params(2)])-scifort_log_two
        value = loggamma_lsf(logy,[1.5_dp,0.0_dp,1.0_dp])
    end function maxwell_lsf

    pure function maxwell_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper interval endpoint, greater than a
        real(dp), intent(in) :: params(:) !! finite location and positive scale
        real(dp) :: value, logwidth, logq, logy, q, t(3), r(3), correction, omitted
        value = negative_infinity(0.0_dp)
        if (b <= params(1)) return
        if (a <= params(1)) then
            value = maxwell_lcdf(b,params)
            return
        end if
        logwidth = log_distance(b,a)
        logq = logwidth-log_distance(a,params(1))
        logy = 2.0_dp*lognormal_logz(a,[1.0_dp,params(1),params(2)])-scifort_log_two
        if (logq+max(0.0_dp,logy) <= log(0.001_dp)) then
            q = exp(logq)
            t = [0.5_dp,0.5_dp-sqrt(3.0_dp/5.0_dp)/2.0_dp, &
                0.5_dp+sqrt(3.0_dp/5.0_dp)/2.0_dp]
            r = 2.0_dp*log1p_safe(t*q)-exp(logy+logq)*t*(2.0_dp+t*q)
            correction = (4.0_dp/9.0_dp)*expm1_safe(r(1)) &
                +(5.0_dp/18.0_dp)*(expm1_safe(r(2))+expm1_safe(r(3)))
            value = maxwell_lpdf(a,params)+logwidth+log1p_safe(correction)
        else
            value = tail_interval(a,b,params,maxwell_lcdf,maxwell_lsf)
            if (value > -scifort_log_two) then
                omitted = bounded_logadd(maxwell_lcdf(a,params),maxwell_lsf(b,params))
                value = log1p_safe(-exp(omitted))
            end if
        end if
    end function maxwell_interval

    pure subroutine betaprime_coordinates(x, params, u, v, logu, logv)
        implicit none
        real(dp), intent(in) :: x !! finite point strictly above location
        real(dp), intent(in) :: params(:) !! positive shapes with finite sum, location, scale
        real(dp), intent(out) :: u !! z/(1+z), possibly rounded to zero or one
        real(dp), intent(out) :: v !! 1/(1+z), retained independently
        real(dp), intent(out) :: logu !! log transformed lower coordinate
        real(dp), intent(out) :: logv !! log transformed upper coordinate
        real(dp) :: logz, e
        logz = lognormal_logz(x,[1.0_dp,params(3),params(4)])
        if (size(params) == 5) logz = logz+params(5)
        logu = -log1pexp(-logz)
        logv = -log1pexp(logz)
        if (logz <= 0.0_dp) then
            e = exp(logz)
            u = e/(1.0_dp+e)
            v = 1.0_dp/(1.0_dp+e)
        else
            e = exp(-logz)
            u = 1.0_dp/(1.0_dp+e)
            v = e/(1.0_dp+e)
        end if
    end subroutine betaprime_coordinates

    pure function betaprime_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! positive shapes with finite sum, location, scale
        real(dp) :: value, u, v, logu, logv
        value = negative_infinity(0.0_dp)
        if (x < params(3)) return
        if (x == params(3)) then
            if (params(1) < 1.0_dp) value = positive_infinity(0.0_dp)
            if (params(1) == 1.0_dp) value = log(params(2))-log(params(4))
            if (size(params) == 5 .and. params(1) == 1.0_dp) then
                value = value+params(5)
            end if
            return
        end if
        call betaprime_coordinates(x,params,u,v,logu,logv)
        if (min(u,v) > 0.0_dp) then
            value = log_beta_kernel(params(1),params(2),u,v)
        else if (u == 0.0_dp) then
            value = betaprime_endpoint_tail(params(1),params(2),logu)+log(params(1))
        else
            value = betaprime_endpoint_tail(params(2),params(1),logv)+log(params(2))
        end if
        value = value-log_distance(x,params(3))
    end function betaprime_lpdf

    pure function betaprime_endpoint_tail(a, b, logu) result(value)
        implicit none
        real(dp), intent(in) :: a !! positive shape at a coordinate rounded to zero
        real(dp), intent(in) :: b !! positive opposite shape; a+b finite
        real(dp), intent(in) :: logu !! finite log coordinate below the subnormal range
        real(dp) :: value, normalization, factor
        if (a < 1.0_dp) then
            value = beta_small_tail(a,b,logu)
        else
            ! Combine normalization before multiplication if a*logu overflows.
            normalization = beta_small_tail(a,b,0.0_dp)
            factor = logu+normalization/a
            value = negative_infinity(0.0_dp)
            if (a > huge(1.0_dp)/abs(factor)) return
            value = a*factor
        end if
    end function betaprime_endpoint_tail

    pure subroutine betaprime_logtails(x, params, logp, logq)
        implicit none
        real(dp), intent(in) :: x !! finite censoring limit
        real(dp), intent(in) :: params(:) !! positive shapes with finite sum, location, scale
        real(dp), intent(out) :: logp !! log lower probability
        real(dp), intent(out) :: logq !! log upper probability
        real(dp) :: u, v, logu, logv, p, q
        logp = negative_infinity(0.0_dp)
        logq = 0.0_dp
        if (x <= params(3)) return
        call betaprime_coordinates(x,params,u,v,logu,logv)
        if (size(params) == 5) then
            if (max(params(1),params(2)) <= epsilon(1.0_dp)/16384.0_dp) then
                logp = -log1pexp(log(params(1))-log(params(2)))
                logq = -log1pexp(log(params(2))-log(params(1)))
                return
            end if
        end if
        if (u == 0.0_dp) then
            logp = betaprime_endpoint_tail(params(1),params(2),logu)
            logq = log_difference(0.0_dp,logp)
        else if (v == 0.0_dp) then
            logq = betaprime_endpoint_tail(params(2),params(1),logv)
            logp = log_difference(0.0_dp,logq)
        else
            call incomplete_beta_xy(params(1),params(2),u,v,p,q,logp,logq)
        end if
    end subroutine betaprime_logtails

    pure function betaprime_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! positive shapes with finite sum, location, scale
        real(dp) :: value, logq
        call betaprime_logtails(x,params,value,logq)
    end function betaprime_lcdf

    pure function betaprime_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! positive shapes with finite sum, location, scale
        real(dp) :: value, logp
        call betaprime_logtails(x,params,logp,value)
    end function betaprime_lsf

    pure function betaprime_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper interval endpoint, greater than a
        real(dp), intent(in) :: params(:) !! positive shapes with finite sum, location, scale
        real(dp) :: value, logwidth, logq, logq2, u, v, logu, logv, q, q2
        real(dp) :: t(3), r(3), correction, omitted, shape_sum
        value = negative_infinity(0.0_dp)
        if (b <= params(3)) return
        if (a <= params(3)) then
            value = betaprime_lcdf(b,params)
            return
        end if
        call betaprime_coordinates(a,params,u,v,logu,logv)
        logwidth = log_distance(b,a)
        logq = logwidth-log_distance(a,params(3))
        logq2 = logq+logu
        shape_sum = params(1)+params(2)
        if (logq+log(max(1.0_dp,abs(params(1)-1.0_dp))) <= log(0.001_dp) &
                .and. logq2+log(shape_sum) <= log(0.001_dp)) then
            q = exp(logq)
            q2 = exp(logq2)
            t = [0.5_dp,0.5_dp-sqrt(3.0_dp/5.0_dp)/2.0_dp, &
                0.5_dp+sqrt(3.0_dp/5.0_dp)/2.0_dp]
            r = (params(1)-1.0_dp)*log1p_safe(t*q)-shape_sum*log1p_safe(t*q2)
            correction = (4.0_dp/9.0_dp)*expm1_safe(r(1)) &
                +(5.0_dp/18.0_dp)*(expm1_safe(r(2))+expm1_safe(r(3)))
            value = betaprime_lpdf(a,params)+logwidth+log1p_safe(correction)
        else
            value = tail_interval(a,b,params,betaprime_lcdf,betaprime_lsf)
            if (value > -scifort_log_two) then
                omitted = bounded_logadd(betaprime_lcdf(a,params),betaprime_lsf(b,params))
                value = log1p_safe(-exp(omitted))
            end if
        end if
    end function betaprime_interval


    pure function invgamma_censored_loglikelihood(data, a, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: a !! finite positive first shape
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! finite positive scale, default one
        real(dp) :: value, params(3)
        params = [a,0.0_dp,1.0_dp]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (min(a,params(3)) <= 0.0_dp) return
        value = censored_loglikelihood(data,params,invgamma_lpdf,invgamma_lcdf, &
            invgamma_lsf,invgamma_interval)
    end function invgamma_censored_loglikelihood

    pure function invgamma_censored_nnlf(data, a, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: a !! finite positive first shape
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! finite positive scale, default one
        real(dp) :: value
        value = -invgamma_censored_loglikelihood(data, a, loc, scale)
    end function invgamma_censored_nnlf

    pure function betaprime_censored_loglikelihood(data, a, b, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: a !! finite positive first shape
        real(dp), intent(in) :: b !! finite positive second shape
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! finite positive scale, default one
        real(dp) :: value, params(4)
        params = [a,b,0.0_dp,1.0_dp]
        if (present(loc)) params(3) = loc
        if (present(scale)) params(4) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (min(a,b,params(4)) <= 0.0_dp) return
        if (0.5_dp*a+0.5_dp*b > 0.5_dp*huge(1.0_dp)) return
        value = censored_loglikelihood(data,params,betaprime_lpdf,betaprime_lcdf, &
            betaprime_lsf,betaprime_interval)
    end function betaprime_censored_loglikelihood

    pure function betaprime_censored_nnlf(data, a, b, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: a !! finite positive first shape
        real(dp), intent(in) :: b !! finite positive second shape
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! finite positive scale, default one
        real(dp) :: value
        value = -betaprime_censored_loglikelihood(data, a, b, loc, scale)
    end function betaprime_censored_nnlf

    pure function maxwell_censored_loglikelihood(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! finite positive scale, default one
        real(dp) :: value, params(2)
        params = [0.0_dp,1.0_dp]
        if (present(loc)) params(1) = loc
        if (present(scale)) params(2) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (params(2) <= 0.0_dp) return
        value = censored_loglikelihood(data,params,maxwell_lpdf,maxwell_lcdf, &
            maxwell_lsf,maxwell_interval)
    end function maxwell_censored_loglikelihood

    pure function maxwell_censored_nnlf(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! finite positive scale, default one
        real(dp) :: value
        value = -maxwell_censored_loglikelihood(data, loc, scale)
    end function maxwell_censored_nnlf


    pure function dgamma_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite observation or censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, finite location, positive scale
        real(dp) :: value, y, p(3)
        y = x
        p = params
        if (x < params(2)) then
            y = -x
            p(2) = -params(2)
        end if
        value = gamma_lpdf(y,p)-scifort_log_two
    end function dgamma_lpdf

    pure function dgamma_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite observation or censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, finite location, positive scale
        real(dp) :: value, y, p(3)
        if (x == params(2)) then
            value = -scifort_log_two
            return
        end if
        y = x
        p = params
        if (x < params(2)) then
            y = -x
            p(2) = -params(2)
        end if
        value = gamma_lsf(y,p)-scifort_log_two
        if (.not. (x < params(2))) value = log1p_safe(-exp(value))
    end function dgamma_lcdf

    pure function dgamma_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite observation or censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, finite location, positive scale
        real(dp) :: value, y, p(3)
        y = x
        p = params
        if (x < params(2)) then
            y = -x
            p(2) = -params(2)
        end if
        value = gamma_lsf(y,p)-scifort_log_two
        if (.not. (x >= params(2))) value = log1p_safe(-exp(value))
    end function dgamma_lsf

    pure function dgamma_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper endpoint, greater than a
        real(dp), intent(in) :: params(:) !! positive shape, finite location, positive scale
        real(dp) :: value, p(3), left, right, omitted
        p = [params(1),-params(2),params(3)]
        if (b <= params(2)) then
            value = gamma_interval(-b,-a,p)-scifort_log_two
        else if (a >= params(2)) then
            value = gamma_interval(a,b,params)-scifort_log_two
        else
            left = gamma_lcdf(-a,p)
            right = gamma_lcdf(b,params)
            value = bounded_logadd(left,right)-scifort_log_two
            ! Preserve both excluded tails for nearly full intervals.
            if (value > -scifort_log_two) then
                omitted = bounded_logadd(gamma_lsf(-a,p), &
                    gamma_lsf(b,params))-scifort_log_two
                value = log1p_safe(-exp(omitted))
            end if
        end if
    end function dgamma_interval

    pure function dweibull_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite observation or censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, finite location, positive scale
        real(dp) :: value, y, p(3)
        y = x
        p = params
        if (x < params(2)) then
            y = -x
            p(2) = -params(2)
        end if
        value = radial_weibull_lpdf(y,p)-scifort_log_two
    end function dweibull_lpdf

    pure function dweibull_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite observation or censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, finite location, positive scale
        real(dp) :: value, y, p(3)
        if (x == params(2)) then
            value = -scifort_log_two
            return
        end if
        y = x
        p = params
        if (x < params(2)) then
            y = -x
            p(2) = -params(2)
        end if
        value = radial_weibull_lsf(y,p)-scifort_log_two
        if (.not. (x < params(2))) value = log1p_safe(-exp(value))
    end function dweibull_lcdf

    pure function dweibull_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite observation or censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, finite location, positive scale
        real(dp) :: value, y, p(3)
        y = x
        p = params
        if (x < params(2)) then
            y = -x
            p(2) = -params(2)
        end if
        value = radial_weibull_lsf(y,p)-scifort_log_two
        if (.not. (x >= params(2))) value = log1p_safe(-exp(value))
    end function dweibull_lsf

    pure function dweibull_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper endpoint, greater than a
        real(dp), intent(in) :: params(:) !! positive shape, finite location, positive scale
        real(dp) :: value, p(3), left, right, omitted
        p = [params(1),-params(2),params(3)]
        if (b <= params(2)) then
            value = weibull_interval(-b,-a,p)-scifort_log_two
        else if (a >= params(2)) then
            value = weibull_interval(a,b,params)-scifort_log_two
        else
            left = radial_weibull_lcdf(-a,p)
            right = radial_weibull_lcdf(b,params)
            value = bounded_logadd(left,right)-scifort_log_two
            ! Preserve both excluded tails for nearly full intervals.
            if (value > -scifort_log_two) then
                omitted = bounded_logadd(radial_weibull_lsf(-a,p), &
                    radial_weibull_lsf(b,params))-scifort_log_two
                value = log1p_safe(-exp(omitted))
            end if
        end if
    end function dweibull_interval

    pure function radial_weibull_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite point at or above location
        real(dp), intent(in) :: params(:) !! positive shape, finite location, positive scale
        real(dp) :: value, logh
        if (x == params(2)) then
            value = weibull_logpdf(x,params(1),params(2),params(3))
            return
        end if
        logh = weibull_log_hazard(x,params)
        value = negative_infinity(0.0_dp)
        if (.not. ieee_is_finite(logh)) return
        if (logh > log(huge(1.0_dp))) return
        value = log(params(1))+logh-log_distance(x,params(2))-exp(logh)
    end function radial_weibull_lpdf

    pure function radial_weibull_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite point at or above location
        real(dp), intent(in) :: params(:) !! positive shape, finite location, positive scale
        real(dp) :: value
        value = negative_infinity(0.0_dp)
        if (x == params(2)) return
        value = log_exponential_mass(weibull_log_hazard(x,params))
    end function radial_weibull_lcdf

    pure function radial_weibull_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite point at or above location
        real(dp), intent(in) :: params(:) !! positive shape, finite location, positive scale
        real(dp) :: value
        value = 0.0_dp
        if (x == params(2)) return
        value = negative_hazard(weibull_log_hazard(x,params))
    end function radial_weibull_lsf

    pure function loggamma_left_kernel(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite point with standardized coordinate below log(tiny)
        real(dp), intent(in) :: params(:) !! positive shape, finite location, positive scale
        real(dp) :: value, z, product
        ! Retain c*z even if z itself overflows for a very small shape.
        z = location_scale_z(x,params(2:3))
        if (ieee_is_finite(z)) then
            product = negative_infinity(0.0_dp)
            if (params(1) <= huge(1.0_dp)/abs(z)) product = params(1)*z
        else
            product = negative_hazard(log(params(1))+log_distance(params(2),x)-log(params(3)))
        end if
        if (params(1) < 10.0_dp) then
            value = product-log_gamma_one_plus(params(1))
        else
            z = location_scale_z(x,params(2:3))
            value = negative_infinity(0.0_dp)
            if (ieee_is_finite(z)) value = gamma_small_kernel(params(1),z)
        end if
    end function loggamma_left_kernel

    pure function loggamma_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! positive shape, finite location, positive scale
        real(dp) :: value, z, u, factor, term
        integer :: k
        z = location_scale_z(x,params(2:3))
        value = negative_infinity(0.0_dp)
        if (z < log(tiny(1.0_dp))) then
            value = log(params(1))+loggamma_left_kernel(x,params)-log(params(3))
        else if (params(1) < 10.0_dp) then
            if (z > log(huge(1.0_dp))) return
            value = log(params(1))+params(1)*z-exp(z)-log_gamma_one_plus(params(1)) &
                -log(params(3))
        else
            ! Gamma saddle form retains the cancellation at the log-gamma mode.
            u = z-log(params(1))
            if (u > log(huge(1.0_dp))) return
            if (abs(u) < 0.1_dp) then
                term = (0.5_dp*u)*u
                factor = term
                do k = 3, 20
                    term = term*u/real(k,dp)
                    factor = factor+term
                end do
            else
                factor = expm1_safe(u)-u
            end if
            if (factor > 1.0_dp) then
                if (params(1) > huge(1.0_dp)/factor) return
            end if
            value = 0.5_dp*log(params(1))-scifort_log_sqrt_two_pi &
                -params(1)*factor-stirling_remainder(params(1))-log(params(3))
        end if
    end function loggamma_lpdf

    pure function loggamma_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, finite location, positive scale
        real(dp) :: value, z
        z = location_scale_z(x,params(2:3))
        if (z < log(tiny(1.0_dp))) then
            value = loggamma_left_kernel(x,params)
        else if (z > log(huge(1.0_dp))) then
            value = 0.0_dp
        else
            value = log_gammainc(params(1),exp(z))
        end if
    end function loggamma_lcdf

    pure function loggamma_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, finite location, positive scale
        real(dp) :: value, z, logp, u
        z = location_scale_z(x,params(2:3))
        if (z < log(tiny(1.0_dp))) then
            logp = loggamma_left_kernel(x,params)
            value = 0.0_dp
            if (ieee_is_finite(logp)) value = log_exponential_mass(bounded_log(-logp))
        else if (z > log(huge(1.0_dp))) then
            ! For exp(z)>huge, the endpoint expansion has negligible relative
            ! correction O(exp(z)/(exp(z)-c)**2). Retain a representable log tail.
            u = z-log(params(1))
            value = loggamma_lpdf(x,params)
            if (ieee_is_finite(value)) value = value+log(params(3)) &
                -log(params(1))-log_expm1_from_log(log(u))
        else
            value = log_gammaincc(params(1),exp(z))
        end if
    end function loggamma_lsf

    pure function loggamma_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper endpoint, greater than a
        real(dp), intent(in) :: params(:) !! positive shape, finite location, positive scale
        real(dp) :: value, z, logwidth, logscaled, d, t(3), r(3), correction, omitted
        z = location_scale_z(a,params(2:3))
        logwidth = log_distance(b,a)
        logscaled = logwidth-log(params(3))
        ! Bound c*d and exp(z)*d independently; the density varies by <0.002.
        if (logscaled+max(log(params(1)),z) <= log(0.001_dp)) then
            d = 0.0_dp
            if (z > log(tiny(1.0_dp))) d = exp(logscaled)
            t = [0.5_dp,0.5_dp-sqrt(3.0_dp/5.0_dp)/2.0_dp, &
                0.5_dp+sqrt(3.0_dp/5.0_dp)/2.0_dp]
            r = t*exp(log(params(1))+logscaled)
            if (z > log(tiny(1.0_dp))) then
                ! expm1(t*d)/d remains t when d underflows.
                if (d > 0.0_dp) then
                    r = r-exp(z+logscaled)*(expm1_safe(t*d)/d)
                else
                    r = r-t*exp(z+logscaled)
                end if
            end if
            correction = (4.0_dp/9.0_dp)*expm1_safe(r(1)) &
                +(5.0_dp/18.0_dp)*(expm1_safe(r(2))+expm1_safe(r(3)))
            value = loggamma_lpdf(a,params)+logwidth+log1p_safe(correction)
        else
            value = tail_interval(a,b,params,loggamma_lcdf,loggamma_lsf)
            if (value > -scifort_log_two) then
                omitted = bounded_logadd(loggamma_lcdf(a,params),loggamma_lsf(b,params))
                value = log1p_safe(-exp(omitted))
            end if
        end if
    end function loggamma_interval


    pure function dgamma_censored_loglikelihood(data, a, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: a !! finite positive shape
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! finite positive scale, default one
        real(dp) :: value, params(3)
        params = [a,0.0_dp,1.0_dp]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (a <= 0.0_dp .or. params(3) <= 0.0_dp) return
        value = censored_loglikelihood(data,params,dgamma_lpdf,dgamma_lcdf, &
            dgamma_lsf,dgamma_interval)
    end function dgamma_censored_loglikelihood

    pure function dgamma_censored_nnlf(data, a, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: a !! finite positive shape
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! finite positive scale, default one
        real(dp) :: value
        value = -dgamma_censored_loglikelihood(data,a,loc,scale)
    end function dgamma_censored_nnlf

    pure function dweibull_censored_loglikelihood(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! finite positive shape
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! finite positive scale, default one
        real(dp) :: value, params(3)
        params = [c,0.0_dp,1.0_dp]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (c <= 0.0_dp .or. params(3) <= 0.0_dp) return
        value = censored_loglikelihood(data,params,dweibull_lpdf,dweibull_lcdf, &
            dweibull_lsf,dweibull_interval)
    end function dweibull_censored_loglikelihood

    pure function dweibull_censored_nnlf(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! finite positive shape
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! finite positive scale, default one
        real(dp) :: value
        value = -dweibull_censored_loglikelihood(data,c,loc,scale)
    end function dweibull_censored_nnlf

    pure function loggamma_censored_loglikelihood(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! finite positive shape
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! finite positive scale, default one
        real(dp) :: value, params(3)
        params = [c,0.0_dp,1.0_dp]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (c <= 0.0_dp .or. params(3) <= 0.0_dp) return
        value = censored_loglikelihood(data,params,loggamma_lpdf,loggamma_lcdf, &
            loggamma_lsf,loggamma_interval)
    end function loggamma_censored_loglikelihood

    pure function loggamma_censored_nnlf(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! finite positive shape
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! finite positive scale, default one
        real(dp) :: value
        value = -loggamma_censored_loglikelihood(data,c,loc,scale)
    end function loggamma_censored_nnlf


    pure function bounded_logadd(a, b) result(value)
        implicit none
        real(dp), intent(in) :: a !! log first nonnegative term, possibly negative infinity
        real(dp), intent(in) :: b !! log second nonnegative term, possibly negative infinity
        real(dp) :: value, big
        big = max(a,b)
        value = big
        if (.not. ieee_is_finite(big)) return
        value = big+log1p_safe(exp(min(a,b)-big))
    end function bounded_logadd

    pure function bounded_log(x) result(value)
        implicit none
        real(dp), intent(in) :: x !! nonnegative dimensionless shape or width
        real(dp) :: value
        value = negative_infinity(0.0_dp)
        if (x > 0.0_dp) value = log(x)
    end function bounded_log

    pure function bounded_left(logz, logw, t) result(left)
        implicit none
        real(dp), intent(in) :: logz !! log coordinate strictly between zero and one
        real(dp), intent(in) :: logw !! log complementary coordinate, retained independently
        real(dp), intent(in) :: t !! breakpoint in [0,1]
        logical :: left
        left = .false.
        if (t <= 0.0_dp) return
        left = .true.
        if (t >= 1.0_dp) return
        if (t <= 0.5_dp) then
            left = logz < log(t)
        else
            left = logw > log1p_safe(-t)
        end if
    end function bounded_left

    pure function bounded_right(logz, logw, t) result(right)
        implicit none
        real(dp), intent(in) :: logz !! log coordinate strictly between zero and one
        real(dp), intent(in) :: logw !! log complementary coordinate, retained independently
        real(dp), intent(in) :: t !! breakpoint in [0,1]
        logical :: right
        right = .false.
        if (t >= 1.0_dp) return
        right = .true.
        if (t <= 0.0_dp) return
        if (t <= 0.5_dp) then
            right = logz > log(t)
        else
            right = logw < log1p_safe(-t)
        end if
    end function bounded_right

    pure function bounded_gap(t, logz, logw) result(value)
        implicit none
        real(dp), intent(in) :: t !! coordinate in (0,1], no smaller than z
        real(dp), intent(in) :: logz !! log lower coordinate, possibly negative infinity
        real(dp), intent(in) :: logw !! log complementary coordinate, retained independently
        real(dp) :: value, large, small
        if (t == 1.0_dp) then
            value = logw
        else
            if (t > 0.5_dp) then
                large = logw
                small = log1p_safe(-t)
            else
                large = log(t)
                small = logz
            end if
            ! Rounded coordinates may straddle a breakpoint by one ulp.
            value = log_difference(large,min(large,small))
        end if
    end function bounded_gap

    pure function bounded_right_gap(t, logz, logw) result(value)
        implicit none
        real(dp), intent(in) :: t !! breakpoint in [0,1), no greater than z
        real(dp), intent(in) :: logz !! log lower coordinate, retained independently
        real(dp), intent(in) :: logw !! log complementary coordinate
        real(dp) :: value
        if (t <= 0.5_dp) then
            value = log_difference(logz,min(logz,bounded_log(t)))
        else
            value = bounded_gap(1.0_dp-t,logw,logz)
        end if
    end function bounded_right_gap

    pure function trapezoid_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! plateau start/end in [0,1], location, width
        real(dp) :: value, left, right, z, w, lz, lw, c, d
        c = params(1)
        d = params(2)
        value = negative_infinity(0.0_dp)
        call beta_distances(x,params,left,right)
        if (left < 0.0_dp .or. right < 0.0_dp) return
        if (left == 0.0_dp .and. c > 0.0_dp) return
        if (right == 0.0_dp .and. d < 1.0_dp) return
        value = scifort_log_two-log1p_safe(d-c)-log(params(4))
        if (left == 0.0_dp .or. right == 0.0_dp) return
        call beta_coordinates(x,params,left,right,z,w,lz,lw)
        if (bounded_left(lz,lw,c)) value = value+lz-log(c)
        if (bounded_right(lz,lw,d)) value = value+lw-log1p_safe(-d)
    end function trapezoid_lpdf

    pure subroutine trapezoid_tails(x, params, logcdf, logsf)
        implicit none
        real(dp), intent(in) :: x !! finite censoring limit
        real(dp), intent(in) :: params(:) !! plateau start/end in [0,1], location, width
        real(dp), intent(out) :: logcdf !! log lower-tail probability
        real(dp), intent(out) :: logsf !! log upper-tail probability
        real(dp) :: left, right, z, w, lz, lw, c, d, lc, lr, lh, gap, terms
        logcdf = negative_infinity(0.0_dp)
        logsf = 0.0_dp
        call beta_distances(x,params,left,right)
        if (left <= 0.0_dp) return
        if (right <= 0.0_dp) then
            logcdf = 0.0_dp
            logsf = negative_infinity(0.0_dp)
            return
        end if
        call beta_coordinates(x,params,left,right,z,w,lz,lw)
        c = params(1)
        d = params(2)
        lc = bounded_log(c)
        lr = bounded_log(1.0_dp-d)
        lh = scifort_log_two-log1p_safe(d-c)
        if (bounded_left(lz,lw,c)) then
            logcdf = lh-scifort_log_two+2.0_dp*lz-lc
            gap = bounded_gap(c,lz,lw)
            terms = bounded_logadd(lr,scifort_log_two+bounded_log(d-c))
            terms = bounded_logadd(terms,gap+log1p_safe(exp(lz-lc)))
            logsf = lh-scifort_log_two+terms
        else if (bounded_right(lz,lw,d)) then
            logsf = lh-scifort_log_two+2.0_dp*lw-lr
            gap = bounded_right_gap(d,lz,lw)
            terms = bounded_logadd(lc,scifort_log_two+bounded_log(d-c))
            terms = bounded_logadd(terms,gap+log1p_safe(exp(lw-lr)))
            logcdf = lh-scifort_log_two+terms
        else
            logcdf = lh+lz+log1p_safe(-0.5_dp*exp(lc-lz))
            logsf = lh+lw+log1p_safe(-0.5_dp*exp(lr-lw))
        end if
        if (logcdf <= logsf) then
            logsf = log1p_safe(-exp(logcdf))
        else
            logcdf = log1p_safe(-exp(logsf))
        end if
    end subroutine trapezoid_tails

    pure function trapezoid_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite left-censoring limit
        real(dp), intent(in) :: params(:) !! plateau start/end in [0,1], location, width
        real(dp) :: value, unused
        call trapezoid_tails(x,params,value,unused)
    end function trapezoid_lcdf

    pure function trapezoid_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite right-censoring limit
        real(dp), intent(in) :: params(:) !! plateau start/end in [0,1], location, width
        real(dp) :: value, unused
        call trapezoid_tails(x,params,unused,value)
    end function trapezoid_lsf

    pure function trapezoid_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper endpoint, greater than a
        real(dp), intent(in) :: params(:) !! plateau start/end in [0,1], location, width
        real(dp) :: value, left, right, z, w, lza, lwa, lzb, lwb, ld, lc, lr, lh
        real(dp) :: omitted, deficit, gap
        value = negative_infinity(0.0_dp)
        call beta_distances(b,params,left,right)
        if (left <= 0.0_dp) return
        call beta_distances(a,params,left,right)
        if (right <= 0.0_dp) return
        if (left <= 0.0_dp) then
            value = trapezoid_lcdf(b,params)
            return
        end if
        call beta_coordinates(a,params,left,right,z,w,lza,lwa)
        call beta_distances(b,params,left,right)
        if (right <= 0.0_dp) then
            value = trapezoid_lsf(a,params)
            return
        end if
        call beta_coordinates(b,params,left,right,z,w,lzb,lwb)
        omitted = exp(trapezoid_lcdf(a,params))+exp(trapezoid_lsf(b,params))
        if (omitted <= 0.5_dp) then
            value = log1p_safe(-omitted)
            return
        end if
        ld = log_distance(b,a)-log(params(4))
        lc = bounded_log(params(1))
        lr = bounded_log(1.0_dp-params(2))
        lh = scifort_log_two-log1p_safe(params(2)-params(1))
        if (.not. bounded_right(lzb,lwb,params(1))) then
            value = lh+ld+bounded_logadd(lza,lzb)-scifort_log_two-lc
        else if (.not. bounded_left(lza,lwa,params(2))) then
            value = lh+ld+bounded_logadd(lwa,lwb)-scifort_log_two-lr
        else
            ! Integrate the constant plateau height minus the two triangular deficits.
            deficit = 0.0_dp
            if (bounded_left(lza,lwa,params(1))) then
                gap = bounded_gap(params(1),lza,lwa)
                deficit = exp(2.0_dp*gap-scifort_log_two-lc-ld)
            end if
            if (bounded_right(lzb,lwb,params(2))) then
                gap = bounded_right_gap(params(2),lzb,lwb)
                deficit = deficit+exp(2.0_dp*gap-scifort_log_two-lr-ld)
            end if
            value = lh+ld+log1p_safe(-deficit)
        end if
    end function trapezoid_interval

    pure function bradford_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! positive shape, unused slot, location, width
        real(dp) :: value, left, right, z, w, lz, lw, lc
        value = negative_infinity(0.0_dp)
        call beta_distances(x,params,left,right)
        if (left < 0.0_dp .or. right < 0.0_dp) return
        lc = log(params(1))
        value = -log(params(4))
        if (lc < log(epsilon(1.0_dp))-2.0_dp) return
        lz = negative_infinity(0.0_dp)
        if (right == 0.0_dp) lz = 0.0_dp
        if (left > 0.0_dp .and. right > 0.0_dp) then
            call beta_coordinates(x,params,left,right,z,w,lz,lw)
        end if
        value = lc-log(params(4))-burr_logsoftplus(lc)-log1pexp(lc+lz)
    end function bradford_lpdf

    pure subroutine bradford_tails(x, params, logcdf, logsf)
        implicit none
        real(dp), intent(in) :: x !! finite censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, unused slot, location, width
        real(dp), intent(out) :: logcdf !! log lower-tail probability
        real(dp), intent(out) :: logsf !! log upper-tail probability
        real(dp) :: left, right, z, w, lz, lw, lc, ln
        logcdf = negative_infinity(0.0_dp)
        logsf = 0.0_dp
        call beta_distances(x,params,left,right)
        if (left <= 0.0_dp) return
        if (right <= 0.0_dp) then
            logcdf = 0.0_dp
            logsf = negative_infinity(0.0_dp)
            return
        end if
        call beta_coordinates(x,params,left,right,z,w,lz,lw)
        lc = log(params(1))
        if (lc < log(epsilon(1.0_dp))-2.0_dp) then
            logcdf = lz
            logsf = lw
        else
            ln = burr_logsoftplus(lc)
            logcdf = burr_logsoftplus(lc+lz)-ln
            logsf = burr_logsoftplus(lc+lw-log1pexp(lc+lz))-ln
        end if
        if (logcdf <= logsf) then
            logsf = log1p_safe(-exp(logcdf))
        else
            logcdf = log1p_safe(-exp(logsf))
        end if
    end subroutine bradford_tails

    pure function bradford_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite left-censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, unused slot, location, width
        real(dp) :: value, unused
        call bradford_tails(x,params,value,unused)
    end function bradford_lcdf

    pure function bradford_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite right-censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, unused slot, location, width
        real(dp) :: value, unused
        call bradford_tails(x,params,unused,value)
    end function bradford_lsf

    pure function bradford_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper endpoint, greater than a
        real(dp), intent(in) :: params(:) !! positive shape, unused slot, location, width
        real(dp) :: value, left, right, z, w, lz, lw, lc, ld, omitted
        value = negative_infinity(0.0_dp)
        call beta_distances(b,params,left,right)
        if (left <= 0.0_dp) return
        call beta_distances(a,params,left,right)
        if (right <= 0.0_dp) return
        if (left <= 0.0_dp) then
            value = bradford_lcdf(b,params)
            return
        end if
        call beta_coordinates(a,params,left,right,z,w,lz,lw)
        call beta_distances(b,params,left,right)
        if (right <= 0.0_dp) then
            value = bradford_lsf(a,params)
            return
        end if
        omitted = exp(bradford_lcdf(a,params))+exp(bradford_lsf(b,params))
        if (omitted <= 0.5_dp) then
            value = log1p_safe(-omitted)
            return
        end if
        lc = log(params(1))
        ld = log_distance(b,a)-log(params(4))
        value = ld
        if (lc < log(epsilon(1.0_dp))-2.0_dp) return
        value = burr_logsoftplus(lc+ld-log1pexp(lc+lz))-burr_logsoftplus(lc)
    end function bradford_interval

    pure function bradford_censored_loglikelihood(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! finite positive Bradford shape
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! finite positive support width, default one
        real(dp) :: value, params(4)
        params = [c,1.0_dp,0.0_dp,1.0_dp]
        if (present(loc)) params(3) = loc
        if (present(scale)) params(4) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (c <= 0.0_dp) return
        if (.not. beta_parameters_valid([1.0_dp,1.0_dp,params(3),params(4)])) return
        value = censored_loglikelihood(data,params,bradford_lpdf,bradford_lcdf, &
            bradford_lsf,bradford_interval)
    end function bradford_censored_loglikelihood

    pure function bradford_censored_nnlf(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! finite positive Bradford shape
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! finite positive support width, default one
        real(dp) :: value
        value = -bradford_censored_loglikelihood(data,c,loc,scale)
    end function bradford_censored_nnlf

    pure function triang_censored_loglikelihood(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! finite standardized mode in [0,1]
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! finite positive support width, default one
        real(dp) :: value
        value = trapezoid_censored_loglikelihood(data,c,c,loc,scale)
    end function triang_censored_loglikelihood

    pure function triang_censored_nnlf(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! finite standardized mode in [0,1]
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! finite positive support width, default one
        real(dp) :: value
        value = -triang_censored_loglikelihood(data,c,loc,scale)
    end function triang_censored_nnlf

    pure function trapezoid_censored_loglikelihood(data, c, d, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! finite plateau start in [0,1]
        real(dp), intent(in) :: d !! finite plateau end in [c,1]
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! finite positive support width, default one
        real(dp) :: value, params(4)
        params = [c,d,0.0_dp,1.0_dp]
        if (present(loc)) params(3) = loc
        if (present(scale)) params(4) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (c < 0.0_dp .or. d < c .or. d > 1.0_dp) return
        if (.not. beta_parameters_valid([1.0_dp,1.0_dp,params(3),params(4)])) return
        value = censored_loglikelihood(data,params,trapezoid_lpdf,trapezoid_lcdf, &
            trapezoid_lsf,trapezoid_interval)
    end function trapezoid_censored_loglikelihood

    pure function trapezoid_censored_nnlf(data, c, d, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! finite plateau start in [0,1]
        real(dp), intent(in) :: d !! finite plateau end in [c,1]
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! finite positive support width, default one
        real(dp) :: value
        value = -trapezoid_censored_loglikelihood(data,c, d,loc,scale)
    end function trapezoid_censored_nnlf


    pure subroutine normal_transform_tails(logt, logcdf, logsf)
        implicit none
        real(dp), intent(in) :: logt !! log positive half-normal coordinate, possibly infinite
        real(dp), intent(out) :: logcdf !! log half-normal CDF
        real(dp), intent(out) :: logsf !! log half-normal survival probability
        real(dp) :: t
        if (logt <= log(0.01_dp)) then
            logcdf = scifort_log_two+normal_small_interval(0.0_dp,logt)
            logsf = log1p_safe(-exp(logcdf))
        else if (logt <= 0.0_dp) then
            t = exp(logt)
            logcdf = log(erf(t/sqrt(2.0_dp)))
            logsf = log1p_safe(-exp(logcdf))
        else
            logsf = negative_infinity(0.0_dp)
            if (logt <= 0.5_dp*(log(huge(1.0_dp))+scifort_log_two)) then
                logsf = scifort_log_two+normal_logsf(exp(logt))
            end if
            logcdf = log1p_safe(-exp(logsf))
        end if
    end subroutine normal_transform_tails

    pure function normal_transform_exponent(logt) result(value)
        implicit none
        real(dp), intent(in) :: logt !! log positive half-normal coordinate, possibly infinite
        real(dp) :: value
        value = negative_infinity(0.0_dp)
        if (logt > 0.5_dp*(log(huge(1.0_dp))+scifort_log_two)) return
        value = 0.0_dp
        if (logt < log(tiny(1.0_dp))) return
        value = negative_hazard(2.0_dp*logt-scifort_log_two)
    end function normal_transform_exponent

    pure function normal_transform_interval(loglo, loghi, logdelta) result(value)
        implicit none
        real(dp), intent(in) :: loglo !! log smaller half-normal coordinate
        real(dp), intent(in) :: loghi !! log larger half-normal coordinate
        real(dp), intent(in) :: logdelta !! log positive difference of the two log coordinates
        real(dp) :: value, t, logwidth
        value = negative_infinity(0.0_dp)
        if (loglo > 0.5_dp*(log(huge(1.0_dp))+scifort_log_two)) return
        t = exp(loglo)
        logwidth = loghi+log_hazard_mass(logdelta)
        ! NaN requests a tail-difference fallback when the interval is not narrow.
        value = normal_small_interval(t,logwidth)
        if (.not. ieee_is_nan(value)) value = value+scifort_log_two
    end function normal_transform_interval

    pure function moyal_logt(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite observation or censoring limit
        real(dp), intent(in) :: params(:) !! finite location and positive scale
        real(dp) :: value, z
        z = location_scale_z(x,params)
        value = -0.5_dp*z
        if (.not. ieee_is_finite(z)) then
            ! The half coordinate can be finite when the full coordinate overflows.
            ! Such overflow requires scale < huge()/2, so doubling scale is safe.
            value = -location_scale_z(x,[params(1),2.0_dp*params(2)])
        end if
    end function moyal_logt

    pure function moyal_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! finite location and positive scale
        real(dp) :: value, logt, exponent
        logt = moyal_logt(x,params)
        exponent = normal_transform_exponent(logt)
        value = negative_infinity(0.0_dp)
        if (.not. ieee_is_finite(exponent)) return
        value = logt+exponent-scifort_log_sqrt_two_pi-log(params(2))
    end function moyal_lpdf

    pure function moyal_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! finite location and positive scale
        real(dp) :: value, unused
        call normal_transform_tails(moyal_logt(x,params),unused,value)
    end function moyal_lcdf

    pure function moyal_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! finite location and positive scale
        real(dp) :: value, unused
        call normal_transform_tails(moyal_logt(x,params),value,unused)
    end function moyal_lsf

    pure function moyal_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper interval endpoint, greater than a
        real(dp), intent(in) :: params(:) !! finite location and positive scale
        real(dp) :: value, omitted, logdelta
        omitted = exp(moyal_lcdf(a,params))+exp(moyal_lsf(b,params))
        if (omitted <= 0.5_dp) then
            value = log1p_safe(-omitted)
            return
        end if
        logdelta = log_distance(b,a)-log(params(2))-scifort_log_two
        value = normal_transform_interval(moyal_logt(b,params),moyal_logt(a,params),logdelta)
        if (ieee_is_nan(value)) value = tail_interval(a,b,params,moyal_lcdf,moyal_lsf)
    end function moyal_interval

    pure function levy_logt(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite point strictly above location
        real(dp), intent(in) :: params(:) !! finite lower support endpoint and positive scale
        real(dp) :: value
        value = 0.5_dp*(log(params(2))-log_distance(x,params(1)))
    end function levy_logt

    pure function levy_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! finite lower support endpoint and positive scale
        real(dp) :: value, logt
        value = negative_infinity(0.0_dp)
        if (x <= params(1)) return
        logt = levy_logt(x,params)
        value = logt-log_distance(x,params(1))-scifort_log_sqrt_two_pi &
            +normal_transform_exponent(logt)
    end function levy_lpdf

    pure function levy_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! finite lower support endpoint and positive scale
        real(dp) :: value, unused
        value = negative_infinity(0.0_dp)
        if (x <= params(1)) return
        call normal_transform_tails(levy_logt(x,params),unused,value)
    end function levy_lcdf

    pure function levy_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! finite lower support endpoint and positive scale
        real(dp) :: value, unused
        value = 0.0_dp
        if (x <= params(1)) return
        call normal_transform_tails(levy_logt(x,params),value,unused)
    end function levy_lsf

    pure function levy_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper interval endpoint, greater than a
        real(dp), intent(in) :: params(:) !! finite lower support endpoint and positive scale
        real(dp) :: value, omitted, logdelta
        value = levy_lcdf(b,params)
        if (a <= params(1)) return
        omitted = exp(levy_lcdf(a,params))+exp(levy_lsf(b,params))
        if (omitted <= 0.5_dp) then
            value = log1p_safe(-omitted)
            return
        end if
        logdelta = burr_logsoftplus(log_distance(b,a)-log_distance(a,params(1))) &
            -scifort_log_two
        value = normal_transform_interval(levy_logt(b,params),levy_logt(a,params),logdelta)
        if (ieee_is_nan(value)) value = tail_interval(a,b,params,levy_lcdf,levy_lsf)
    end function levy_interval

    pure function levy_l_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! finite upper support endpoint and positive scale
        real(dp) :: value
        value = levy_lpdf(-x,[-params(1),params(2)])
    end function levy_l_lpdf

    pure function levy_l_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! finite upper support endpoint and positive scale
        real(dp) :: value
        value = levy_lsf(-x,[-params(1),params(2)])
    end function levy_l_lcdf

    pure function levy_l_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! finite upper support endpoint and positive scale
        real(dp) :: value
        value = levy_lcdf(-x,[-params(1),params(2)])
    end function levy_l_lsf

    pure function levy_l_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper interval endpoint, greater than a
        real(dp), intent(in) :: params(:) !! finite upper support endpoint and positive scale
        real(dp) :: value
        value = levy_interval(-b,-a,[-params(1),params(2)])
    end function levy_l_interval

    pure function moyal_censored_loglikelihood(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(2)
        params = [0.0_dp, 1.0_dp]
        if (present(loc)) params(1) = loc
        if (present(scale)) params(2) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (params(2) <= 0.0_dp) return
        value = censored_loglikelihood(data, params, moyal_lpdf, moyal_lcdf, &
            moyal_lsf, moyal_interval)
    end function moyal_censored_loglikelihood

    pure function moyal_censored_nnlf(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -moyal_censored_loglikelihood(data, loc, scale)
    end function moyal_censored_nnlf

    pure function levy_censored_loglikelihood(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(2)
        params = [0.0_dp, 1.0_dp]
        if (present(loc)) params(1) = loc
        if (present(scale)) params(2) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (params(2) <= 0.0_dp) return
        value = censored_loglikelihood(data, params, levy_lpdf, levy_lcdf, &
            levy_lsf, levy_interval)
    end function levy_censored_loglikelihood

    pure function levy_censored_nnlf(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -levy_censored_loglikelihood(data, loc, scale)
    end function levy_censored_nnlf

    pure function levy_l_censored_loglikelihood(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite upper support endpoint, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(2)
        params = [0.0_dp, 1.0_dp]
        if (present(loc)) params(1) = loc
        if (present(scale)) params(2) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (params(2) <= 0.0_dp) return
        value = censored_loglikelihood(data, params, levy_l_lpdf, levy_l_lcdf, &
            levy_l_lsf, levy_l_interval)
    end function levy_l_censored_loglikelihood

    pure function levy_l_censored_nnlf(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite upper support endpoint, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -levy_l_censored_loglikelihood(data, loc, scale)
    end function levy_l_censored_nnlf


    pure function hazard_complement(logh, nh) result(value)
        implicit none
        real(dp), intent(in) :: logh !! log positive hazard, possibly infinite
        real(dp), intent(in) :: nh !! negative hazard, possibly negative infinity
        real(dp) :: value
        value = log_hazard_mass(logh)
        if (nh <= -scifort_log_two) value = log1p_safe(-exp(nh))
    end function hazard_complement

    pure subroutine genlogistic_hazard(x, params, z, logh, nh)
        implicit none
        real(dp), intent(in) :: x !! finite observation or censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, finite location, positive scale
        real(dp), intent(out) :: z !! standardized coordinate, possibly infinite
        real(dp), intent(out) :: logh !! log(c*softplus(-z)), possibly infinite
        real(dp), intent(out) :: nh !! negative c*softplus(-z), possibly infinite
        real(dp) :: base
        z = location_scale_z(x,params(2:3))
        logh = log(params(1))+burr_logsoftplus(-z)
        if (.not. ieee_is_finite(z)) then
            if (z < 0.0_dp) then
                logh = log(params(1))+log_distance(params(2),x)-log(params(3))
            end if
        end if
        nh = negative_hazard(logh)
        if (.not. ieee_is_finite(z)) return
        base = log1pexp(-z)
        if (base < tiny(1.0_dp)) return
        if (base > 1.0_dp) then
            if (params(1) > huge(1.0_dp)/base) return
        end if
        nh = -params(1)*base
    end subroutine genlogistic_hazard

    pure function genlogistic_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! positive shape, finite location, positive scale
        real(dp) :: value, z, logh, nh
        call genlogistic_hazard(x,params,z,logh,nh)
        value = log(params(1))-log(params(3))-log1pexp(z)+nh
    end function genlogistic_lpdf

    pure function genlogistic_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, finite location, positive scale
        real(dp) :: value, z, logh
        call genlogistic_hazard(x,params,z,logh,value)
    end function genlogistic_lcdf

    pure function genlogistic_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, finite location, positive scale
        real(dp) :: value, z, logh, nh
        call genlogistic_hazard(x,params,z,logh,nh)
        value = hazard_complement(logh,nh)
    end function genlogistic_lsf

    pure function genlogistic_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper interval endpoint, greater than a
        real(dp), intent(in) :: params(:) !! positive shape, finite location, positive scale
        real(dp) :: value, za, zb, logh, anchor, nh, omitted, ld, delta, lu, diff
        omitted = exp(genlogistic_lcdf(a,params))+exp(genlogistic_lsf(b,params))
        if (omitted <= 0.5_dp) then
            value = log1p_safe(-omitted)
            return
        end if
        call genlogistic_hazard(a,params,za,logh,nh)
        call genlogistic_hazard(b,params,zb,logh,anchor)
        value = anchor
        if (.not. ieee_is_finite(anchor)) return
        ld = log_distance(b,a)-log(params(3))
        ! Difference softplus(-za)-softplus(-zb), retaining the original width.
        if (ld < log(0.5_dp)) then
            lu = ld
            if (ld >= log(epsilon(1.0_dp))-2.0_dp) lu = log(expm1_safe(exp(ld)))
            diff = burr_logsoftplus(-log1pexp(zb)+lu)
        else if (zb <= 0.0_dp) then
            diff = ld
            if (ld < log(huge(1.0_dp))) then
                delta = exp(ld)
                diff = log(delta+log1pexp(za)-log1pexp(zb))
            end if
        else if (za >= 0.0_dp) then
            lu = -log1pexp(za)+log_hazard_mass(ld)
            diff = lu
            if (lu >= log(epsilon(1.0_dp))-2.0_dp) diff = log(-log1p_safe(-exp(lu)))
        else if (ieee_is_finite(za)) then
            diff = log(log1pexp(-za)-log1pexp(-zb))
        else
            diff = log_distance(params(2),a)-log(params(3))
        end if
        value = anchor+log_hazard_mass(log(params(1))+diff)
    end function genlogistic_interval

    pure function log_expm1_from_log(logz) result(value)
        implicit none
        real(dp), intent(in) :: logz !! log positive standardized distance, possibly infinite
        real(dp) :: value, z
        value = logz
        if (logz < log(epsilon(1.0_dp))-2.0_dp) return
        value = positive_infinity(0.0_dp)
        if (logz > log(huge(1.0_dp))) return
        z = exp(logz)
        value = z+log1p_safe(-exp(-z))
        if (z < 0.5_dp) value = log(expm1_safe(z))
    end function log_expm1_from_log

    pure subroutine gompertz_hazard(x, params, z, logh, nh)
        implicit none
        real(dp), intent(in) :: x !! finite point strictly above location
        real(dp), intent(in) :: params(:) !! positive shape, finite endpoint, positive scale
        real(dp), intent(out) :: z !! nonnegative standardized coordinate, possibly infinite
        real(dp), intent(out) :: logh !! log(c*expm1(z)), possibly infinite
        real(dp), intent(out) :: nh !! negative c*expm1(z), possibly negative infinity
        real(dp) :: logz, base
        z = location_scale_z(x,params(2:3))
        logz = log_distance(x,params(2))-log(params(3))
        logh = log(params(1))+log_expm1_from_log(logz)
        ! Avoid a log/exp round trip for representable z, especially near unit mass.
        if (z >= 0.5_dp) logh = log(params(1))+z+log1p_safe(-exp(-z))
        nh = negative_hazard(logh)
        if (z <= 0.0_dp .or. z > log(huge(1.0_dp))) return
        base = expm1_safe(z)
        if (base < tiny(1.0_dp)) return
        if (base > 1.0_dp) then
            if (params(1) > huge(1.0_dp)/base) return
        end if
        nh = -params(1)*base
    end subroutine gompertz_hazard

    pure function gompertz_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! positive shape, finite endpoint, positive scale
        real(dp) :: value, z, logh, nh
        value = negative_infinity(0.0_dp)
        if (x < params(2)) return
        if (x == params(2)) then
            value = log(params(1))-log(params(3))
            return
        end if
        call gompertz_hazard(x,params,z,logh,nh)
        if (.not. ieee_is_finite(nh)) return
        value = log(params(1))-log(params(3))+z+nh
    end function gompertz_lpdf

    pure function gompertz_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, finite endpoint, positive scale
        real(dp) :: value, z, logh, nh
        value = negative_infinity(0.0_dp)
        if (x <= params(2)) return
        call gompertz_hazard(x,params,z,logh,nh)
        value = hazard_complement(logh,nh)
    end function gompertz_lcdf

    pure function gompertz_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, finite endpoint, positive scale
        real(dp) :: value, z, logh
        value = 0.0_dp
        if (x <= params(2)) return
        call gompertz_hazard(x,params,z,logh,value)
    end function gompertz_lsf

    pure function gompertz_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper interval endpoint, greater than a
        real(dp), intent(in) :: params(:) !! positive shape, finite endpoint, positive scale
        real(dp) :: value, z, logh, anchor, omitted, logwidth
        value = gompertz_lcdf(b,params)
        if (a <= params(2)) return
        omitted = exp(gompertz_lcdf(a,params))+exp(gompertz_lsf(b,params))
        if (omitted <= 0.5_dp) then
            value = log1p_safe(-omitted)
            return
        end if
        call gompertz_hazard(a,params,z,logh,anchor)
        value = anchor
        if (.not. ieee_is_finite(anchor)) return
        logwidth = log_distance(b,a)-log(params(3))
        value = anchor+log_hazard_mass(log(params(1))+z+log_expm1_from_log(logwidth))
    end function gompertz_interval

    pure function invweibull_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! positive shape, finite endpoint, positive scale
        real(dp) :: value, logh, nh
        value = negative_infinity(0.0_dp)
        if (x <= params(2)) return
        logh = -weibull_log_hazard(x,params)
        if (.not. ieee_is_finite(logh)) return
        nh = negative_hazard(logh)
        value = log(params(1))-log_distance(x,params(2))+logh+nh
    end function invweibull_lpdf

    pure function invweibull_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, finite endpoint, positive scale
        real(dp) :: value
        value = negative_infinity(0.0_dp)
        if (x <= params(2)) return
        value = negative_hazard(-weibull_log_hazard(x,params))
    end function invweibull_lcdf

    pure function invweibull_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, finite endpoint, positive scale
        real(dp) :: value, logh
        value = 0.0_dp
        if (x <= params(2)) return
        logh = -weibull_log_hazard(x,params)
        value = hazard_complement(logh,negative_hazard(logh))
    end function invweibull_lsf

    pure function invweibull_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper interval endpoint, greater than a
        real(dp), intent(in) :: params(:) !! positive shape, finite endpoint, positive scale
        real(dp) :: value, omitted, anchor, logdelta, logh
        value = invweibull_lcdf(b,params)
        if (a <= params(2)) return
        omitted = exp(invweibull_lcdf(a,params))+exp(invweibull_lsf(b,params))
        if (omitted <= 0.5_dp) then
            value = log1p_safe(-omitted)
            return
        end if
        anchor = value
        if (.not. ieee_is_finite(anchor)) return
        logdelta = log(params(1))+burr_logsoftplus( &
            log_distance(b,a)-log_distance(a,params(2)))
        logh = -weibull_log_hazard(a,params)
        value = anchor+log_hazard_mass(logh+log_hazard_mass(logdelta))
    end function invweibull_interval

    pure function genlogistic_censored_loglikelihood(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! finite positive shape
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! finite positive scale, default one
        real(dp) :: value, params(3)
        params = [c,0.0_dp,1.0_dp]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (c <= 0.0_dp .or. params(3) <= 0.0_dp) return
        value = censored_loglikelihood(data,params,genlogistic_lpdf,genlogistic_lcdf, &
            genlogistic_lsf,genlogistic_interval)
    end function genlogistic_censored_loglikelihood

    pure function genlogistic_censored_nnlf(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! finite positive shape
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! finite positive scale, default one
        real(dp) :: value
        value = -genlogistic_censored_loglikelihood(data,c,loc,scale)
    end function genlogistic_censored_nnlf

    pure function gompertz_censored_loglikelihood(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! finite positive shape
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! finite positive scale, default one
        real(dp) :: value, params(3)
        params = [c,0.0_dp,1.0_dp]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (c <= 0.0_dp .or. params(3) <= 0.0_dp) return
        value = censored_loglikelihood(data,params,gompertz_lpdf,gompertz_lcdf, &
            gompertz_lsf,gompertz_interval)
    end function gompertz_censored_loglikelihood

    pure function gompertz_censored_nnlf(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! finite positive shape
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! finite positive scale, default one
        real(dp) :: value
        value = -gompertz_censored_loglikelihood(data,c,loc,scale)
    end function gompertz_censored_nnlf

    pure function invweibull_censored_loglikelihood(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! finite positive shape
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! finite positive scale, default one
        real(dp) :: value, params(3)
        params = [c,0.0_dp,1.0_dp]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (c <= 0.0_dp .or. params(3) <= 0.0_dp) return
        value = censored_loglikelihood(data,params,invweibull_lpdf,invweibull_lcdf, &
            invweibull_lsf,invweibull_interval)
    end function invweibull_censored_loglikelihood

    pure function invweibull_censored_nnlf(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! finite positive shape
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! finite positive scale, default one
        real(dp) :: value
        value = -invweibull_censored_loglikelihood(data,c,loc,scale)
    end function invweibull_censored_nnlf


    pure function uniform_censored_loglikelihood(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! finite positive support width, default one
        real(dp) :: value
        value = powerlaw_censored_loglikelihood(data,1.0_dp,loc,scale)
    end function uniform_censored_loglikelihood

    pure function uniform_censored_nnlf(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! finite positive support width, default one
        real(dp) :: value
        value = -uniform_censored_loglikelihood(data,loc,scale)
    end function uniform_censored_nnlf

    pure function arcsine_censored_loglikelihood(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! finite positive support width, default one
        real(dp) :: value, params(4)
        params = [0.5_dp,0.5_dp,0.0_dp,1.0_dp]
        if (present(loc)) params(3) = loc
        if (present(scale)) params(4) = scale
        value = quiet_nan(0.0_dp)
        if (.not. beta_parameters_valid(params)) return
        value = censored_loglikelihood(data,params,beta_lpdf,beta_lcdf,beta_lsf,arcsine_interval)
    end function arcsine_censored_loglikelihood

    pure function arcsine_censored_nnlf(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! finite positive support width, default one
        real(dp) :: value
        value = -arcsine_censored_loglikelihood(data,loc,scale)
    end function arcsine_censored_nnlf

    pure function arcsine_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper endpoint, greater than a
        real(dp), intent(in) :: params(:) !! half beta shapes, finite endpoint, positive width
        real(dp) :: value, omitted
        omitted = exp(beta_lcdf(a,params))+exp(beta_lsf(b,params))
        if (omitted <= 0.5_dp) then
            value = log1p_safe(-omitted)
        else
            value = beta_interval(a,b,params)
        end if
    end function arcsine_interval

    pure function powerlaw_censored_loglikelihood(data, a, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: a !! finite positive power-function shape
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! finite positive support width, default one
        real(dp) :: value, params(4)
        params = [a,1.0_dp,0.0_dp,1.0_dp]
        if (present(loc)) params(3) = loc
        if (present(scale)) params(4) = scale
        value = quiet_nan(0.0_dp)
        if (.not. beta_parameters_valid(params)) return
        value = censored_loglikelihood(data,params,powerlaw_lpdf,powerlaw_lcdf, &
            powerlaw_lsf,powerlaw_interval)
    end function powerlaw_censored_loglikelihood

    pure function powerlaw_censored_nnlf(data, a, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: a !! finite positive power-function shape
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! finite positive support width, default one
        real(dp) :: value
        value = -powerlaw_censored_loglikelihood(data,a,loc,scale)
    end function powerlaw_censored_nnlf

    pure function powerlaw_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! positive shape, unit beta shape, location, width
        real(dp) :: value, left, right, z, w, logz, logw, factor
        value = negative_infinity(0.0_dp)
        call beta_distances(x,params,left,right)
        if (left < 0.0_dp .or. right < 0.0_dp) return
        if (left == 0.0_dp) then
            if (params(1) < 1.0_dp) value = positive_infinity(0.0_dp)
            if (params(1) == 1.0_dp) value = -log(params(4))
            return
        end if
        if (right == 0.0_dp .or. params(1) == 1.0_dp) then
            value = log(params(1))-log(params(4))
            return
        end if
        call beta_coordinates(x,params,left,right,z,w,logz,logw)
        factor = params(1)-1.0_dp
        if (factor > 1.0_dp) then
            if (-logz > huge(1.0_dp)/factor) return
        end if
        value = log(params(1))-log(params(4))+factor*logz
    end function powerlaw_lpdf

    pure subroutine powerlaw_hazard(x, params, nh, logh)
        implicit none
        real(dp), intent(in) :: x !! finite point strictly inside support
        real(dp), intent(in) :: params(:) !! positive shape, unit beta shape, location, width
        real(dp), intent(out) :: nh !! log CDF, possibly negative infinity
        real(dp), intent(out) :: logh !! log of negative log CDF
        real(dp) :: left, right, z, w, logz, logw
        call beta_distances(x,params,left,right)
        call beta_coordinates(x,params,left,right,z,w,logz,logw)
        if (w < tiny(1.0_dp)) then
            ! -log(z)=w+O(w*w); retain log(w) when division rounds to zero
            ! or to a subnormal with poor relative precision.
            logh = log(params(1))+logw
            nh = negative_hazard(logh)
            return
        end if
        logh = log(params(1))+log(-logz)
        nh = negative_infinity(0.0_dp)
        if (params(1) > 1.0_dp) then
            if (-logz > huge(1.0_dp)/params(1)) return
        end if
        nh = params(1)*logz
    end subroutine powerlaw_hazard

    pure function powerlaw_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite left censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, unit beta shape, location, width
        real(dp) :: value, left, right, logh
        value = negative_infinity(0.0_dp)
        call beta_distances(x,params,left,right)
        if (left <= 0.0_dp) return
        value = 0.0_dp
        if (right <= 0.0_dp) return
        call powerlaw_hazard(x,params,value,logh)
    end function powerlaw_lcdf

    pure function powerlaw_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite right censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, unit beta shape, location, width
        real(dp) :: value, left, right, nh, logh
        value = 0.0_dp
        call beta_distances(x,params,left,right)
        if (left <= 0.0_dp) return
        value = negative_infinity(0.0_dp)
        if (right <= 0.0_dp) return
        call powerlaw_hazard(x,params,nh,logh)
        if (nh <= -scifort_log_two) then
            value = log1p_safe(-exp(nh))
        else
            value = log_hazard_mass(logh)
        end if
    end function powerlaw_lsf

    pure function powerlaw_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper endpoint, greater than a
        real(dp), intent(in) :: params(:) !! positive shape, unit beta shape, location, width
        real(dp) :: value, left, right, logr, logdelta, omitted
        value = negative_infinity(0.0_dp)
        call beta_distances(b,params,left,right)
        if (left <= 0.0_dp) return
        call beta_distances(a,params,left,right)
        if (right <= 0.0_dp) return
        if (left <= 0.0_dp) then
            value = powerlaw_lcdf(b,params)
            return
        end if
        call beta_distances(b,params,left,right)
        if (right <= 0.0_dp) then
            value = powerlaw_lsf(a,params)
            return
        end if
        omitted = exp(powerlaw_lcdf(a,params))+exp(powerlaw_lsf(b,params))
        if (omitted <= 0.5_dp) then
            value = log1p_safe(-omitted)
        else if (params(1) == 1.0_dp) then
            value = log_distance(b,a)-log(params(4))
        else
            logr = log_distance(b,a)-log_distance(a,params(3))
            logdelta = log(params(1))+burr_logsoftplus(logr)
            value = powerlaw_lcdf(b,params)+log_hazard_mass(logdelta)
        end if
    end function powerlaw_interval


    pure function fisk_censored_loglikelihood(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! finite positive shape
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! finite positive scale, default one
        real(dp) :: value
        value = burr12_censored_loglikelihood(data,c,1.0_dp,loc,scale)
    end function fisk_censored_loglikelihood

    pure function fisk_censored_nnlf(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! finite positive shape
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! finite positive scale, default one
        real(dp) :: value
        value = -fisk_censored_loglikelihood(data,c,loc,scale)
    end function fisk_censored_nnlf

    pure function burr_censored_loglikelihood(data, c, d, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! finite positive first shape
        real(dp), intent(in) :: d !! finite positive second shape
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! finite positive scale, default one
        real(dp) :: value, params(4)
        params = [c,d,0.0_dp,1.0_dp]
        if (present(loc)) params(3) = loc
        if (present(scale)) params(4) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (c <= 0.0_dp .or. d <= 0.0_dp .or. params(4) <= 0.0_dp) return
        value = censored_loglikelihood(data,params,burr_lpdf,burr_lcdf, &
            burr_lsf,burr_interval)
    end function burr_censored_loglikelihood

    pure function burr_censored_nnlf(data, c, d, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! finite positive first shape
        real(dp), intent(in) :: d !! finite positive second shape
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! finite positive scale, default one
        real(dp) :: value
        value = -burr_censored_loglikelihood(data,c,d,loc,scale)
    end function burr_censored_nnlf

    pure function burr12_censored_loglikelihood(data, c, d, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! finite positive first shape
        real(dp), intent(in) :: d !! finite positive second shape
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! finite positive scale, default one
        real(dp) :: value, params(4)
        params = [c,d,0.0_dp,1.0_dp]
        if (present(loc)) params(3) = loc
        if (present(scale)) params(4) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (c <= 0.0_dp .or. d <= 0.0_dp .or. params(4) <= 0.0_dp) return
        value = censored_loglikelihood(data,params,burr12_lpdf,burr12_lcdf, &
            burr12_lsf,burr12_interval)
    end function burr12_censored_loglikelihood

    pure function burr12_censored_nnlf(data, c, d, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! finite positive first shape
        real(dp), intent(in) :: d !! finite positive second shape
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! finite positive scale, default one
        real(dp) :: value
        value = -burr12_censored_loglikelihood(data,c,d,loc,scale)
    end function burr12_censored_nnlf

    pure function burr_logsoftplus(x) result(value)
        implicit none
        real(dp), intent(in) :: x !! log argument of log1p(exp(x)), possibly infinite
        real(dp) :: value
        ! For exp(x) below roundoff, log(log1p(exp(x))) = x + O(exp(x)).
        value = x
        if (x >= log(epsilon(1.0_dp))-2.0_dp) value = log(log1pexp(x))
    end function burr_logsoftplus

    pure subroutine burr_hazard(x, params, reverse, q, logh, nh)
        implicit none
        real(dp), intent(in) :: x !! finite point strictly above location
        real(dp), intent(in) :: params(:) !! c, d, finite endpoint, positive scale
        logical, intent(in) :: reverse !! true for Type III, false for Type XII
        real(dp), intent(out) :: q !! signed c*log(z), possibly infinite
        real(dp), intent(out) :: logh !! log(d*log1p(exp(q))), possibly infinite
        real(dp), intent(out) :: nh !! negative d*log1p(exp(q)), possibly infinite
        real(dp) :: v, logbase, base, product
        v = lognormal_logz(x,[params(1),params(3),params(4)])
        q = weibull_log_hazard(x,[params(1),params(3),params(4)])
        if (reverse) q = -q
        if (ieee_is_finite(q)) then
            logbase = burr_logsoftplus(q)
        else if (q > 0.0_dp) then
            ! d can compensate c even when c*log(z) exceeds huge().
            logbase = log(params(1))+log(abs(v))
        else
            logbase = negative_infinity(0.0_dp)
        end if
        logh = log(params(2))+logbase
        nh = negative_hazard(logh)
        if (.not. ieee_is_finite(q)) then
            if (q > 0.0_dp) then
                ! Here c is large and abs(log(z)) > 1. Group c*d first.
                if (params(2) <= huge(1.0_dp)/params(1)) then
                    product = params(1)*params(2)
                    logh = log(product)+log(abs(v))
                    nh = negative_infinity(0.0_dp)
                    if (product <= huge(1.0_dp)/abs(v)) nh = -product*abs(v)
                end if
            end if
            return
        end if
        base = log1pexp(q)
        if (base < tiny(1.0_dp)) return
        if (base > 1.0_dp) then
            if (params(2) > huge(1.0_dp)/base) return
        end if
        ! Direct multiplication avoids amplifying log/exp rounding near unit mass.
        nh = -params(2)*base
    end subroutine burr_hazard

    pure function burr_density(x, params, reverse) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! c, d, finite endpoint, positive scale
        logical, intent(in) :: reverse !! true for Type III, false for Type XII
        real(dp) :: value, exponent, q, logh, nh
        value = negative_infinity(0.0_dp)
        if (x < params(3)) return
        if (x == params(3)) then
            exponent = params(1)
            if (reverse) then
                if (params(1) > 1.0_dp) then
                    if (params(2) > huge(1.0_dp)/params(1)) return
                end if
                exponent = params(1)*params(2)
            end if
            if (exponent < 1.0_dp) value = positive_infinity(0.0_dp)
            if (exponent == 1.0_dp) then
                value = log(params(1))+log(params(2))-log(params(4))
            end if
            return
        end if
        call burr_hazard(x,params,reverse,q,logh,nh)
        value = log(params(1))+log(params(2))-log_distance(x,params(3)) &
            -log1pexp(-q)+nh
    end function burr_density

    pure function burr_tail(x, params, reverse, complement) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite censoring limit
        real(dp), intent(in) :: params(:) !! c, d, finite endpoint, positive scale
        logical, intent(in) :: reverse !! true for Type III, false for Type XII
        logical, intent(in) :: complement !! true to complement the direct exp(-hazard) tail
        real(dp) :: value, q, logh, nh
        value = 0.0_dp
        if (reverse .neqv. complement) value = negative_infinity(0.0_dp)
        if (x <= params(3)) return
        call burr_hazard(x,params,reverse,q,logh,nh)
        value = nh
        if (complement) then
            if (nh <= -scifort_log_two) then
                value = log1p_safe(-exp(nh))
            else
                value = log_hazard_mass(logh)
            end if
        end if
    end function burr_tail

    pure function burr_mass(a, b, params, reverse) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper endpoint, greater than a
        real(dp), intent(in) :: params(:) !! c, d, finite endpoint, positive scale
        logical, intent(in) :: reverse !! true for Type III, false for Type XII
        real(dp) :: value, qa, qb, qlo, qhi, la, lb, logdelta, logr, span, delta
        real(dp) :: logdiff, logu, omitted, anchor, unused, nh
        value = burr_tail(b,params,reverse,.not. reverse)
        if (a <= params(3)) return
        omitted = exp(burr_tail(a,params,reverse,.not. reverse)) &
            +exp(burr_tail(b,params,reverse,reverse))
        if (omitted <= 0.5_dp) then
            value = log1p_safe(-omitted)
            return
        end if
        call burr_hazard(a,params,reverse,qa,unused,nh)
        anchor = nh
        call burr_hazard(b,params,reverse,qb,unused,nh)
        if (reverse) anchor = nh
        qlo = min(qa,qb)
        qhi = max(qa,qb)
        logr = log_distance(b,a)-log_distance(a,params(3))
        logdelta = log(params(1))+burr_logsoftplus(logr)
        if (logdelta < log(0.5_dp)) then
            ! softplus(qhi)-softplus(qlo)=softplus(log(sigmoid(qlo)*expm1(delta))).
            logu = logdelta
            if (logdelta >= log(epsilon(1.0_dp))-2.0_dp) then
                delta = exp(logdelta)
                logu = log(expm1_safe(delta))
            end if
            logdiff = burr_logsoftplus(-log1pexp(-qlo)+logu)
        else if (qlo >= 0.0_dp) then
            logdiff = logdelta
            if (logdelta < log(huge(1.0_dp))) then
                delta = exp(logdelta)
                span = delta+log1pexp(-qhi)-log1pexp(-qlo)
                logdiff = log(span)
            end if
        else if (qhi <= 0.0_dp) then
            ! 1-exp(-difference)=sigmoid(qhi)*(1-exp(-delta)).
            logu = -log1pexp(-qhi)+log_hazard_mass(logdelta)
            logdiff = logu
            if (logu >= log(epsilon(1.0_dp))-2.0_dp) then
                logdiff = log(-log1p_safe(-exp(logu)))
            end if
        else if (ieee_is_finite(qhi)) then
            logdiff = log(log1pexp(qhi)-log1pexp(qlo))
        else
            ! Crossing interval with an unrepresentably large positive endpoint.
            la = abs(lognormal_logz(a,[params(1),params(3),params(4)]))
            lb = abs(lognormal_logz(b,[params(1),params(3),params(4)]))
            span = lb
            if (reverse) span = la
            logdiff = log(params(1))+log(span)
        end if
        value = anchor+log_hazard_mass(log(params(2))+logdiff)
    end function burr_mass

    pure function burr_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite observation or censoring limit
        real(dp), intent(in) :: params(:) !! c, d, finite endpoint, positive scale
        real(dp) :: value
        value = burr_density(x,params,.true.)
    end function burr_lpdf

    pure function burr_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite observation or censoring limit
        real(dp), intent(in) :: params(:) !! c, d, finite endpoint, positive scale
        real(dp) :: value
        value = burr_tail(x,params,.true.,.false.)
    end function burr_lcdf

    pure function burr_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite observation or censoring limit
        real(dp), intent(in) :: params(:) !! c, d, finite endpoint, positive scale
        real(dp) :: value
        value = burr_tail(x,params,.true.,.true.)
    end function burr_lsf

    pure function burr_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper endpoint, greater than a
        real(dp), intent(in) :: params(:) !! c, d, finite endpoint, positive scale
        real(dp) :: value
        value = burr_mass(a,b,params,.true.)
    end function burr_interval

    pure function burr12_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite observation or censoring limit
        real(dp), intent(in) :: params(:) !! c, d, finite endpoint, positive scale
        real(dp) :: value
        value = burr_density(x,params,.false.)
    end function burr12_lpdf

    pure function burr12_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite observation or censoring limit
        real(dp), intent(in) :: params(:) !! c, d, finite endpoint, positive scale
        real(dp) :: value
        value = burr_tail(x,params,.false.,.true.)
    end function burr12_lcdf

    pure function burr12_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite observation or censoring limit
        real(dp), intent(in) :: params(:) !! c, d, finite endpoint, positive scale
        real(dp) :: value
        value = burr_tail(x,params,.false.,.false.)
    end function burr12_lsf

    pure function burr12_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper endpoint, greater than a
        real(dp), intent(in) :: params(:) !! c, d, finite endpoint, positive scale
        real(dp) :: value
        value = burr_mass(a,b,params,.false.)
    end function burr12_interval


    pure function lomax_censored_loglikelihood(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! finite positive shape
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! finite positive scale, default one
        real(dp) :: value, params(3)
        params = [c, 0.0_dp, 1.0_dp]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (c <= 0.0_dp .or. params(3) <= 0.0_dp) return
        value = censored_loglikelihood(data,params,lomax_lpdf,lomax_lcdf, &
            lomax_lsf,lomax_interval)
    end function lomax_censored_loglikelihood

    pure function lomax_censored_nnlf(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! finite positive shape
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! finite positive scale, default one
        real(dp) :: value
        value = -lomax_censored_loglikelihood(data,c,loc,scale)
    end function lomax_censored_nnlf

    pure function lomax_base(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite point strictly above location
        real(dp), intent(in) :: params(:) !! positive shape, finite endpoint, positive scale
        real(dp) :: value, z
        z = location_scale_z(x,params(2:3))
        if (ieee_is_finite(z)) then
            value = log1p_safe(z)
        else
            value = log1pexp(lognormal_logz(x,params))
        end if
    end function lomax_base

    pure function lomax_negative_hazard(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite point strictly above location
        real(dp), intent(in) :: params(:) !! positive shape, finite endpoint, positive scale
        real(dp) :: value, base
        base = lomax_base(x,params)
        if (base < tiny(1.0_dp)) then
            ! Subnormal displacements lose relative precision; shape can amplify it.
            value = negative_hazard(log(params(1))+lognormal_logz(x,params))
        else
            value = negative_infinity(0.0_dp)
            if (base > 1.0_dp) then
                if (params(1) > huge(1.0_dp)/base) return
            end if
            value = -params(1)*base
        end if
    end function lomax_negative_hazard

    pure function lomax_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! positive shape, finite endpoint, positive scale
        real(dp) :: value
        value = negative_infinity(0.0_dp)
        if (x < params(2)) return
        value = log(params(1))-log(params(3))
        if (x == params(2)) return
        value = value-lomax_base(x,params)+lomax_negative_hazard(x,params)
    end function lomax_lpdf

    pure function lomax_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite left censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, finite endpoint, positive scale
        real(dp) :: value, h, base, logh
        value = negative_infinity(0.0_dp)
        if (x <= params(2)) return
        h = lomax_negative_hazard(x,params)
        if (h <= -scifort_log_two) then
            value = log1p_safe(-exp(h))
        else
            base = lomax_base(x,params)
            if (base < tiny(1.0_dp)) then
                logh = log(params(1))+lognormal_logz(x,params)
            else
                logh = log(params(1))+log(base)
            end if
            value = log_hazard_mass(logh)
        end if
    end function lomax_lcdf

    pure function lomax_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite right censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, finite endpoint, positive scale
        real(dp) :: value
        value = 0.0_dp
        if (x <= params(2)) return
        value = lomax_negative_hazard(x,params)
    end function lomax_lsf

    pure function lomax_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper interval endpoint, greater than a
        real(dp), intent(in) :: params(:) !! positive shape, finite endpoint, positive scale
        real(dp) :: value, logs, logd, logdenom, logr, logspan
        value = lomax_lcdf(b,params)
        if (a <= params(2)) return
        ! Denominator scale+a-loc is a positive sum which may exceed huge().
        logs = log(params(3))
        logd = log_distance(a,params(2))
        logdenom = max(logs,logd)+log1pexp(-abs(logs-logd))
        logr = log_distance(b,a)-logdenom
        logspan = logr
        if (logr >= log(epsilon(1.0_dp))-2.0_dp) logspan = log(log1pexp(logr))
        value = lomax_lsf(a,params)+log_hazard_mass(log(params(1))+logspan)
    end function lomax_interval

    pure function pareto_censored_loglikelihood(data, b, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: b !! finite positive shape
        real(dp), intent(in), optional :: loc !! finite location shift, default zero
        real(dp), intent(in), optional :: scale !! finite positive scale, default one
        real(dp) :: value, params(3)
        params = [b, 0.0_dp, 1.0_dp]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (b <= 0.0_dp .or. params(3) <= 0.0_dp) return
        value = censored_loglikelihood(data,params,pareto_lpdf,pareto_lcdf, &
            pareto_lsf,pareto_interval)
    end function pareto_censored_loglikelihood

    pure function pareto_censored_nnlf(data, b, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: b !! finite positive shape
        real(dp), intent(in), optional :: loc !! finite location shift, default zero
        real(dp), intent(in), optional :: scale !! finite positive scale, default one
        real(dp) :: value
        value = -pareto_censored_loglikelihood(data,b,loc,scale)
    end function pareto_censored_nnlf

    pure function pareto_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! positive shape, finite location, positive scale
        real(dp) :: value, logz
        value = negative_infinity(0.0_dp)
        if (x <= params(2)) return
        logz = lognormal_logz(x,params)
        if (logz < 0.0_dp) return
        value = log(params(1))-log_distance(x,params(2))-weibull_log_hazard(x,params)
    end function pareto_lpdf

    pure function pareto_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite left censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, finite location, positive scale
        real(dp) :: value, logz
        value = negative_infinity(0.0_dp)
        if (x <= params(2)) return
        logz = lognormal_logz(x,params)
        if (logz <= 0.0_dp) return
        ! Keep small hazards in log form, but avoid a log/exp round trip for
        ! large hazards: its error would be amplified in a nearly unit CDF.
        value = weibull_log_hazard(x,params)
        if (value >= scifort_log_two) then
            value = log1p_safe(-exp(-value))
        else
            value = log_hazard_mass(log(params(1))+log(logz))
        end if
    end function pareto_lcdf

    pure function pareto_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite right censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, finite location, positive scale
        real(dp) :: value
        value = 0.0_dp
        if (x <= params(2)) return
        value = -max(0.0_dp,weibull_log_hazard(x,params))
    end function pareto_lsf

    pure function pareto_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper interval endpoint, greater than a
        real(dp), intent(in) :: params(:) !! positive shape, finite location, positive scale
        real(dp) :: value, logr, logspan
        ! Clip intervals crossing the support boundary without forming loc+scale.
        value = pareto_lcdf(b,params)
        if (a <= params(2)) return
        if (lognormal_logz(a,params) <= 0.0_dp) return
        logr = log_distance(b,a)-log_distance(a,params(2))
        logspan = logr
        if (logr >= log(epsilon(1.0_dp))-2.0_dp) logspan = log(log1pexp(logr))
        value = pareto_lsf(a,params)+log_hazard_mass(log(params(1))+logspan)
    end function pareto_interval

    pure function loglaplace_censored_loglikelihood(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! finite positive shape
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! finite positive scale, default one
        real(dp) :: value, params(3)
        params = [c, 0.0_dp, 1.0_dp]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (c <= 0.0_dp .or. params(3) <= 0.0_dp) return
        value = censored_loglikelihood(data,params,loglaplace_lpdf,loglaplace_lcdf, &
            loglaplace_lsf,loglaplace_interval)
    end function loglaplace_censored_loglikelihood

    pure function loglaplace_censored_nnlf(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! finite positive shape
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! finite positive scale, default one
        real(dp) :: value
        value = -loglaplace_censored_loglikelihood(data,c,loc,scale)
    end function loglaplace_censored_nnlf

    pure function loglaplace_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! positive shape, finite endpoint, positive scale
        real(dp) :: value
        value = negative_infinity(0.0_dp)
        if (x < params(2)) return
        if (x == params(2)) then
            if (params(1) < 1.0_dp) value = positive_infinity(0.0_dp)
            if (params(1) == 1.0_dp) value = -scifort_log_two-log(params(3))
            return
        end if
        value = log(params(1))-scifort_log_two-log_distance(x,params(2)) &
            -abs(weibull_log_hazard(x,params))
    end function loglaplace_lpdf

    pure function loglaplace_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite left censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, finite endpoint, positive scale
        real(dp) :: value, z
        value = negative_infinity(0.0_dp)
        if (x <= params(2)) return
        z = weibull_log_hazard(x,params)
        if (z <= 0.0_dp) then
            value = z-scifort_log_two
        else
            value = log1p_safe(-0.5_dp*exp(-z))
        end if
    end function loglaplace_lcdf

    pure function loglaplace_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite right censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, finite endpoint, positive scale
        real(dp) :: value, z
        value = 0.0_dp
        if (x <= params(2)) return
        z = weibull_log_hazard(x,params)
        if (z >= 0.0_dp) then
            value = -z-scifort_log_two
        else
            value = log1p_safe(-0.5_dp*exp(z))
        end if
    end function loglaplace_lsf

    pure function loglaplace_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper interval endpoint, greater than a
        real(dp), intent(in) :: params(:) !! positive shape, finite endpoint, positive scale
        real(dp) :: value, za, zb, logr, logspan, omitted, left, right, large, small
        value = negative_infinity(0.0_dp)
        if (b <= params(2)) return
        if (a <= params(2)) then
            value = loglaplace_lcdf(b,params)
            return
        end if
        za = lognormal_logz(a,params)
        zb = lognormal_logz(b,params)
        ! Work with the log of c*log((b-loc)/(a-loc)); retain adjacent widths.
        logr = log_distance(b,a)-log_distance(a,params(2))
        logspan = logr
        if (logr >= log(epsilon(1.0_dp))-2.0_dp) logspan = log(log1pexp(logr))
        if (zb <= 0.0_dp) then
            value = loglaplace_lcdf(b,params)+log_hazard_mass(log(params(1))+logspan)
        else if (za >= 0.0_dp) then
            value = loglaplace_lsf(a,params)+log_hazard_mass(log(params(1))+logspan)
        else
            omitted = exp(loglaplace_lcdf(a,params))+exp(loglaplace_lsf(b,params))
            if (omitted <= 0.5_dp) then
                value = log1p_safe(-omitted)
            else
                left = -scifort_log_two+log_hazard_mass(log(params(1))+log(-za))
                right = -scifort_log_two+log_hazard_mass(log(params(1))+log(zb))
                large = max(left,right)
                small = min(left,right)
                value = large+log1p_safe(exp(small-large))
            end if
        end if
    end function loglaplace_interval


    pure function laplace_asymmetric_censored_loglikelihood(data, kappa, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: kappa !! finite positive asymmetry parameter
        real(dp), intent(in), optional :: loc !! finite cusp location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(3)
        params = [kappa, 0.0_dp, 1.0_dp]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (kappa <= 0.0_dp .or. params(3) <= 0.0_dp) return
        value = censored_loglikelihood(data, params, laplace_asymmetric_lpdf, &
            laplace_asymmetric_lcdf, laplace_asymmetric_lsf, laplace_asymmetric_interval)
    end function laplace_asymmetric_censored_loglikelihood

    pure function laplace_asymmetric_censored_nnlf(data, kappa, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: kappa !! finite positive asymmetry parameter
        real(dp), intent(in), optional :: loc !! finite cusp location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -laplace_asymmetric_censored_loglikelihood(data, kappa, loc, scale)
    end function laplace_asymmetric_censored_nnlf

    pure function laplace_asymmetric_log_hazard(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite point on either side of the cusp
        real(dp), intent(in) :: params(:) !! positive shape, finite cusp, positive scale
        real(dp) :: value, t, logs
        value = negative_infinity(0.0_dp)
        if (x == params(2)) return
        t = log(params(1))
        logs = log(params(3))
        value = log_distance(max(x,params(2)),min(x,params(2)))
        if (x > params(2)) then
            value = value+(t-logs)
        else
            value = value-(t+logs)
        end if
    end function laplace_asymmetric_log_hazard

    pure function laplace_asymmetric_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! positive shape, finite cusp, positive scale
        real(dp) :: value, t
        t = abs(log(params(1)))
        value = negative_hazard(laplace_asymmetric_log_hazard(x,params)) &
            -t-log1p_safe(exp(-t)**2)-log(params(3))
    end function laplace_asymmetric_lpdf

    pure function laplace_asymmetric_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, finite cusp, positive scale
        real(dp) :: value, t, logh, left, right, large, small
        t = log(params(1))
        logh = laplace_asymmetric_log_hazard(x,params)
        left = -log1pexp(-2.0_dp*t)
        right = -log1pexp(2.0_dp*t)
        if (x <= params(2)) then
            value = left+negative_hazard(logh)
        else
            value = right+negative_hazard(logh)
            if (value <= -scifort_log_two) then
                value = log1p_safe(-exp(value))
            else
                ! F = left weight + right weight*(1-exp(-H)); retain underflowed F.
                right = right+log_hazard_mass(logh)
                large = max(left,right)
                small = min(left,right)
                value = large+log1p_safe(exp(small-large))
            end if
        end if
    end function laplace_asymmetric_lcdf

    pure function laplace_asymmetric_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, finite cusp, positive scale
        real(dp) :: value, t, logh, left, right, large, small
        t = log(params(1))
        logh = laplace_asymmetric_log_hazard(x,params)
        left = -log1pexp(-2.0_dp*t)
        right = -log1pexp(2.0_dp*t)
        if (x >= params(2)) then
            value = right+negative_hazard(logh)
        else
            value = left+negative_hazard(logh)
            if (value <= -scifort_log_two) then
                value = log1p_safe(-exp(value))
            else
                ! S = right weight + left weight*(1-exp(-H)); retain underflowed S.
                left = left+log_hazard_mass(logh)
                large = max(left,right)
                small = min(left,right)
                value = large+log1p_safe(exp(small-large))
            end if
        end if
    end function laplace_asymmetric_lsf

    pure function laplace_asymmetric_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper interval endpoint, greater than a
        real(dp), intent(in) :: params(:) !! positive shape, finite cusp, positive scale
        real(dp) :: value, t, logs, left, right, large, small, omitted
        t = log(params(1))
        logs = log(params(3))
        if (b <= params(2)) then
            value = laplace_asymmetric_lcdf(b,params) &
                + log_hazard_mass(log_distance(b,a)-(logs+t))
        else if (a >= params(2)) then
            value = laplace_asymmetric_lsf(a,params) &
                + log_hazard_mass(log_distance(b,a)+(t-logs))
        else
            omitted = exp(laplace_asymmetric_lcdf(a,params)) &
                + exp(laplace_asymmetric_lsf(b,params))
            if (omitted <= 0.5_dp) then
                value = log1p_safe(-omitted)
            else
                ! Add positive masses from each side rather than subtracting near-unit tails.
                left = -log1pexp(-2.0_dp*t) &
                    + log_hazard_mass(log_distance(params(2),a)-(logs+t))
                right = -log1pexp(2.0_dp*t) &
                    + log_hazard_mass(log_distance(b,params(2))+(t-logs))
                large = max(left,right)
                small = min(left,right)
                value = large+log1p_safe(exp(small-large))
            end if
        end if
    end function laplace_asymmetric_interval

    pure function hypsecant_censored_loglikelihood(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite center, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(2)
        params = [0.0_dp, 1.0_dp]
        if (present(loc)) params(1) = loc
        if (present(scale)) params(2) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (params(2) <= 0.0_dp) return
        value = censored_loglikelihood(data, params, hypsecant_lpdf, hypsecant_lcdf, &
            hypsecant_lsf, hypsecant_interval)
    end function hypsecant_censored_loglikelihood

    pure function hypsecant_censored_nnlf(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite center, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -hypsecant_censored_loglikelihood(data, loc, scale)
    end function hypsecant_censored_nnlf

    pure function hypsecant_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! finite center and positive finite scale
        real(dp) :: value, z, t
        z = abs(location_scale_z(x,params))
        t = exp(-z)
        value = -z-log1p_safe(t*t)+(scifort_log_two-log(scifort_pi))-log(params(2))
    end function hypsecant_lpdf

    pure function hypsecant_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! finite center and positive finite scale
        real(dp) :: value, z
        z = location_scale_z(x,params)
        value = (scifort_log_two-log(scifort_pi))+cauchy_logatan_exp(-abs(z))
        if (z > 0.0_dp) value = log1p_safe(-exp(value))
    end function hypsecant_lcdf

    pure function hypsecant_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! finite center and positive finite scale
        real(dp) :: value, z
        z = location_scale_z(x,params)
        value = (scifort_log_two-log(scifort_pi))+cauchy_logatan_exp(-abs(z))
        if (z < 0.0_dp) value = log1p_safe(-exp(value))
    end function hypsecant_lsf

    pure function hypsecant_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper interval endpoint, greater than a
        real(dp), intent(in) :: params(:) !! finite center and positive finite scale
        real(dp) :: value, za, zb, logmass, logratio, omitted, large, small
        value = negative_infinity(0.0_dp)
        za = location_scale_z(a,params)
        zb = location_scale_z(b,params)
        logmass = log_exponential_mass(log_distance(b,a)-log(params(2)))
        if (a >= params(1) .or. b <= params(1)) then
            ! Reflect a negative interval; 0 <= za <= zb thereafter.
            if (b <= params(1)) then
                large = -za
                za = -zb
                zb = large
            end if
            if (.not. ieee_is_finite(za)) return
            ! atan(exp(-za))-atan(exp(-zb)) has a positive difference formula.
            logratio = -za+logmass-log1p_safe(exp(-za)*exp(-zb))
        else
            omitted = exp(hypsecant_lcdf(a,params))+exp(hypsecant_lsf(b,params))
            if (omitted <= 0.5_dp) then
                value = log1p_safe(-omitted)
                return
            end if
            ! tan(angle)=(1-exp(-(zb-za)))/(exp(za)+exp(-zb)).
            large = max(za,-zb)
            small = min(za,-zb)
            logratio = logmass-(large+log1p_safe(exp(small-large)))
        end if
        value = (scifort_log_two-log(scifort_pi))+cauchy_logatan_exp(logratio)
    end function hypsecant_interval

    pure function halfnorm_censored_loglikelihood(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! positive finite normal scale, default one
        real(dp) :: value, params(2)
        params = [0.0_dp, 1.0_dp]
        if (present(loc)) params(1) = loc
        if (present(scale)) params(2) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (params(2) <= 0.0_dp) return
        value = censored_loglikelihood(data, params, halfnorm_lpdf, halfnorm_lcdf, &
            halfnorm_lsf, halfnorm_interval)
    end function halfnorm_censored_loglikelihood

    pure function halfnorm_censored_nnlf(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! positive finite normal scale, default one
        real(dp) :: value
        value = -halfnorm_censored_loglikelihood(data, loc, scale)
    end function halfnorm_censored_nnlf

    pure function halfnorm_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! finite support endpoint and positive finite scale
        real(dp) :: value, z
        value = negative_infinity(0.0_dp)
        if (x < params(1)) return
        z = location_scale_z(x, params)
        if (z > sqrt(huge(1.0_dp))*sqrt(2.0_dp)) return
        value = (-0.5_dp*z)*z + (scifort_log_two-scifort_log_sqrt_two_pi) &
            - log(params(2))
    end function halfnorm_lpdf

    pure function halfnorm_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! finite support endpoint and positive finite scale
        real(dp) :: value, other
        call halfnorm_tails(x, params, value, other)
    end function halfnorm_lcdf

    pure function halfnorm_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! finite support endpoint and positive finite scale
        real(dp) :: value, other
        call halfnorm_tails(x, params, other, value)
    end function halfnorm_lsf

    pure subroutine halfnorm_tails(x, params, logcdf, logsf)
        implicit none
        real(dp), intent(in) :: x !! finite observation or censoring limit
        real(dp), intent(in) :: params(:) !! finite support endpoint and positive finite scale
        real(dp), intent(out) :: logcdf !! log lower-tail probability, -infinity at/below support
        real(dp), intent(out) :: logsf !! log upper-tail probability, zero at/below support
        real(dp) :: logz, z
        logcdf = negative_infinity(0.0_dp)
        logsf = 0.0_dp
        if (x <= params(1)) return
        logz = log_distance(x, params(1))-log(params(2))
        if (logz <= 0.0_dp) then
            if (logz <= log(0.01_dp)) then
                ! F(x)=2*P(loc < N(loc,scale) <= x); retain log widths below tiny.
                logcdf = scifort_log_two + normal_small_interval(0.0_dp, logz)
            else
                z = location_scale_z(x, params)
                logcdf = log(erf(z/sqrt(2.0_dp)))
            end if
            logsf = log1p_safe(-exp(logcdf))
        else
            z = location_scale_z(x, params)
            logsf = scifort_log_two + normal_logsf(z)
            logcdf = log1p_safe(-exp(logsf))
        end if
    end subroutine halfnorm_tails

    pure function halfnorm_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper interval endpoint, greater than a
        real(dp), intent(in) :: params(:) !! finite support endpoint and positive finite scale
        real(dp) :: value, lower_sf, upper_sf, z, logscaled
        value = negative_infinity(0.0_dp)
        if (b <= params(1)) return
        if (a <= params(1)) then
            value = halfnorm_lcdf(b, params)
            return
        end if
        lower_sf = halfnorm_lsf(a, params)
        if (.not. ieee_is_finite(lower_sf)) return
        upper_sf = halfnorm_lsf(b, params)
        if (upper_sf-lower_sf <= -scifort_log_two) then
            ! Preserve tiny excluded probabilities for nearly full-support intervals.
            value = lower_sf + log1p_safe(-exp(upper_sf-lower_sf))
        else
            z = location_scale_z(a, params)
            logscaled = log_distance(b,a)-log(params(2))
            value = normal_small_interval(z,logscaled)
            if (ieee_is_nan(value)) then
                value = tail_interval(a,b,params,halfnorm_lcdf,halfnorm_lsf)
            else
                value = value+scifort_log_two
            end if
        end if
    end function halfnorm_interval

    pure function gibrat_censored_loglikelihood(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = lognormal_censored_loglikelihood(data, 1.0_dp, loc, scale)
    end function gibrat_censored_loglikelihood

    pure function gibrat_censored_nnlf(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -gibrat_censored_loglikelihood(data, loc, scale)
    end function gibrat_censored_nnlf

    ! Erlang is gamma, including the existing positive noninteger-shape convention.
    pure function erlang_censored_loglikelihood(data, a, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: a !! finite positive shape, conventionally an integer
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = gamma_censored_loglikelihood(data, a, loc, scale)
    end function erlang_censored_loglikelihood

    pure function erlang_censored_nnlf(data, a, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: a !! finite positive shape, conventionally an integer
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -erlang_censored_loglikelihood(data, a, loc, scale)
    end function erlang_censored_nnlf

    pure function rayleigh_censored_loglikelihood(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive finite scale (default 1)
        real(dp) :: value, params(2)
        params = [0.0_dp, 1.0_dp]
        if (present(loc)) params(1) = loc
        if (present(scale)) params(2) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (params(2) <= 0.0_dp) return
        value = censored_loglikelihood(data, params, rayleigh_lpdf, rayleigh_lcdf, &
            rayleigh_lsf, rayleigh_interval)
    end function rayleigh_censored_loglikelihood

    pure function rayleigh_censored_nnlf(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive finite scale (default 1)
        real(dp) :: value
        value = -rayleigh_censored_loglikelihood(data, loc, scale)
    end function rayleigh_censored_nnlf

    pure function weibull_max_censored_loglikelihood(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! finite shape, > 0
        real(dp), intent(in), optional :: loc !! upper support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! scale, > 0 (default 1)
        real(dp) :: value, params(3)

        params = [c, 0.0_dp, 1.0_dp]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        value = quiet_nan(0.0_dp)
        if (.not. ieee_is_finite(c)) return
        if (c <= 0.0_dp .or. .not. valid_loc_scale(params(2), params(3))) return
        value = censored_loglikelihood(data, params, weibull_max_lpdf, weibull_max_lcdf, &
            weibull_max_lsf, weibull_max_interval)
    end function weibull_max_censored_loglikelihood

    pure function weibull_max_censored_nnlf(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! finite shape, > 0
        real(dp), intent(in), optional :: loc !! upper support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! scale, > 0 (default 1)
        real(dp) :: value
        value = -weibull_max_censored_loglikelihood(data, c, loc, scale)
    end function weibull_max_censored_nnlf

    pure function negative_hazard(logh) result(value)
        implicit none
        real(dp), intent(in) :: logh !! log nonnegative cumulative hazard, possibly infinite
        real(dp) :: value
        value = negative_infinity(0.0_dp)
        if (logh > log(huge(1.0_dp))) return
        value = -exp(logh)
    end function negative_hazard

    pure function log_hazard_mass(logh) result(value)
        implicit none
        real(dp), intent(in) :: logh !! log positive cumulative-hazard increment, possibly infinite
        real(dp) :: value
        if (logh <= log(scifort_log_two)) then
            value = log_exponential_mass(logh)
        else
            value = log1p_safe(-exp(negative_hazard(logh)))
        end if
    end function log_hazard_mass

    pure function rayleigh_log_hazard(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite point strictly above the support endpoint
        real(dp), intent(in) :: params(:) !! finite lower support endpoint and positive scale
        real(dp) :: value
        value = 2.0_dp*(log_distance(x,params(1))-log(params(2))) - scifort_log_two
    end function rayleigh_log_hazard

    pure function rayleigh_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! finite lower support endpoint and positive scale
        real(dp) :: value
        value = negative_infinity(0.0_dp)
        if (x <= params(1)) return
        value = negative_hazard(rayleigh_log_hazard(x,params))
        value = value + log_distance(x,params(1)) - 2.0_dp*log(params(2))
    end function rayleigh_lpdf

    pure function rayleigh_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! finite lower support endpoint and positive scale
        real(dp) :: value
        value = negative_infinity(0.0_dp)
        if (x <= params(1)) return
        value = log_hazard_mass(rayleigh_log_hazard(x,params))
    end function rayleigh_lcdf

    pure function rayleigh_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! finite lower support endpoint and positive scale
        real(dp) :: value
        value = 0.0_dp
        if (x <= params(1)) return
        value = negative_hazard(rayleigh_log_hazard(x,params))
    end function rayleigh_lsf

    pure function rayleigh_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper interval endpoint, > a
        real(dp), intent(in) :: params(:) !! finite lower support endpoint and positive scale
        real(dp) :: value, loga, logb, logdelta
        value = negative_infinity(0.0_dp)
        if (b <= params(1)) return
        if (a <= params(1)) then
            value = rayleigh_lcdf(b,params)
        else
            value = rayleigh_lsf(a,params)
            if (.not. ieee_is_finite(value)) return
            loga = log_distance(a,params(1))
            logb = log_distance(b,params(1))
            ! Delta H = (b-a)*((b-loc)+(a-loc))/(2*scale**2).
            logdelta = log_distance(b,a)+logb+log1p_safe(exp(loga-logb)) &
                -scifort_log_two-2.0_dp*log(params(2))
            value = value + log_hazard_mass(logdelta)
        end if
    end function rayleigh_interval

    pure function weibull_max_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! positive shape, finite upper support endpoint, scale
        real(dp) :: value, logh
        value = negative_infinity(0.0_dp)
        if (x > params(2)) return
        if (x == params(2)) then
            value = weibull_logpdf(0.0_dp,params(1),scale=params(3))
            return
        end if
        logh = weibull_log_hazard(-x,[params(1),-params(2),params(3)])
        if (logh > log(huge(1.0_dp))) return
        value = log(params(1))-log_distance(params(2),x)+logh-exp(logh)
    end function weibull_max_lpdf

    pure function weibull_max_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, finite upper support endpoint, scale
        real(dp) :: value
        value = 0.0_dp
        if (x >= params(2)) return
        value = negative_hazard(weibull_log_hazard(-x,[params(1),-params(2),params(3)]))
    end function weibull_max_lcdf

    pure function weibull_max_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, finite upper support endpoint, scale
        real(dp) :: value
        value = negative_infinity(0.0_dp)
        if (x >= params(2)) return
        value = log_hazard_mass(weibull_log_hazard(-x,[params(1),-params(2),params(3)]))
    end function weibull_max_lsf

    pure function weibull_max_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper interval endpoint, > a
        real(dp), intent(in) :: params(:) !! positive shape, finite upper support endpoint, scale
        real(dp) :: value, lower_cdf, upper_cdf
        value = negative_infinity(0.0_dp)
        if (a >= params(2)) return
        if (b >= params(2)) then
            value = weibull_max_lsf(a,params)
        else
            upper_cdf = weibull_max_lcdf(b,params)
            if (.not. ieee_is_finite(upper_cdf)) return
            lower_cdf = weibull_max_lcdf(a,params)
            if (lower_cdf-upper_cdf <= -scifort_log_two) then
                value = upper_cdf+log1p_safe(-exp(lower_cdf-upper_cdf))
            else
                value = weibull_interval(-b,-a,[params(1),-params(2),params(3)])
            end if
        end if
    end function weibull_max_interval


    pure function halfcauchy_censored_loglikelihood(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive finite scale (default 1)
        real(dp) :: value, params(2)
        params = [0.0_dp, 1.0_dp]
        if (present(loc)) params(1) = loc
        if (present(scale)) params(2) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (params(2) <= 0.0_dp) return
        value = censored_loglikelihood(data, params, halfcauchy_lpdf, halfcauchy_lcdf, &
            halfcauchy_lsf, halfcauchy_interval)
    end function halfcauchy_censored_loglikelihood

    pure function halfcauchy_censored_nnlf(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive finite scale (default 1)
        real(dp) :: value
        value = -halfcauchy_censored_loglikelihood(data, loc, scale)
    end function halfcauchy_censored_nnlf

    pure function halfcauchy_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! finite lower support endpoint and positive finite scale
        real(dp) :: value
        value = negative_infinity(0.0_dp)
        if (x < params(1)) return
        value = scifort_log_two + cauchy_lpdf(x, params)
    end function halfcauchy_lpdf

    pure function halfcauchy_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! finite lower support endpoint and positive finite scale
        real(dp) :: value, other
        call halfcauchy_tails(x, params, value, other)
    end function halfcauchy_lcdf

    pure function halfcauchy_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! finite lower support endpoint and positive finite scale
        real(dp) :: value, other
        call halfcauchy_tails(x, params, other, value)
    end function halfcauchy_lsf

    pure function halfcauchy_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper interval endpoint, > a
        real(dp), intent(in) :: params(:) !! finite lower support endpoint and positive finite scale
        real(dp) :: value, lower_sf, upper_sf
        value = negative_infinity(0.0_dp)
        if (b <= params(1)) return
        if (a <= params(1)) then
            value = halfcauchy_lcdf(b, params)
        else
            lower_sf = halfcauchy_lsf(a, params)
            if (.not. ieee_is_finite(lower_sf)) return
            upper_sf = halfcauchy_lsf(b, params)
            if (upper_sf - lower_sf <= -scifort_log_two) then
                ! Wide intervals retain tiny omitted probabilities near unit mass.
                value = lower_sf + log1p_safe(-exp(upper_sf - lower_sf))
            else
                value = scifort_log_two + cauchy_interval(a, b, params)
            end if
        end if
    end function halfcauchy_interval

    pure function halflogistic_censored_loglikelihood(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive finite scale (default 1)
        real(dp) :: value, params(2)
        params = [0.0_dp, 1.0_dp]
        if (present(loc)) params(1) = loc
        if (present(scale)) params(2) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (params(2) <= 0.0_dp) return
        value = censored_loglikelihood(data, params, halflogistic_lpdf, halflogistic_lcdf, &
            halflogistic_lsf, halflogistic_interval)
    end function halflogistic_censored_loglikelihood

    pure function halflogistic_censored_nnlf(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive finite scale (default 1)
        real(dp) :: value
        value = -halflogistic_censored_loglikelihood(data, loc, scale)
    end function halflogistic_censored_nnlf

    pure function halflogistic_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! finite lower support endpoint and positive finite scale
        real(dp) :: value
        value = negative_infinity(0.0_dp)
        if (x < params(1)) return
        if (x == params(1)) then
            value = -scifort_log_two - log(params(2))
        else
            value = scifort_log_two + logistic_lpdf(x, params)
        end if
    end function halflogistic_lpdf

    pure function halflogistic_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! finite lower support endpoint and positive finite scale
        real(dp) :: value, other
        call halflogistic_tails(x, params, value, other)
    end function halflogistic_lcdf

    pure function halflogistic_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! finite lower support endpoint and positive finite scale
        real(dp) :: value, other
        call halflogistic_tails(x, params, other, value)
    end function halflogistic_lsf

    pure function halflogistic_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper interval endpoint, > a
        real(dp), intent(in) :: params(:) !! finite lower support endpoint and positive finite scale
        real(dp) :: value, lower_sf, upper_sf
        value = negative_infinity(0.0_dp)
        if (b <= params(1)) return
        if (a <= params(1)) then
            value = halflogistic_lcdf(b, params)
        else
            lower_sf = halflogistic_lsf(a, params)
            if (.not. ieee_is_finite(lower_sf)) return
            upper_sf = halflogistic_lsf(b, params)
            if (upper_sf - lower_sf <= -scifort_log_two) then
                ! Wide intervals retain tiny omitted probabilities near unit mass.
                value = lower_sf + log1p_safe(-exp(upper_sf - lower_sf))
            else
                value = scifort_log_two + logistic_interval(a, b, params)
            end if
        end if
    end function halflogistic_interval

    pure subroutine halfcauchy_tails(x, params, logcdf, logsf)
        implicit none
        real(dp), intent(in) :: x !! finite observation or censoring limit
        real(dp), intent(in) :: params(:) !! finite lower support endpoint and positive finite scale
        real(dp), intent(out) :: logcdf !! log lower-tail probability, -infinity at/below support
        real(dp), intent(out) :: logsf !! log upper-tail probability, zero at/below support
        real(dp) :: logz
        logcdf = negative_infinity(0.0_dp)
        logsf = 0.0_dp
        if (x <= params(1)) return
        logz = cauchy_log_distance(x, params)
        if (logz <= 0.0_dp) then
            logcdf = scifort_log_two + cauchy_logatan_exp(logz) - log(scifort_pi)
            logsf = log1p_safe(-exp(logcdf))
        else
            logsf = scifort_log_two + cauchy_log_small_tail(logz)
            logcdf = log1p_safe(-exp(logsf))
        end if
    end subroutine halfcauchy_tails

    pure subroutine halflogistic_tails(x, params, logcdf, logsf)
        implicit none
        real(dp), intent(in) :: x !! finite observation or censoring limit
        real(dp), intent(in) :: params(:) !! finite lower support endpoint and positive finite scale
        real(dp), intent(out) :: logcdf !! log lower-tail probability, -infinity at/below support
        real(dp), intent(out) :: logsf !! log upper-tail probability, zero at/below support
        real(dp) :: logz
        logcdf = negative_infinity(0.0_dp)
        logsf = 0.0_dp
        if (x <= params(1)) return
        logz = cauchy_log_distance(x, params)
        if (logz <= 0.0_dp) then
            ! F(z) = (1-exp(-z))/(1+exp(-z)); retain log(z) if z underflows.
            logcdf = log_exponential_mass(logz) - log1pexp(-exp(logz))
            logsf = log1p_safe(-exp(logcdf))
        else
            logsf = scifort_log_two + logistic_lsf(x, params)
            logcdf = log1p_safe(-exp(logsf))
        end if
    end subroutine halflogistic_tails


    pure function cauchy_censored_loglikelihood(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite median (default 0)
        real(dp), intent(in), optional :: scale !! positive finite half width at half maximum (default 1)
        real(dp) :: value, params(2)
        params = [0.0_dp, 1.0_dp]
        if (present(loc)) params(1) = loc
        if (present(scale)) params(2) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (params(2) <= 0.0_dp) return
        value = censored_loglikelihood(data, params, cauchy_lpdf, cauchy_lcdf, &
            cauchy_lsf, cauchy_interval)
    end function cauchy_censored_loglikelihood

    pure function cauchy_censored_nnlf(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite median (default 0)
        real(dp), intent(in), optional :: scale !! positive finite half width at half maximum (default 1)
        real(dp) :: value
        value = -cauchy_censored_loglikelihood(data, loc, scale)
    end function cauchy_censored_nnlf

    pure function cauchy_log_distance(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite observation or censoring limit
        real(dp), intent(in) :: params(:) !! finite median and positive finite scale
        real(dp) :: value
        value = negative_infinity(0.0_dp)
        if (x == params(1)) return
        value = log_distance(max(x, params(1)), min(x, params(1))) - log(params(2))
    end function cauchy_log_distance

    pure function cauchy_logatan_exp(logx) result(value)
        implicit none
        real(dp), intent(in) :: logx !! log positive argument; negative infinity means zero
        real(dp) :: value
        if (logx < 0.5_dp * log(epsilon(1.0_dp))) then
            ! atan(t)/t differs from 1 by at most t**2/3 for 0 <= t <= 1.
            value = logx
        else if (logx <= 0.0_dp) then
            value = log(atan(exp(logx)))
        else
            ! Reciprocal identity keeps EXP bounded and subtraction well conditioned.
            value = log(0.5_dp * scifort_pi - atan(exp(-logx)))
        end if
    end function cauchy_logatan_exp

    pure function cauchy_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! finite median and positive finite scale
        real(dp) :: value, logz
        logz = cauchy_log_distance(x, params)
        value = -log(scifort_pi) - log(params(2)) - log1pexp(2.0_dp * logz)
    end function cauchy_lpdf

    pure function cauchy_log_small_tail(logz) result(value)
        implicit none
        real(dp), intent(in) :: logz !! log absolute standardized displacement, possibly -infinity
        real(dp) :: value
        if (logz >= 0.0_dp) then
            value = cauchy_logatan_exp(-logz) - log(scifort_pi)
        else
            value = log(0.5_dp - atan(exp(logz)) / scifort_pi)
        end if
    end function cauchy_log_small_tail

    pure function cauchy_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! finite median and positive finite scale
        real(dp) :: value
        value = cauchy_log_small_tail(cauchy_log_distance(x, params))
        if (x > params(1)) value = log1p_safe(-exp(value))
    end function cauchy_lcdf

    pure function cauchy_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! finite median and positive finite scale
        real(dp) :: value
        value = cauchy_log_small_tail(cauchy_log_distance(x, params))
        if (x < params(1)) value = log1p_safe(-exp(value))
    end function cauchy_lsf

    pure function cauchy_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper interval endpoint, > a
        real(dp), intent(in) :: params(:) !! finite median and positive finite scale
        real(dp) :: value, loga, logb, logratio, left, right, large, small
        loga = cauchy_log_distance(a, params)
        logb = cauchy_log_distance(b, params)
        if (b <= params(1) .or. a >= params(1)) then
            ! Same-side arctangent difference has a positive denominator:
            ! tan(angle) = ((b-a)/scale) / (1 + za*zb).
            logratio = log_distance(b, a) - log(params(2)) - log1pexp(loga + logb)
            value = cauchy_logatan_exp(logratio) - log(scifort_pi)
        else if (min(loga, logb) >= 0.0_dp) then
            ! Wide crossing interval: retain even tiny missing tail probabilities.
            left = exp(cauchy_log_small_tail(loga))
            right = exp(cauchy_log_small_tail(logb))
            value = log1p_safe(-(left + right))
        else
            ! Narrow crossing interval: add positive angles, never subtract CDFs.
            left = cauchy_logatan_exp(loga)
            right = cauchy_logatan_exp(logb)
            large = max(left, right)
            small = min(left, right)
            value = min(0.0_dp, large + log1p_safe(exp(small - large)) - log(scifort_pi))
        end if
    end function cauchy_interval


    pure function laplace_censored_loglikelihood(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! center (default 0)
        real(dp), intent(in), optional :: scale !! positive finite mean absolute deviation (default 1)
        real(dp) :: value, params(2)

        params = [0.0_dp, 1.0_dp]
        if (present(loc)) params(1) = loc
        if (present(scale)) params(2) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (params(2) <= 0.0_dp) return
        value = censored_loglikelihood(data, params, laplace_lpdf, laplace_lcdf, &
            laplace_lsf, laplace_interval)
    end function laplace_censored_loglikelihood

    pure function laplace_censored_nnlf(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! center (default 0)
        real(dp), intent(in), optional :: scale !! positive finite mean absolute deviation (default 1)
        real(dp) :: value
        value = -laplace_censored_loglikelihood(data, loc, scale)
    end function laplace_censored_nnlf

    pure function laplace_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! finite center and positive finite scale
        real(dp) :: value
        value = -abs(location_scale_z(x, params)) - scifort_log_two - log(params(2))
    end function laplace_lpdf

    pure function laplace_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! upper censoring limit
        real(dp), intent(in) :: params(:) !! finite center and positive finite scale
        real(dp) :: value, z
        z = location_scale_z(x, params)
        if (z <= 0.0_dp) then
            value = z - scifort_log_two
        else
            value = log1p_safe(-0.5_dp * exp(-z))
        end if
    end function laplace_lcdf

    pure function laplace_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! lower censoring limit
        real(dp), intent(in) :: params(:) !! finite center and positive finite scale
        real(dp) :: value, z
        z = location_scale_z(x, params)
        if (z >= 0.0_dp) then
            value = -z - scifort_log_two
        else
            value = log1p_safe(-0.5_dp * exp(z))
        end if
    end function laplace_lsf

    pure function laplace_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper interval endpoint, > a
        real(dp), intent(in) :: params(:) !! finite center and positive finite scale
        real(dp) :: value, left, right, large, small, logscale
        logscale = log(params(2))
        if (b <= params(1)) then
            ! Left-side exponential mass, anchored at the endpoint nearest the center.
            value = laplace_lcdf(b, params) + &
                log_exponential_mass(log_distance(b, a) - logscale)
        else if (a >= params(1)) then
            value = laplace_lsf(a, params) + &
                log_exponential_mass(log_distance(b, a) - logscale)
        else
            ! Split at the cusp; log-add the two positive half-masses.
            left = log_exponential_mass(log_distance(params(1), a) - logscale)
            right = log_exponential_mass(log_distance(b, params(1)) - logscale)
            large = max(left, right)
            small = min(left, right)
            value = large + log1p_safe(0.5_dp * expm1_safe(small - large))
        end if
    end function laplace_interval


    pure function gumbel_r_censored_loglikelihood(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! location (default 0)
        real(dp), intent(in), optional :: scale !! positive finite scale (default 1)
        real(dp) :: value, params(3)

        params = [0.0_dp, 1.0_dp, 1.0_dp]
        if (present(loc)) params(1) = loc
        if (present(scale)) params(2) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (params(2) <= 0.0_dp) return
        value = censored_loglikelihood(data, params, gumbel_lpdf, gumbel_lcdf, &
            gumbel_lsf, gumbel_interval)
    end function gumbel_r_censored_loglikelihood

    pure function gumbel_r_censored_nnlf(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! location (default 0)
        real(dp), intent(in), optional :: scale !! positive finite scale (default 1)
        real(dp) :: value
        value = -gumbel_r_censored_loglikelihood(data, loc, scale)
    end function gumbel_r_censored_nnlf

    pure function gumbel_l_censored_loglikelihood(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! location (default 0)
        real(dp), intent(in), optional :: scale !! positive finite scale (default 1)
        real(dp) :: value, params(3)

        params = [0.0_dp, 1.0_dp, -1.0_dp]
        if (present(loc)) params(1) = loc
        if (present(scale)) params(2) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (params(2) <= 0.0_dp) return
        value = censored_loglikelihood(data, params, gumbel_lpdf, gumbel_lcdf, &
            gumbel_lsf, gumbel_interval)
    end function gumbel_l_censored_loglikelihood

    pure function gumbel_l_censored_nnlf(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! location (default 0)
        real(dp), intent(in), optional :: scale !! positive finite scale (default 1)
        real(dp) :: value
        value = -gumbel_l_censored_loglikelihood(data, loc, scale)
    end function gumbel_l_censored_nnlf

    pure function gumbel_right_logcdf(z) result(value)
        implicit none
        real(dp), intent(in) :: z !! standardized point for the right-skewed family
        real(dp) :: value
        value = negative_infinity(0.0_dp)
        if (-z > log(huge(1.0_dp))) return
        value = -exp(-z)
    end function gumbel_right_logcdf

    pure function gumbel_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! exact observation
        real(dp), intent(in) :: params(:) !! location, positive scale, direction (+1 or -1)
        real(dp) :: value, z
        z = params(3) * location_scale_z(x, params)
        value = gumbel_right_logcdf(z)
        if (.not. ieee_is_finite(value)) return
        value = value - z - log(params(2))
    end function gumbel_lpdf

    pure function gumbel_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! upper censoring limit
        real(dp), intent(in) :: params(:) !! location, positive scale, direction (+1 or -1)
        real(dp) :: value, z
        z = params(3) * location_scale_z(x, params)
        if (params(3) > 0.0_dp) then
            value = gumbel_right_logcdf(z)
        else
            value = log_exponential_mass(-z)
        end if
    end function gumbel_lcdf

    pure function gumbel_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! lower censoring limit
        real(dp), intent(in) :: params(:) !! location, positive scale, direction (+1 or -1)
        real(dp) :: value, z
        z = params(3) * location_scale_z(x, params)
        if (params(3) > 0.0_dp) then
            value = log_exponential_mass(-z)
        else
            value = gumbel_right_logcdf(z)
        end if
    end function gumbel_lsf

    pure function gumbel_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower endpoint
        real(dp), intent(in) :: b !! finite upper endpoint, > a
        real(dp), intent(in) :: params(:) !! location, positive scale, direction (+1 or -1)
        real(dp) :: value, za, zb, logwidth, logincrement
        if (params(3) > 0.0_dp) then
            za = location_scale_z(a, params)
            zb = location_scale_z(b, params)
        else
            ! Reflection interchanges the interval endpoints and the two tails.
            za = -location_scale_z(b, params)
            zb = -location_scale_z(a, params)
        end if
        logwidth = log_distance(b, a) - log(params(2))
        ! exp(-za)-exp(-zb) = exp(-za)*(1-exp(-(b-a)/scale)).
        logincrement = -za + log_exponential_mass(logwidth)
        value = gumbel_right_logcdf(zb) + log_exponential_mass(logincrement)
    end function gumbel_interval

    pure function logistic_censored_loglikelihood(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! center (default 0)
        real(dp), intent(in), optional :: scale !! positive finite scale (default 1)
        real(dp) :: value, params(2)

        params = [0.0_dp, 1.0_dp]
        if (present(loc)) params(1) = loc
        if (present(scale)) params(2) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (params(2) <= 0.0_dp) return
        value = censored_loglikelihood(data, params, logistic_lpdf, logistic_lcdf, &
            logistic_lsf, logistic_interval)
    end function logistic_censored_loglikelihood

    pure function logistic_censored_nnlf(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! center (default 0)
        real(dp), intent(in), optional :: scale !! positive finite scale (default 1)
        real(dp) :: value
        value = -logistic_censored_loglikelihood(data, loc, scale)
    end function logistic_censored_nnlf

    pure function location_scale_z(x, params) result(z)
        implicit none
        real(dp), intent(in) :: x !! evaluation point, possibly infinite
        real(dp), intent(in) :: params(:) !! finite center and positive finite scale
        real(dp) :: z, delta, factor, ratio
        integer :: power

        if (.not. ieee_is_finite(x)) then
            z = x
            return
        end if
        factor = 1.0_dp
        if ((x > 0.0_dp .and. params(1) < 0.0_dp) .or. &
            (x < 0.0_dp .and. params(1) > 0.0_dp)) then
            if (abs(x) >= huge(x) - abs(params(1))) factor = 2.0_dp
        end if
        ! Halve near subtraction overflow; retain nearby differences elsewhere.
        delta = x / factor - params(1) / factor
        if (delta == 0.0_dp) then
            z = delta
            return
        end if
        z = sign(positive_infinity(0.0_dp), delta)
        ! Divide bounded significands, then scale by a checked binary exponent.
        ! A rounded huge*scale threshold can otherwise miss division overflow.
        ratio = fraction(delta) / fraction(params(2))
        power = exponent(delta) - exponent(params(2)) + exponent(ratio)
        if (factor == 2.0_dp) power = power + 1
        if (power > maxexponent(z)) return
        z = scale(fraction(ratio), power)
    end function location_scale_z

    pure function logistic_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! exact observation
        real(dp), intent(in) :: params(:) !! finite center and positive finite scale
        real(dp) :: value, a
        a = abs(location_scale_z(x, params))
        value = -a - 2.0_dp * log1p_safe(exp(-a)) - log(params(2))
    end function logistic_lpdf

    pure function logistic_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! upper censoring limit
        real(dp), intent(in) :: params(:) !! finite center and positive finite scale
        real(dp) :: value
        value = -log1pexp(-location_scale_z(x, params))
    end function logistic_lcdf

    pure function logistic_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! lower censoring limit
        real(dp), intent(in) :: params(:) !! finite center and positive finite scale
        real(dp) :: value
        value = -log1pexp(location_scale_z(x, params))
    end function logistic_lsf

    pure function logistic_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower endpoint
        real(dp), intent(in) :: b !! finite upper endpoint, > a
        real(dp), intent(in) :: params(:) !! finite center and positive finite scale
        real(dp) :: value
        ! F(b)-F(a) = F(b)*S(a)*(1-exp(-(b-a)/scale)).
        value = logistic_lcdf(b, params) + logistic_lsf(a, params) + &
            log_exponential_mass(log_distance(b, a) - log(params(2)))
    end function logistic_interval

    function make_censored_data(uncensored, left, right, interval, status) result(data)
        implicit none
        real(dp), intent(in), optional :: uncensored(:) !! finite exact observations
        real(dp), intent(in), optional :: left(:) !! finite upper censoring limits
        real(dp), intent(in), optional :: right(:) !! finite lower censoring limits
        real(dp), intent(in), optional :: interval(:,:) !! (n,2) lower/upper limits
        integer, intent(out), optional :: status !! 0 valid, 1 invalid; invalid data is retained
        type(censored_data) :: data

        if (present(uncensored)) data%uncensored = uncensored
        if (present(left)) data%left = left
        if (present(right)) data%right = right
        if (present(interval)) data%interval = interval
        if (present(status)) then
            status = censored_status_invalid_input
            if (censored_data_valid(data)) status = censored_status_success
        end if
    end function make_censored_data

    pure function censored_data_valid(data) result(valid)
        implicit none
        type(censored_data), intent(in) :: data !! observations; unallocated groups are empty
        logical :: valid
        real(dp) :: a, b
        integer :: i

        valid = .false.
        if (allocated(data%uncensored)) then
            if (.not. all(ieee_is_finite(data%uncensored))) return
        end if
        if (allocated(data%left)) then
            if (.not. all(ieee_is_finite(data%left))) return
        end if
        if (allocated(data%right)) then
            if (.not. all(ieee_is_finite(data%right))) return
        end if
        if (allocated(data%interval)) then
            if (size(data%interval, 2) /= 2) return
            do i = 1, size(data%interval, 1)
                a = data%interval(i,1)
                b = data%interval(i,2)
                if (ieee_is_nan(a) .or. ieee_is_nan(b)) return
                if (a > b) return
                if (.not. ieee_is_finite(a) .and. .not. ieee_is_finite(b)) return
            end do
        end if
        valid = .true.
    end function censored_data_valid

    pure function censored_size(data) result(n)
        implicit none
        type(censored_data), intent(in) :: data !! observations; unallocated groups are empty
        integer :: n
        n = 0
        if (allocated(data%uncensored)) n = n + size(data%uncensored)
        if (allocated(data%left)) n = n + size(data%left)
        if (allocated(data%right)) n = n + size(data%right)
        if (allocated(data%interval)) n = n + size(data%interval, 1)
    end function censored_size

    pure function censored_loglikelihood(data, params, logpdf, logcdf, logsf, &
            loginterval) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: params(:) !! finite parameters; callbacks validate their domain
        procedure(censored_log_function) :: logpdf !! log density callback
        procedure(censored_log_function) :: logcdf !! stable log lower-tail callback
        procedure(censored_log_function) :: logsf !! stable log upper-tail callback
        procedure(censored_log_interval), optional :: loginterval !! optional stable interval mass
        real(dp) :: value
        real(dp) :: a, b, term
        integer :: i

        value = quiet_nan(0.0_dp)
        if (.not. censored_data_valid(data)) return
        if (.not. all(ieee_is_finite(params))) return
        value = 0.0_dp
        if (allocated(data%uncensored)) then
            do i = 1, size(data%uncensored)
                value = add_log_term(value, logpdf(data%uncensored(i), params))
            end do
        end if
        if (allocated(data%left)) then
            do i = 1, size(data%left)
                value = add_log_term(value, logcdf(data%left(i), params))
            end do
        end if
        if (allocated(data%right)) then
            do i = 1, size(data%right)
                value = add_log_term(value, logsf(data%right(i), params))
            end do
        end if
        if (allocated(data%interval)) then
            do i = 1, size(data%interval, 1)
                a = data%interval(i,1)
                b = data%interval(i,2)
                if (a == b) then
                    term = logpdf(a, params)
                else if (.not. ieee_is_finite(a)) then
                    term = logcdf(b, params)
                else if (.not. ieee_is_finite(b)) then
                    term = logsf(a, params)
                else if (present(loginterval)) then
                    term = loginterval(a, b, params)
                else
                    term = tail_interval(a, b, params, logcdf, logsf)
                end if
                value = add_log_term(value, term)
            end do
        end if
    end function censored_loglikelihood

    pure function add_log_term(total, term) result(value)
        implicit none
        real(dp), intent(in) :: total !! accumulated log likelihood, possibly infinite
        real(dp), intent(in) :: term !! next contribution, possibly infinite
        real(dp) :: value
        value = quiet_nan(0.0_dp)
        if (ieee_is_nan(total) .or. ieee_is_nan(term)) return
        if (.not. ieee_is_finite(total) .and. .not. ieee_is_finite(term)) then
            if (total /= term) return
        end if
        value = total + term
    end function add_log_term

    pure function log_difference(large, small) result(value)
        implicit none
        real(dp), intent(in) :: large !! log of larger nonnegative probability
        real(dp), intent(in) :: small !! log of smaller nonnegative probability
        real(dp) :: value

        value = quiet_nan(large)
        if (ieee_is_nan(large) .or. ieee_is_nan(small)) return
        if (small > large) return
        if (large == small) then
            value = negative_infinity(large)
        else
            value = large + log(-expm1_safe(small - large))
        end if
    end function log_difference

    pure function tail_interval(a, b, params, logcdf, logsf) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower endpoint
        real(dp), intent(in) :: b !! finite upper endpoint, > a
        real(dp), intent(in) :: params(:) !! distribution parameters
        procedure(censored_log_function) :: logcdf !! stable log lower-tail callback
        procedure(censored_log_function) :: logsf !! stable log upper-tail callback
        real(dp) :: value, upper_cdf, lower_sf

        upper_cdf = logcdf(b, params)
        lower_sf = logsf(a, params)
        if (upper_cdf <= lower_sf) then
            value = log_difference(upper_cdf, logcdf(a, params))
        else
            value = log_difference(lower_sf, logsf(b, params))
        end if
    end function tail_interval

    pure function log_distance(b, a) result(value)
        implicit none
        real(dp), intent(in) :: b !! finite upper endpoint, > a
        real(dp), intent(in) :: a !! finite lower endpoint
        real(dp) :: value, big, small

        if (a < 0.0_dp .and. b > 0.0_dp) then
            big = max(abs(a), abs(b))
            small = min(abs(a), abs(b))
            value = log(big) + log1p_safe(small / big)
        else
            value = log(b - a)
        end if
    end function log_distance

    pure function log_exponential_mass(logwidth) result(value)
        implicit none
        real(dp), intent(in) :: logwidth !! log positive cumulative-hazard increment
        real(dp) :: value

        if (logwidth < log(epsilon(1.0_dp)) - 2.0_dp) then
            value = logwidth
        else if (logwidth > log(huge(1.0_dp))) then
            value = 0.0_dp
        else
            value = log(-expm1_safe(-exp(logwidth)))
        end if
    end function log_exponential_mass

    pure function exponential_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower endpoint
        real(dp), intent(in) :: b !! finite upper endpoint, > a
        real(dp), intent(in) :: params(:) !! location and positive scale
        real(dp) :: value, lower

        value = negative_infinity(0.0_dp)
        if (b <= params(1)) return
        lower = max(a, params(1))
        value = exponential_logsf(lower, params(1), params(2)) + &
            log_exponential_mass(log_distance(b, lower) - log(params(2)))
    end function exponential_interval

    pure function normal_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower endpoint
        real(dp), intent(in) :: b !! finite upper endpoint, > a
        real(dp), intent(in) :: params(:) !! mean and positive standard deviation
        real(dp) :: value, z, logwidth, logscaled

        z = (a - params(1)) / params(2)
        logwidth = log_distance(b, a)
        logscaled = logwidth - log(params(2))
        value = normal_small_interval(z, logscaled)
        if (ieee_is_nan(value)) value = tail_interval(a, b, params, normal_lcdf, normal_lsf)
    end function normal_interval

    pure function normal_small_interval(lower_z, logscaled) result(value)
        implicit none
        real(dp), intent(in) :: lower_z !! standardized lower endpoint
        real(dp), intent(in) :: logscaled !! log standardized interval width
        real(dp) :: value, z, halfwidth, offset, correction

        ! Three-point Gauss-Legendre density integral; NaN requests tail fallback.
        ! Nodes 0,+/-sqrt(3/5), weights 8/9,5/9,5/9 (DLMF 3.5(v)).
        value = quiet_nan(0.0_dp)
        z = lower_z
        if (ieee_is_finite(z)) then
            if (logscaled + log(max(1.0_dp, abs(z))) <= log(0.01_dp)) then
                if (abs(z) > sqrt(huge(1.0_dp))*sqrt(2.0_dp)) then
                    value = negative_infinity(0.0_dp)
                    return
                end if
                halfwidth = 0.5_dp*exp(logscaled)
                z = z + halfwidth
                offset = halfwidth*sqrt(3.0_dp/5.0_dp)
                correction = (5.0_dp/18.0_dp)*( &
                    expm1_safe(-z*offset - 0.5_dp*offset**2) + &
                    expm1_safe(z*offset - 0.5_dp*offset**2))
                value = (-0.5_dp*z)*z - scifort_log_sqrt_two_pi + &
                    logscaled + log1p_safe(correction)
                return
            end if
        end if
    end function normal_small_interval

    pure function lognormal_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower endpoint
        real(dp), intent(in) :: b !! finite upper endpoint, > a
        real(dp), intent(in) :: params(:) !! positive shape, location, positive scale
        real(dp) :: value, logr, logspan, z

        value = negative_infinity(0.0_dp)
        if (b <= params(2)) return
        if (a <= params(2)) then
            value = lognormal_lcdf(b, params)
            return
        end if
        ! log(b-loc)-log(a-loc) loses adjacent endpoint widths. Instead use
        ! log1p((b-a)/(a-loc)), retaining its logarithm if the ratio underflows.
        logr = log_distance(b, a) - log_distance(a, params(2))
        if (logr < log(epsilon(1.0_dp)) - 2.0_dp) then
            logspan = logr
        else
            logspan = log(log1pexp(logr))
        end if
        z = lognormal_standardized(a, params)
        value = normal_small_interval(z, logspan - log(params(1)))
        if (ieee_is_nan(value)) then
            value = tail_interval(a, b, params, lognormal_lcdf, lognormal_lsf)
        end if
    end function lognormal_interval

    pure function lognormal_standardized(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite point above location
        real(dp), intent(in) :: params(:) !! positive shape, location, positive scale
        real(dp) :: value, logz
        logz = lognormal_logz(x, params)
        if (params(1) < 1.0_dp) then
            if (abs(logz) > huge(1.0_dp)*params(1)) then
                value = sign(positive_infinity(0.0_dp), logz)
                return
            end if
        end if
        value = logz/params(1)
    end function lognormal_standardized

    pure function lognormal_logz(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite point above location
        real(dp), intent(in) :: params(:) !! positive shape, location, positive scale
        real(dp) :: value
        ! Reuse the stable log distance/scale calculation, with unit exponent.
        value = weibull_log_hazard(x, [1.0_dp, params(2), params(3)])
    end function lognormal_logz

    pure function weibull_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower endpoint
        real(dp), intent(in) :: b !! finite upper endpoint, > a
        real(dp), intent(in) :: params(:) !! positive shape, location, positive scale
        real(dp) :: value, lower, loglo, loghi, logratio, logr, delta, logdelta

        value = negative_infinity(0.0_dp)
        if (b <= params(2)) return
        lower = max(a, params(2))
        loghi = weibull_log_hazard(b, params)
        if (lower == params(2)) then
            value = log_exponential_mass(loghi)
            return
        end if
        loglo = weibull_log_hazard(lower, params)
        if (loglo > log(huge(1.0_dp))) return
        logratio = log_distance(b, lower) - log_distance(lower, params(2))
        if (logratio < log(epsilon(1.0_dp)) - 2.0_dp) then
            logr = log(params(1)) + logratio
        else
            logr = log(params(1)) + log(log1pexp(logratio))
        end if
        ! Here the omitted factor exp(-delta) is below exp(-50) < 2e-22.
        if (logr > log(50.0_dp)) then
            logdelta = loghi
        else if (logr < log(epsilon(1.0_dp)) - 2.0_dp) then
            logdelta = loglo + logr
        else
            delta = exp(logr)
            logdelta = loglo + log(expm1_safe(delta))
        end if
        value = -exp(loglo) + log_exponential_mass(logdelta)
    end function weibull_interval

    pure function weibull_log_hazard(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite point strictly above location
        real(dp), intent(in) :: params(:) !! positive shape, location, positive scale
        real(dp) :: value, logz, distance
        logical :: finite_distance
        logz = log_distance(x, params(2)) - log(params(3))
        ! Avoid subtracting nearly equal large logarithms when the normalized
        ! distance is near one: a large shape would amplify that lost precision.
        finite_distance = .true.
        if (params(2) < 0.0_dp .and. x > 0.0_dp) then
            finite_distance = x <= huge(1.0_dp) + params(2)
        end if
        if (finite_distance) then
            distance = x - params(2)
            if (0.5_dp*distance <= params(3) .and. distance >= 0.5_dp*params(3)) then
                logz = log1p_safe((distance - params(3))/params(3))
            end if
        end if
        if (params(1) > 1.0_dp) then
            if (abs(logz) > huge(1.0_dp)/params(1)) then
                value = sign(positive_infinity(0.0_dp), logz)
                return
            end if
        end if
        value = params(1)*logz
    end function weibull_log_hazard

    pure function normal_censored_loglikelihood(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! mean (default 0)
        real(dp), intent(in), optional :: scale !! standard deviation, > 0 (default 1)
        real(dp) :: value, params(2)

        params = [0.0_dp, 1.0_dp]
        if (present(loc)) params(1) = loc
        if (present(scale)) params(2) = scale
        value = quiet_nan(0.0_dp)
        if (.not. valid_loc_scale(params(1), params(2))) return
        value = censored_loglikelihood(data, params, normal_lpdf, normal_lcdf, &
            normal_lsf, normal_interval)
    end function normal_censored_loglikelihood

    pure function exponential_censored_loglikelihood(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! inverse rate, > 0 (default 1)
        real(dp) :: value, params(2)

        params = [0.0_dp, 1.0_dp]
        if (present(loc)) params(1) = loc
        if (present(scale)) params(2) = scale
        value = quiet_nan(0.0_dp)
        if (.not. valid_loc_scale(params(1), params(2))) return
        value = censored_loglikelihood(data, params, exponential_lpdf, exponential_lcdf, &
            exponential_lsf, exponential_interval)
    end function exponential_censored_loglikelihood

    pure function weibull_censored_loglikelihood(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! finite shape, > 0
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! scale, > 0 (default 1)
        real(dp) :: value, params(3)

        params = [c, 0.0_dp, 1.0_dp]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        value = quiet_nan(0.0_dp)
        if (.not. ieee_is_finite(c)) return
        if (c <= 0.0_dp .or. .not. valid_loc_scale(params(2), params(3))) return
        value = censored_loglikelihood(data, params, weibull_lpdf, weibull_lcdf, &
            weibull_lsf, weibull_interval)
    end function weibull_censored_loglikelihood

    pure function normal_censored_nnlf(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! mean (default 0)
        real(dp), intent(in), optional :: scale !! standard deviation, > 0 (default 1)
        real(dp) :: value
        value = -normal_censored_loglikelihood(data, loc, scale)
    end function normal_censored_nnlf

    pure function exponential_censored_nnlf(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! inverse rate, > 0 (default 1)
        real(dp) :: value
        value = -exponential_censored_loglikelihood(data, loc, scale)
    end function exponential_censored_nnlf

    pure function weibull_censored_nnlf(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! finite shape, > 0
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! scale, > 0 (default 1)
        real(dp) :: value
        value = -weibull_censored_loglikelihood(data, c, loc, scale)
    end function weibull_censored_nnlf

    pure function normal_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! evaluation point
        real(dp), intent(in) :: params(:) !! mean and standard deviation
        real(dp) :: value
        value = normal_logpdf(x, params(1), params(2))
    end function normal_lpdf

    pure function normal_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! evaluation point
        real(dp), intent(in) :: params(:) !! mean and standard deviation
        real(dp) :: value
        value = normal_logcdf(x, params(1), params(2))
    end function normal_lcdf

    pure function normal_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! evaluation point
        real(dp), intent(in) :: params(:) !! mean and standard deviation
        real(dp) :: value
        value = normal_logsf(x, params(1), params(2))
    end function normal_lsf

    pure function exponential_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! evaluation point
        real(dp), intent(in) :: params(:) !! location and scale
        real(dp) :: value
        value = exponential_logpdf(x, params(1), params(2))
    end function exponential_lpdf

    pure function exponential_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! evaluation point
        real(dp), intent(in) :: params(:) !! location and scale
        real(dp) :: value
        value = exponential_logcdf(x, params(1), params(2))
    end function exponential_lcdf

    pure function exponential_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! evaluation point
        real(dp), intent(in) :: params(:) !! location and scale
        real(dp) :: value
        value = exponential_logsf(x, params(1), params(2))
    end function exponential_lsf

    pure function weibull_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! evaluation point
        real(dp), intent(in) :: params(:) !! shape, location, scale
        real(dp) :: value
        value = weibull_logpdf(x, params(1), params(2), params(3))
    end function weibull_lpdf

    pure function weibull_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! evaluation point
        real(dp), intent(in) :: params(:) !! shape, location, scale
        real(dp) :: value
        value = weibull_logcdf(x, params(1), params(2), params(3))
    end function weibull_lcdf

    pure function weibull_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! evaluation point
        real(dp), intent(in) :: params(:) !! shape, location, scale
        real(dp) :: value
        value = weibull_logsf(x, params(1), params(2), params(3))
    end function weibull_lsf

    pure function lognormal_censored_loglikelihood(data, s, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: s !! finite log standard deviation, > 0
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! finite scale, > 0 (default 1)
        real(dp) :: value, params(3)

        params = [s, 0.0_dp, 1.0_dp]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (s <= 0.0_dp .or. .not. valid_loc_scale(params(2), params(3))) return
        value = censored_loglikelihood(data, params, lognormal_lpdf, lognormal_lcdf, &
            lognormal_lsf, lognormal_interval)
    end function lognormal_censored_loglikelihood

    pure function lognormal_censored_nnlf(data, s, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: s !! finite log standard deviation, > 0
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! finite scale, > 0 (default 1)
        real(dp) :: value
        value = -lognormal_censored_loglikelihood(data, s, loc, scale)
    end function lognormal_censored_nnlf

    pure function lognormal_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! positive shape, location, positive scale
        real(dp) :: value, z
        value = negative_infinity(0.0_dp)
        if (x <= params(2)) return
        z = lognormal_standardized(x, params)
        if (abs(z) > sqrt(huge(1.0_dp))*sqrt(2.0_dp)) return
        value = normal_logpdf(z) - log(params(1)) - log_distance(x, params(2))
    end function lognormal_lpdf

    pure function lognormal_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, location, positive scale
        real(dp) :: value
        value = negative_infinity(0.0_dp)
        if (x <= params(2)) return
        value = normal_logcdf(lognormal_standardized(x, params))
    end function lognormal_lcdf

    pure function lognormal_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, location, positive scale
        real(dp) :: value
        value = 0.0_dp
        if (x <= params(2)) return
        value = normal_logsf(lognormal_standardized(x, params))
    end function lognormal_lsf

    pure function gamma_censored_loglikelihood(data, a, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: a !! finite shape, > 0
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! finite inverse rate, > 0 (default 1)
        real(dp) :: value, params(3)
        params = [a, 0.0_dp, 1.0_dp]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (a <= 0.0_dp .or. params(3) <= 0.0_dp) return
        if (a == 1.0_dp) then
            value = exponential_censored_loglikelihood(data, params(2), params(3))
        else
            value = censored_loglikelihood(data, params, gamma_lpdf, gamma_lcdf, &
                gamma_lsf, gamma_interval)
        end if
    end function gamma_censored_loglikelihood

    pure function gamma_censored_nnlf(data, a, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: a !! finite shape, > 0
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! finite inverse rate, > 0 (default 1)
        real(dp) :: value
        value = -gamma_censored_loglikelihood(data, a, loc, scale)
    end function gamma_censored_nnlf

    pure function gamma_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! positive shape, location, positive scale
        real(dp) :: value, logz, distance_log, z
        value = negative_infinity(0.0_dp)
        if (x < params(2)) return
        if (x == params(2)) then
            value = gamma_logpdf(x, params(1), params(2), params(3))
            return
        end if
        distance_log = log_distance(x, params(2))
        logz = distance_log - log(params(3))
        if (logz > log(huge(1.0_dp))) return
        if (logz < log(tiny(1.0_dp))) then
            ! The omitted exp(-z) correction is below the smallest normal real.
            value = log(params(1))+gamma_small_kernel(params(1),logz)-distance_log
        else
            z = gamma_scaled(x,params,logz)
            if (params(1) < 10.0_dp) then
                value = (params(1)-1.0_dp)*logz-z-log_gamma(params(1))-log(params(3))
            else
                value = log(params(1))+log_gamma_kernel(params(1),z)-distance_log
            end if
        end if
    end function gamma_lpdf

    pure function gamma_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, location, positive scale
        real(dp) :: value, logz
        value = negative_infinity(0.0_dp)
        if (x <= params(2)) return
        logz = log_distance(x, params(2))-log(params(3))
        if (logz > log(huge(1.0_dp))) then
            value = 0.0_dp
        else if (logz < log(tiny(1.0_dp))) then
            ! P(a,z) = z**a/Gamma(a+1) * (1 + O(z)); DLMF 8.7.
            value = gamma_small_kernel(params(1),logz)
        else
            value = log_gammainc(params(1), gamma_scaled(x,params,logz))
        end if
    end function gamma_lcdf

    pure function gamma_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, location, positive scale
        real(dp) :: value, logz, logp
        value = 0.0_dp
        if (x <= params(2)) return
        logz = log_distance(x, params(2))-log(params(3))
        if (logz > log(huge(1.0_dp))) then
            value = negative_infinity(0.0_dp)
        else if (logz < log(tiny(1.0_dp))) then
            logp = gamma_lcdf(x, params)
            value = log_exponential_mass(log(-logp))
        else
            value = log_gammaincc(params(1), gamma_scaled(x,params,logz))
        end if
    end function gamma_lsf

    pure function gamma_small_kernel(a, logz) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite positive shape
        real(dp), intent(in) :: logz !! log standardized argument below log(tiny)
        real(dp) :: value, factor
        if (a < 10.0_dp) then
            value = a*logz-log_gamma_one_plus(a)
        else
            ! Stirling form avoids overflow in Gamma(a+1) and in a*log(z).
            factor = -logz+log(a)-1.0_dp
            value = negative_infinity(0.0_dp)
            if (a > huge(1.0_dp)/factor) return
            value = -a*factor-0.5_dp*log(a)-scifort_log_sqrt_two_pi-stirling_remainder(a)
        end if
    end function gamma_small_kernel

    pure function gamma_scaled(x, params, logz) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite point strictly above location
        real(dp), intent(in) :: params(:) !! positive shape, location, positive scale
        real(dp), intent(in) :: logz !! log standardized distance within normal real range
        real(dp) :: value
        ! Direct division retains near-mode accuracy for large shapes. Use the
        ! logarithm only if the unscaled distance or ratio could overflow.
        value = exp(logz)
        if (logz > log(huge(1.0_dp))-1.0_dp) return
        if (params(2) < 0.0_dp .and. x > 0.0_dp) then
            if (x > huge(1.0_dp)+params(2)) return
        end if
        value = (x-params(2))/params(3)
    end function gamma_scaled

    pure function gamma_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper endpoint, > a
        real(dp), intent(in) :: params(:) !! positive shape, location, positive scale
        real(dp) :: value, logwidth, logrelative, logscaled, q, d, correction, t(3), r(3)
        value = negative_infinity(0.0_dp)
        if (b <= params(2)) return
        if (a <= params(2)) then
            value = gamma_lcdf(b, params)
            return
        end if
        logwidth = log_distance(b,a)
        logrelative = logwidth-log_distance(a,params(2))
        logscaled = logwidth-log(params(3))
        ! Bound both relative endpoint separation and density variation.
        ! Three-point Gauss-Legendre is exact through degree five (DLMF 3.5).
        if (logrelative+log(max(1.0_dp,abs(params(1)-1.0_dp))) <= log(0.001_dp) &
                .and. logscaled <= log(0.001_dp)) then
            q = exp(logrelative)
            d = exp(logscaled)
            t = [0.5_dp, 0.5_dp-sqrt(3.0_dp/5.0_dp)/2.0_dp, &
                0.5_dp+sqrt(3.0_dp/5.0_dp)/2.0_dp]
            r = (params(1)-1.0_dp)*log1p_safe(t*q)-t*d
            correction = (4.0_dp/9.0_dp)*expm1_safe(r(1)) + &
                (5.0_dp/18.0_dp)*(expm1_safe(r(2))+expm1_safe(r(3)))
            value = gamma_lpdf(a,params)+logwidth+log1p_safe(correction)
        else
            value = tail_interval(a,b,params,gamma_lcdf,gamma_lsf)
        end if
    end function gamma_interval

    pure function beta_censored_loglikelihood(data, a, b, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: a !! finite first shape, > 0
        real(dp), intent(in) :: b !! finite second shape, > 0
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! finite support width, > 0 (default 1)
        real(dp) :: value, params(4)
        params = [a, b, 0.0_dp, 1.0_dp]
        if (present(loc)) params(3) = loc
        if (present(scale)) params(4) = scale
        value = quiet_nan(0.0_dp)
        if (.not. beta_parameters_valid(params)) return
        value = censored_loglikelihood(data,params,beta_lpdf,beta_lcdf,beta_lsf,beta_interval)
    end function beta_censored_loglikelihood

    pure function beta_censored_nnlf(data, a, b, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: a !! finite first shape, > 0
        real(dp), intent(in) :: b !! finite second shape, > 0
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! finite support width, > 0 (default 1)
        real(dp) :: value
        value = -beta_censored_loglikelihood(data,a,b,loc,scale)
    end function beta_censored_nnlf

    pure function beta_parameters_valid(params) result(valid)
        implicit none
        real(dp), intent(in) :: params(:) !! first/second shapes, location, support width
        logical :: valid
        valid = .false.
        if (.not. all(ieee_is_finite(params))) return
        if (min(params(1),params(2),params(4)) <= 0.0_dp) return
        ! This initial interface requires a finite rounded upper support edge.
        if (params(3) > 0.0_dp) then
            if (params(4) > huge(1.0_dp)-params(3)) return
        end if
        valid = .true.
    end function beta_parameters_valid

    pure subroutine beta_distances(x, params, left, right)
        implicit none
        real(dp), intent(in) :: x !! finite observation or censoring limit
        real(dp), intent(in) :: params(:) !! valid first/second shapes, location, width
        real(dp), intent(out) :: left !! distance above lower edge; negative outside support
        real(dp), intent(out) :: right !! distance below upper edge; negative outside support
        real(dp) :: upper, correction
        ! FastTwoSum with operands ordered by magnitude (Ogita et al. 2005).
        upper = params(3)+params(4)
        if (abs(params(3)) >= params(4)) then
            correction = (params(3)-upper)+params(4)
        else
            correction = (params(4)-upper)+params(3)
        end if
        left = 0.0_dp
        right = params(4)
        if (x == params(3)) return
        if (x < params(3)) then
            left = -1.0_dp
            return
        end if
        left = positive_infinity(0.0_dp)
        if (params(3) < 0.0_dp .and. x > 0.0_dp) then
            if (x <= huge(1.0_dp)+params(3)) left = x-params(3)
        else
            left = x-params(3)
        end if
        if (x > 0.0_dp .and. upper < 0.0_dp) then
            if (x > huge(1.0_dp)+upper) then
                right = -1.0_dp
                return
            end if
        else if (x < 0.0_dp .and. upper > 0.0_dp) then
            if (upper > huge(1.0_dp)+x) then
                right = huge(1.0_dp)
                return
            end if
        end if
        right = (upper-x)+correction
    end subroutine beta_distances

    pure subroutine beta_coordinates(x, params, left, right, z, w, logz, logw)
        implicit none
        real(dp), intent(in) :: x !! finite point strictly inside support
        real(dp), intent(in) :: params(:) !! valid first/second shapes, location, width
        real(dp), intent(in) :: left !! positive distance above lower edge, possibly infinite
        real(dp), intent(in) :: right !! positive distance below upper edge
        real(dp), intent(out) :: z !! standardized lower distance, possibly rounded to 0 or 1
        real(dp), intent(out) :: w !! standardized upper distance, formed independently of z
        real(dp), intent(out) :: logz !! finite log standardized lower distance
        real(dp), intent(out) :: logw !! finite log standardized upper distance
        if (ieee_is_finite(left)) then
            z = left/params(4)
        else
            z = x/params(4)-params(3)/params(4)
        end if
        w = right/params(4)
        logz = log_distance(x,params(3))-log(params(4))
        logw = log(right)-log(params(4))
        if (w < 0.5_dp) logz = log1p_safe(-w)
        if (z < 0.5_dp) logw = log1p_safe(-z)
    end subroutine beta_coordinates

    pure function beta_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! valid first/second shapes, location, width
        real(dp) :: value, left, right, z, w, logz, logw, shape, other
        value = negative_infinity(0.0_dp)
        call beta_distances(x,params,left,right)
        if (left < 0.0_dp .or. right < 0.0_dp) return
        if (left == 0.0_dp .or. right == 0.0_dp) then
            shape = params(1)
            other = params(2)
            if (right == 0.0_dp) then
                shape = params(2)
                other = params(1)
            end if
            if (shape < 1.0_dp) value = positive_infinity(0.0_dp)
            if (shape == 1.0_dp) value = log(other)-log(params(4))
            return
        end if
        call beta_coordinates(x,params,left,right,z,w,logz,logw)
        if (min(params(1),params(2)) >= stirling_threshold .and. min(z,w) > 0.0_dp) then
            value = log_beta_kernel(params(1),params(2),z,w)-logz-logw-log(params(4))
        else
            value = (params(1)-1.0_dp)*logz+(params(2)-1.0_dp)*logw - &
                log_beta(params(1),params(2))-log(params(4))
        end if
    end function beta_lpdf

    pure subroutine beta_logtails(x, params, logp, logq)
        implicit none
        real(dp), intent(in) :: x !! finite censoring limit
        real(dp), intent(in) :: params(:) !! valid first/second shapes, location, width
        real(dp), intent(out) :: logp !! stable log lower probability
        real(dp), intent(out) :: logq !! stable log upper probability
        real(dp) :: left, right, z, w, logz, logw, p, q
        logp = negative_infinity(0.0_dp)
        logq = 0.0_dp
        call beta_distances(x,params,left,right)
        if (left <= 0.0_dp) return
        if (right <= 0.0_dp) then
            logp = 0.0_dp
            logq = negative_infinity(0.0_dp)
            return
        end if
        call beta_coordinates(x,params,left,right,z,w,logz,logw)
        if (z == 0.0_dp) then
            logp = quiet_nan(0.0_dp)
            logq = logp
            if (log(params(2))+logz > log(epsilon(1.0_dp))-2.0_dp) return
            logp = beta_small_tail(params(1),params(2),logz)
            logq = log_difference(0.0_dp,logp)
        else if (w == 0.0_dp) then
            logp = quiet_nan(0.0_dp)
            logq = logp
            if (log(params(1))+logw > log(epsilon(1.0_dp))-2.0_dp) return
            logq = beta_small_tail(params(2),params(1),logw)
            logp = log_difference(0.0_dp,logq)
        else
            call incomplete_beta_xy(params(1),params(2),z,w,p,q,logp,logq)
        end if
    end subroutine beta_logtails

    pure function beta_small_tail(a, b, logz) result(value)
        implicit none
        real(dp), intent(in) :: a !! positive shape at the nearby endpoint
        real(dp), intent(in) :: b !! positive shape at the far endpoint
        real(dp), intent(in) :: logz !! log coordinate rounded to zero in ordinary arithmetic
        real(dp) :: value, normalization
        ! The gamma-difference form retains the O(a) normalization when a is tiny.
        if (a < sqrt(epsilon(1.0_dp))*min(b,1.0_dp)) then
            if (b > 1.0_dp/epsilon(1.0_dp)) then
                ! For a < sqrt(epsilon), the gamma-ratio correction is O(a/b),
                ! below roundoff in this tail; avoid powers of enormous b.
                normalization = a*log(b)-log_gamma_one_plus(a)
            else if (b >= stirling_threshold) then
                normalization = a*log(b)+log_gamma_ratio_scaled(b,a)-log_gamma_one_plus(a)
            else
                normalization = log_gamma_difference(b,a)-log_gamma_one_plus(a)
            end if
        else
            normalization = -log(a)-log_beta(a,b)
        end if
        value = a*logz+normalization
    end function beta_small_tail

    pure function beta_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! valid first/second shapes, location, width
        real(dp) :: value, logq
        call beta_logtails(x,params,value,logq)
    end function beta_lcdf

    pure function beta_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! valid first/second shapes, location, width
        real(dp) :: value, logp
        call beta_logtails(x,params,logp,value)
    end function beta_lsf

    pure function beta_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper interval endpoint, > a
        real(dp), intent(in) :: params(:) !! valid first/second shapes, location, width
        real(dp) :: value, left, right, logwidth, logql, logqr, ql, qr, t(3), r(3), correction
        value = negative_infinity(0.0_dp)
        call beta_distances(b,params,left,right)
        if (left <= 0.0_dp) return
        call beta_distances(a,params,left,right)
        if (right <= 0.0_dp) return
        if (left <= 0.0_dp) then
            value = beta_lcdf(b,params)
            return
        end if
        call beta_distances(b,params,left,right)
        if (right <= 0.0_dp) then
            value = beta_lsf(a,params)
            return
        end if
        call beta_distances(a,params,left,right)
        logwidth = log_distance(b,a)
        logql = logwidth-log_distance(a,params(3))
        logqr = logwidth-log(right)
        if (logql+log(max(1.0_dp,abs(params(1)-1.0_dp))) <= log(0.001_dp) .and. &
                logqr+log(max(1.0_dp,abs(params(2)-1.0_dp))) <= log(0.001_dp)) then
            ql = exp(logql)
            qr = exp(logqr)
            t = [0.5_dp, 0.5_dp-sqrt(3.0_dp/5.0_dp)/2.0_dp, &
                0.5_dp+sqrt(3.0_dp/5.0_dp)/2.0_dp]
            r = (params(1)-1.0_dp)*log1p_safe(t*ql) + &
                (params(2)-1.0_dp)*log1p_safe(-t*qr)
            correction = (4.0_dp/9.0_dp)*expm1_safe(r(1)) + &
                (5.0_dp/18.0_dp)*(expm1_safe(r(2))+expm1_safe(r(3)))
            value = beta_lpdf(a,params)+logwidth+log1p_safe(correction)
        else
            value = tail_interval(a,b,params,beta_lcdf,beta_lsf)
        end if
    end function beta_interval

end module scifort_censored
