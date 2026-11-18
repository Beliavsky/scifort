! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
! Input-driven comparison driver, intentionally outside unattended FPM tests.
program censored_reference_driver
    use scifort_censored, only : censored_data, normal_censored_loglikelihood, &
        exponential_censored_loglikelihood, weibull_censored_loglikelihood, &
        lognormal_censored_loglikelihood, gamma_censored_loglikelihood, &
        genlogistic_censored_loglikelihood, &
        gompertz_censored_loglikelihood, &
        invweibull_censored_loglikelihood, &
        moyal_censored_loglikelihood, &
        levy_censored_loglikelihood, &
        levy_l_censored_loglikelihood, &
        dgamma_censored_loglikelihood, dweibull_censored_loglikelihood, &
        loggamma_censored_loglikelihood, &
        invgamma_censored_loglikelihood, betaprime_censored_loglikelihood, &
        maxwell_censored_loglikelihood, &
        bradford_censored_loglikelihood, &
        triang_censored_loglikelihood, &
        trapezoid_censored_loglikelihood, &
        erlang_censored_loglikelihood, &
        gibrat_censored_loglikelihood, &
        halfnorm_censored_loglikelihood, &
        uniform_censored_loglikelihood, powerlaw_censored_loglikelihood, &
        arcsine_censored_loglikelihood, &
        lomax_censored_loglikelihood, fisk_censored_loglikelihood, &
        burr_censored_loglikelihood, burr12_censored_loglikelihood, &
        hypsecant_censored_loglikelihood, pareto_censored_loglikelihood, &
        laplace_asymmetric_censored_loglikelihood, loglaplace_censored_loglikelihood, &
        beta_censored_loglikelihood, logistic_censored_loglikelihood, &
        laplace_censored_loglikelihood, cauchy_censored_loglikelihood, &
        halfcauchy_censored_loglikelihood, halflogistic_censored_loglikelihood, &
        rayleigh_censored_loglikelihood, weibull_max_censored_loglikelihood, &
        gumbel_r_censored_loglikelihood, gumbel_l_censored_loglikelihood
    use scifort_censored, only : powernorm_censored_loglikelihood, &
        powerlognorm_censored_loglikelihood, skewcauchy_censored_loglikelihood
    use scifort_censored, only : exponpow_censored_loglikelihood, &
        exponweib_censored_loglikelihood, mielke_censored_loglikelihood
    use scifort_censored, only : fatiguelife_censored_loglikelihood, &
        johnsonsu_censored_loglikelihood, genexpon_censored_loglikelihood
    use scifort_censored, only : cosine_censored_loglikelihood, &
        anglit_censored_loglikelihood, semicircular_censored_loglikelihood
    use scifort_censored, only : truncexpon_censored_loglikelihood, &
        johnsonsb_censored_loglikelihood, kappa3_censored_loglikelihood
    use scifort_censored, only : loguniform_censored_loglikelihood, &
        truncpareto_censored_loglikelihood, genhalflogistic_censored_loglikelihood
    use scifort_censored, only : foldnorm_censored_loglikelihood, &
        foldcauchy_censored_loglikelihood, genpareto_censored_loglikelihood
    use scifort_censored, only : gennorm_censored_loglikelihood, &
        halfgennorm_censored_loglikelihood, rdist_censored_loglikelihood
    use scifort_censored, only : chi2_censored_loglikelihood
    use scifort_censored, only : chi_censored_loglikelihood
    use scifort_censored, only : nakagami_censored_loglikelihood
    use scifort_censored, only : t_censored_loglikelihood
    use scifort_censored, only : f_censored_loglikelihood
    use scifort_censored, only : jf_skew_t_censored_loglikelihood
    use scifort_censored, only : invgauss_censored_loglikelihood
    use scifort_censored, only : wald_censored_loglikelihood
    use scifort_censored, only : recipinvgauss_censored_loglikelihood
    use scifort_censored, only : gengamma_censored_loglikelihood
    use scifort_censored, only : truncnorm_censored_loglikelihood
    use scifort_censored, only : truncweibull_min_censored_loglikelihood
    use scifort_censored, only : alpha_censored_loglikelihood
    use scifort_censored, only : genextreme_censored_loglikelihood
    use scifort_censored, only : crystalball_censored_loglikelihood
    use scifort_censored, only : argus_censored_loglikelihood
    use scifort_censored, only : pearson3_censored_loglikelihood
    use scifort_censored, only : kappa4_censored_loglikelihood
    use scifort_censored, only : gausshyper_censored_loglikelihood
    use scifort_censored, only : irwinhall_censored_loglikelihood
    use scifort_censored, only : vonmises_line_censored_loglikelihood
    use scifort_censored, only : exponnorm_censored_loglikelihood
    use scifort_censored, only : skewnorm_censored_loglikelihood
    use scifort_censored, only : dpareto_lognorm_censored_loglikelihood
    use scifort_censored, only : ncx2_censored_loglikelihood
    use scifort_censored, only : rice_censored_loglikelihood
    use scifort_censored, only : ncf_censored_loglikelihood
    use scifort_censored, only : kstwobign_censored_loglikelihood
    use scifort_censored, only : wrapcauchy_censored_loglikelihood
    use scifort_censored, only : rel_breitwigner_censored_loglikelihood
    use scifort_kinds, only : dp
    implicit none
    type(censored_data) :: data
    integer :: family, nu, nl, nr, ni, ios, i
    real(dp) :: c, second_shape, third_shape, fourth_shape, loc, scale, value
    do
        read(*, *, iostat=ios) family, c, loc, scale, nu, nl, nr, ni
        if (ios < 0) exit
        if (ios /= 0) error stop 'invalid comparison input'
        if (any(family == [6,25,26,38,43,46,49,50,52,53,58,60,61,73,74,78,79,80,83,86,87,92,93,95])) &
            read(*, *) second_shape
        if (family == 53 .or. family == 80 .or. family == 87 .or. family == 92 .or. family == 95) read(*, *) third_shape
        if (family == 87 .or. family == 92) read(*, *) fourth_shape
        allocate(data%uncensored(nu), data%left(nl), data%right(nr), data%interval(ni,2))
        if (nu > 0) read(*, *) data%uncensored
        if (nl > 0) read(*, *) data%left
        if (nr > 0) read(*, *) data%right
        do i = 1, ni
            read(*, *) data%interval(i,:)
        end do
        select case (family)
        case (81)
            value = alpha_censored_loglikelihood(data,c,loc,scale)
        case (82)
            value = genextreme_censored_loglikelihood(data,c,loc,scale)
        case (83)
            value = crystalball_censored_loglikelihood(data,c,second_shape,loc,scale)
        case (96)
            value = kstwobign_censored_loglikelihood(data,loc,scale)
        case (97)
            value = wrapcauchy_censored_loglikelihood(data,c,loc,scale)
        case (98)
            value = rel_breitwigner_censored_loglikelihood(data,c,loc,scale)
        case (93)
            value = ncx2_censored_loglikelihood(data,c,second_shape,loc,scale)
        case (94)
            value = rice_censored_loglikelihood(data,c,loc,scale)
        case (95)
            value = ncf_censored_loglikelihood(data,c,second_shape,third_shape,loc,scale)
        case (90)
            value = exponnorm_censored_loglikelihood(data,c,loc,scale)
        case (91)
            value = skewnorm_censored_loglikelihood(data,c,loc,scale)
        case (92)
            value = dpareto_lognorm_censored_loglikelihood(data,c,second_shape, &
                third_shape,fourth_shape,loc,scale)
        case (87)
            value = gausshyper_censored_loglikelihood(data,c,second_shape,third_shape,fourth_shape,loc,scale)
        case (88)
            value = irwinhall_censored_loglikelihood(data,c,loc,scale)
        case (89)
            value = vonmises_line_censored_loglikelihood(data,c,loc,scale)
        case (84)
            value = argus_censored_loglikelihood(data,c,loc,scale)
        case (85)
            value = pearson3_censored_loglikelihood(data,c,loc,scale)
        case (86)
            value = kappa4_censored_loglikelihood(data,c,second_shape,loc,scale)
        case (78)
            value = gengamma_censored_loglikelihood(data,c,second_shape,loc,scale)
        case (79)
            value = truncnorm_censored_loglikelihood(data,c,second_shape,loc,scale)
        case (80)
            value = truncweibull_min_censored_loglikelihood(data,c,second_shape,third_shape,loc,scale)
        case (75)
            value = invgauss_censored_loglikelihood(data,c,loc,scale)
        case (76)
            value = wald_censored_loglikelihood(data,loc,scale)
        case (77)
            value = recipinvgauss_censored_loglikelihood(data,c,loc,scale)
        case (72)
            value = t_censored_loglikelihood(data,c,loc,scale)
        case (73)
            value = f_censored_loglikelihood(data,c,second_shape,loc,scale)
        case (74)
            value = jf_skew_t_censored_loglikelihood(data,c,second_shape,loc,scale)
        case (69)
            value = chi2_censored_loglikelihood(data,c,loc,scale)
        case (70)
            value = chi_censored_loglikelihood(data,c,loc,scale)
        case (71)
            value = nakagami_censored_loglikelihood(data,c,loc,scale)
        case (66)
            value = gennorm_censored_loglikelihood(data,c,loc,scale)
        case (67)
            value = halfgennorm_censored_loglikelihood(data,c,loc,scale)
        case (68)
            value = rdist_censored_loglikelihood(data,c,loc,scale)
        case (63)
            value = foldnorm_censored_loglikelihood(data,c,loc,scale)
        case (64)
            value = foldcauchy_censored_loglikelihood(data,c,loc,scale)
        case (65)
            value = genpareto_censored_loglikelihood(data,c,loc,scale)
        case (60)
            value = loguniform_censored_loglikelihood(data,c,second_shape,loc,scale)
        case (61)
            value = truncpareto_censored_loglikelihood(data,c,second_shape,loc,scale)
        case (62)
            value = genhalflogistic_censored_loglikelihood(data,c,loc,scale)
        case (57)
            value = truncexpon_censored_loglikelihood(data,c,loc,scale)
        case (58)
            value = johnsonsb_censored_loglikelihood(data,c,second_shape,loc,scale)
        case (59)
            value = kappa3_censored_loglikelihood(data,c,loc,scale)
        case (54)
            value = cosine_censored_loglikelihood(data,loc,scale)
        case (55)
            value = anglit_censored_loglikelihood(data,loc,scale)
        case (56)
            value = semicircular_censored_loglikelihood(data,loc,scale)
        case (45)
            value = powernorm_censored_loglikelihood(data,c,loc,scale)
        case (46)
            value = powerlognorm_censored_loglikelihood(data,c,second_shape,loc,scale)
        case (47)
            value = skewcauchy_censored_loglikelihood(data,c,loc,scale)
        case (48)
            value = exponpow_censored_loglikelihood(data,c,loc,scale)
        case (49)
            value = exponweib_censored_loglikelihood(data,c,second_shape,loc,scale)
        case (50)
            value = mielke_censored_loglikelihood(data,c,second_shape,loc,scale)
        case (51)
            value = fatiguelife_censored_loglikelihood(data,c,loc,scale)
        case (52)
            value = johnsonsu_censored_loglikelihood(data,c,second_shape,loc,scale)
        case (53)
            value = genexpon_censored_loglikelihood(data,c,second_shape,third_shape,loc,scale)
        case (1)
            value = normal_censored_loglikelihood(data, loc, scale)
        case (2)
            value = exponential_censored_loglikelihood(data, loc, scale)
        case (3)
            value = weibull_censored_loglikelihood(data, c, loc, scale)
        case (4)
            value = lognormal_censored_loglikelihood(data, c, loc, scale)
        case (5)
            value = gamma_censored_loglikelihood(data, c, loc, scale)
        case (6)
            value = beta_censored_loglikelihood(data,c,second_shape,loc,scale)
        case (7)
            value = logistic_censored_loglikelihood(data,loc,scale)
        case (8)
            value = gumbel_r_censored_loglikelihood(data,loc,scale)
        case (9)
            value = gumbel_l_censored_loglikelihood(data,loc,scale)
        case (10)
            value = laplace_censored_loglikelihood(data,loc,scale)
        case (11)
            value = cauchy_censored_loglikelihood(data,loc,scale)
        case (12)
            value = halfcauchy_censored_loglikelihood(data,loc,scale)
        case (13)
            value = halflogistic_censored_loglikelihood(data,loc,scale)
        case (14)
            value = rayleigh_censored_loglikelihood(data,loc,scale)
        case (15)
            value = weibull_max_censored_loglikelihood(data,c,loc,scale)
        case (16)
            value = erlang_censored_loglikelihood(data,c,loc,scale)
        case (17)
            value = gibrat_censored_loglikelihood(data,loc,scale)
        case (18)
            value = halfnorm_censored_loglikelihood(data,loc,scale)
        case (19)
            value = hypsecant_censored_loglikelihood(data,loc,scale)
        case (20)
            value = laplace_asymmetric_censored_loglikelihood(data,c,loc,scale)
        case (21)
            value = loglaplace_censored_loglikelihood(data,c,loc,scale)
        case (22)
            value = pareto_censored_loglikelihood(data,c,loc,scale)
        case (23)
            value = lomax_censored_loglikelihood(data,c,loc,scale)
        case (24)
            value = fisk_censored_loglikelihood(data,c,loc,scale)
        case (25)
            value = burr_censored_loglikelihood(data,c,second_shape,loc,scale)
        case (26)
            value = burr12_censored_loglikelihood(data,c,second_shape,loc,scale)
        case (27)
            value = uniform_censored_loglikelihood(data,loc,scale)
        case (28)
            value = powerlaw_censored_loglikelihood(data,c,loc,scale)
        case (29)
            value = arcsine_censored_loglikelihood(data,loc,scale)
        case (30)
            value = genlogistic_censored_loglikelihood(data,c,loc,scale)
        case (31)
            value = gompertz_censored_loglikelihood(data,c,loc,scale)
        case (32)
            value = invweibull_censored_loglikelihood(data,c,loc,scale)
        case (33)
            value = moyal_censored_loglikelihood(data,loc,scale)
        case (34)
            value = levy_censored_loglikelihood(data,loc,scale)
        case (35)
            value = levy_l_censored_loglikelihood(data,loc,scale)
        case (36)
            value = bradford_censored_loglikelihood(data,c,loc,scale)
        case (37)
            value = triang_censored_loglikelihood(data,c,loc,scale)
        case (38)
            value = trapezoid_censored_loglikelihood(data,c,second_shape,loc,scale)
        case (39)
            value = dgamma_censored_loglikelihood(data,c,loc,scale)
        case (40)
            value = dweibull_censored_loglikelihood(data,c,loc,scale)
        case (41)
            value = loggamma_censored_loglikelihood(data,c,loc,scale)
        case (42)
            value = invgamma_censored_loglikelihood(data,c,loc,scale)
        case (43)
            value = betaprime_censored_loglikelihood(data,c,second_shape,loc,scale)
        case (44)
            value = maxwell_censored_loglikelihood(data,loc,scale)
        case default
            error stop 'unknown comparison family'
        end select
        write(*, '(es26.17e3)') value
        deallocate(data%uncensored, data%left, data%right, data%interval)
    end do
end program censored_reference_driver
