! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
! Input-driven comparison driver, intentionally outside unattended FPM tests.
program censored_reference_driver
    use scifort_censored, only : censored_data, normal_censored_loglikelihood, &
        exponential_censored_loglikelihood, weibull_censored_loglikelihood, &
        lognormal_censored_loglikelihood, gamma_censored_loglikelihood, &
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
    use scifort_kinds, only : dp
    implicit none
    type(censored_data) :: data
    integer :: family, nu, nl, nr, ni, ios, i
    real(dp) :: c, second_shape, loc, scale, value
    do
        read(*, *, iostat=ios) family, c, loc, scale, nu, nl, nr, ni
        if (ios < 0) exit
        if (ios /= 0) error stop 'invalid comparison input'
        if (family == 6 .or. family == 25 .or. family == 26) read(*, *) second_shape
        allocate(data%uncensored(nu), data%left(nl), data%right(nr), data%interval(ni,2))
        if (nu > 0) read(*, *) data%uncensored
        if (nl > 0) read(*, *) data%left
        if (nr > 0) read(*, *) data%right
        do i = 1, ni
            read(*, *) data%interval(i,:)
        end do
        select case (family)
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
        case default
            error stop 'unknown comparison family'
        end select
        write(*, '(es26.17e3)') value
        deallocate(data%uncensored, data%left, data%right, data%interval)
    end do
end program censored_reference_driver
