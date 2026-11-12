! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
! Input-driven comparison driver, intentionally outside unattended FPM tests.
program censored_reference_driver
    use scifort_censored, only : censored_data, normal_censored_loglikelihood, &
        exponential_censored_loglikelihood, weibull_censored_loglikelihood, &
        lognormal_censored_loglikelihood, gamma_censored_loglikelihood, &
        beta_censored_loglikelihood
    use scifort_kinds, only : dp
    implicit none
    type(censored_data) :: data
    integer :: family, nu, nl, nr, ni, ios, i
    real(dp) :: c, second_shape, loc, scale, value
    do
        read(*, *, iostat=ios) family, c, loc, scale, nu, nl, nr, ni
        if (ios < 0) exit
        if (ios /= 0) error stop 'invalid comparison input'
        if (family == 6) read(*, *) second_shape
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
        case default
            error stop 'unknown comparison family'
        end select
        write(*, '(es26.17e3)') value
        deallocate(data%uncensored, data%left, data%right, data%interval)
    end do
end program censored_reference_driver
