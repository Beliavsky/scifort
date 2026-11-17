! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_censored_radial_gamma
    use scifort_stats, only : dp, censored_data, make_censored_data, fit_result, &
        fit_status_invalid_input, gamma_censored_loglikelihood, &
        halfnorm_censored_loglikelihood, rayleigh_censored_loglikelihood, &
        maxwell_censored_loglikelihood, exponential_censored_loglikelihood, &
        chi2_censored_loglikelihood, chi2_censored_nnlf, chi2_fit, &
        chi2_logpdf, chi2_logcdf, chi2_logsf, chi2_ppf, chi2_isf, &
        chi_censored_loglikelihood, chi_censored_nnlf, chi_fit, &
        chi_logpdf, chi_logcdf, chi_logsf, chi_ppf, chi_isf, &
        nakagami_censored_loglikelihood, nakagami_censored_nnlf, nakagami_fit, &
        nakagami_logpdf, nakagami_logcdf, nakagami_logsf, nakagami_ppf, nakagami_isf
    use scifort_math, only : positive_infinity, quiet_nan
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    implicit none
    type(censored_data) :: data, one, empty, mixed
    type(fit_result) :: fit
    real(dp) :: x, shape, actual, expected, p, inf, nan, previous, scale, loc
    real(dp) :: lower(3), upper(3), anchor, tiny_shape
    real(dp), parameter :: scipy_shape(3) = &
        [0.9099424378451767_dp,0.5690179321507661_dp,0.3_dp]
    real(dp), parameter :: scipy_objective(3) = &
        [11.137830726064342_dp,10.012294312027775_dp,10.649162735292245_dp]
    real(dp), parameter :: minimum_shape_mixed(3) = &
        [-1489.1038694251238_dp,-1488.016999932435_dp,-1481.2175065125543_dp]
    real(dp), parameter :: maximum_shape_center(3) = &
        [-9.906135507935513_dp,-1.9586594581134578_dp,5.295669411148652_dp]
    real(dp), parameter :: maximum_shape_adjacent(3) = &
        [-30.700550853998337_dp,-29.684546609775737_dp,-29.36168954611308_dp]
    integer :: family, i, j, k, checks, failures
    checks = 0
    failures = 0
    inf = positive_infinity(0.0_dp)
    nan = quiet_nan(0.0_dp)
    empty = make_censored_data()
    mixed = make_censored_data(uncensored=[0.2_dp,0.4_dp,0.6_dp],left=[0.1_dp], &
        right=[0.7_dp],interval=reshape([0.25_dp,0.5_dp,0.35_dp,0.55_dp],[2,2]))
    do family = 1, 3
        do j = 1, 5
            shape = real(j,dp)/3.0_dp
            previous = -inf
            do i = 1, 20
                x = real(i,dp)/5.0_dp
                one = make_censored_data(uncensored=[x])
                call close(ll(one,family,shape),scalar(x,family,shape,1),2e-12_dp,'density')
                one = make_censored_data(left=[x])
                actual = ll(one,family,shape)
                call close(actual,scalar(x,family,shape,2),2e-12_dp,'CDF')
                call check(actual >= previous .and. actual <= 0.0_dp,'monotone CDF')
                previous = actual
                one = make_censored_data(right=[x])
                call close(ll(one,family,shape),scalar(x,family,shape,3),2e-12_dp,'SF')
                one = make_censored_data(interval=reshape([x,nearest(x,1.0_dp)],[1,2]))
                expected = scalar(x,family,shape,1)+log(nearest(x,1.0_dp)-x)
                call close(ll(one,family,shape),expected,3e-13_dp,'adjacent interval')
            end do
            do i = 1, 9
                p = real(i,dp)/10.0_dp
                select case (family)
                case (1)
                    x = chi2_ppf(p,shape)
                case (2)
                    x = chi_ppf(p,shape)
                case (3)
                    x = nakagami_ppf(p,shape)
                end select
                one = make_censored_data(left=[x])
                call close(ll(one,family,shape),log(p),3e-12_dp,'CDF/PPF')
                select case (family)
                case (1)
                    x = chi2_isf(p,shape)
                case (2)
                    x = chi_isf(p,shape)
                case (3)
                    x = nakagami_isf(p,shape)
                end select
                one = make_censored_data(right=[x])
                call close(ll(one,family,shape),log(p),3e-12_dp,'SF/ISF')
            end do
        end do
        do k = -1, 1
            scale = 10.0_dp**(200*k)
            loc = -0.3_dp*scale
            x = 0.7_dp*scale
            data = make_censored_data(uncensored=[x])
            one = make_censored_data(uncensored=[1.0_dp])
            call close(ll(data,family,1.3_dp,loc,scale),ll(one,family,1.3_dp)-log(scale), &
                3e-13_dp,'affine density')
            data = make_censored_data(interval=reshape([x,nearest(x,1.0_dp)],[1,2]))
            expected = ll(one,family,1.3_dp)-log(scale)+log(nearest(x,1.0_dp)-x)
            call close(ll(data,family,1.3_dp,loc,scale),expected,3e-13_dp,'affine adjacent')
        end do
        call close(ll(empty,family,1.0_dp),0.0_dp,0.0_dp,'valid empty')
        call check(ieee_is_nan(ll(empty,family,0.0_dp)),'zero shape')
        call check(ieee_is_nan(ll(empty,family,-1.0_dp)),'negative shape')
        call check(ieee_is_nan(ll(empty,family,nan)),'NaN shape')
        call check(ieee_is_nan(ll(empty,family,inf)),'infinite shape')
        anchor = 1.0e6_dp
        if (family < 3) anchor = 2.0e6_dp
        call close(ll(empty,family,anchor),0.0_dp,0.0_dp,'largest supported shape')
        call check(ieee_is_nan(ll(empty,family,nearest(anchor,1.0_dp))),'shape above cap')
        call check(ieee_is_nan(ll(empty,family,huge(1.0_dp))),'huge shape rejected')
        call check(ieee_is_nan(ll(empty,family,1.0_dp,nan)),'NaN location')
        call check(ieee_is_nan(ll(empty,family,1.0_dp,scale=inf)),'infinite scale')
        call check(ieee_is_nan(ll(empty,family,1.0_dp,scale=0.0_dp)),'zero scale')
        one = make_censored_data(uncensored=[nan])
        call check(ieee_is_nan(ll(one,family,1.0_dp)),'NaN exact')
        one = make_censored_data(uncensored=[inf])
        call check(ieee_is_nan(ll(one,family,1.0_dp)),'infinite exact')
        one = make_censored_data(left=[-inf])
        call check(ieee_is_nan(ll(one,family,1.0_dp)),'invalid infinite left')
        one = make_censored_data(right=[inf])
        call check(ieee_is_nan(ll(one,family,1.0_dp)),'invalid infinite right')
        one = make_censored_data(interval=reshape([-inf,inf],[1,2]))
        call check(ieee_is_nan(ll(one,family,1.0_dp)),'invalid doubly infinite interval')
        one = make_censored_data(interval=reshape([0.1_dp,0.1_dp],[1,2]))
        call close(ll(one,family,1.0_dp),scalar(0.1_dp,family,1.0_dp,1), &
            2e-13_dp,'equal interval is exact')
        one = make_censored_data(left=[0.0_dp])
        call check(ll(one,family,1.0_dp) == -inf,'lower CDF')
        one = make_censored_data(interval=reshape([-inf,0.1_dp],[1,2]))
        call close(ll(one,family,1.0_dp),scalar(0.1_dp,family,1.0_dp,2), &
            2e-12_dp,'left-infinite interval')
        one = make_censored_data(interval=reshape([0.1_dp,inf],[1,2]))
        call close(ll(one,family,1.0_dp),scalar(0.1_dp,family,1.0_dp,3), &
            2e-12_dp,'right-infinite interval')
        one = make_censored_data(right=[0.0_dp])
        call close(ll(one,family,1.0_dp),0.0_dp,0.0_dp,'lower SF')
        one = make_censored_data(uncensored=[-1.0_dp])
        call check(ll(one,family,1.0_dp) == -inf,'outside density')
        one = make_censored_data(uncensored=[0.0_dp])
        anchor = 1.0_dp
        if (family == 1) anchor = 2.0_dp
        if (family == 3) anchor = 0.5_dp
        call close(ll(one,family,anchor),scalar(0.0_dp,family,anchor,1),2e-13_dp, &
            'finite endpoint density')
        call check(ll(one,family,anchor/2.0_dp) == inf,'singular endpoint')
        call check(ll(one,family,2.0_dp*anchor) == -inf,'zero endpoint density')
        tiny_shape = nearest(0.0_dp,1.0_dp)
        if (family < 3) then
            call check(ieee_is_nan(ll(empty,family,tiny_shape)),'half-shape underflow')
            tiny_shape = 2.0_dp*tiny_shape
        end if
        call close(ll(empty,family,tiny_shape),0.0_dp,0.0_dp,'minimum valid shape')
        one = make_censored_data(uncensored=[0.7_dp],left=[0.7_dp],right=[0.7_dp])
        call close(ll(one,family,tiny_shape),minimum_shape_mixed(family),3e-13_dp, &
            'minimum shape reference')
        anchor = 1.0e6_dp
        x = 1.0_dp
        if (family < 3) then
            anchor = 2.0e6_dp
            x = anchor
            if (family == 2) x = sqrt(anchor)
        end if
        one = make_censored_data(uncensored=[x],left=[x],right=[x])
        call close(ll(one,family,anchor),maximum_shape_center(family),2e-11_dp, &
            'maximum shape center reference')
        one = make_censored_data(interval=reshape([x,nearest(x,1.0_dp)],[1,2]))
        call close(ll(one,family,anchor),maximum_shape_adjacent(family),3e-13_dp, &
            'maximum shape adjacent reference')
        one = make_censored_data(left=[1e-200_dp])
        call check(ieee_is_finite(ll(one,family,tiny_shape)),'tiny-shape CDF')
        lower = [1.3_dp,0.0_dp,1.0_dp]
        call fit_family(mixed,family,lower,lower,fit)
        call check(fit%success,'fixed fit')
        call close(fit%nllf,-ll(mixed,family,1.3_dp),2e-13_dp,'fixed objective')
        select case (family)
        case (1)
            call close(chi2_censored_nnlf(mixed,1.3_dp),fit%nllf,0.0_dp,'NNLF')
            call chi2_fit([0.2_dp,0.4_dp,0.6_dp],lower,lower,fit)
        case (2)
            call close(chi_censored_nnlf(mixed,1.3_dp),fit%nllf,0.0_dp,'NNLF')
            call chi_fit([0.2_dp,0.4_dp,0.6_dp],lower,lower,fit)
        case (3)
            call close(nakagami_censored_nnlf(mixed,1.3_dp),fit%nllf,0.0_dp,'NNLF')
            call nakagami_fit([0.2_dp,0.4_dp,0.6_dp],lower,lower,fit)
        end select
        call check(fit%success,'array fit compatibility')
        lower = [0.3_dp,0.0_dp,1.3_dp]
        upper = [5.0_dp,0.0_dp,1.3_dp]
        call fit_family(mixed,family,lower,upper,fit)
        call check(fit%success,'shape fit')
        call close(fit%params(1),scipy_shape(family),3e-5_dp,'SciPy optimum shape')
        call close(fit%nllf,scipy_objective(family),2e-8_dp,'SciPy optimum objective')
        call check(fit%nllf <= -ll(mixed,family,1.3_dp,scale=1.3_dp),'fit improves')
        lower = [0.3_dp,-0.2_dp,0.5_dp]
        upper = [5.0_dp,0.05_dp,2.0_dp]
        expected = -ll(mixed,family,2.65_dp,-0.075_dp,1.25_dp)
        call fit_family(mixed,family,lower,upper,fit)
        call check(ieee_is_finite(fit%nllf),'finite all-free fit')
        call check(fit%nllf <= expected,'all-free improves midpoint')
        call fit_family(empty,family,lower,upper,fit)
        call check(fit%status == fit_status_invalid_input,'empty fit invalid')
    end do
    call close(ll(mixed,1,2.0_dp),exponential_censored_loglikelihood(mixed,scale=2.0_dp), &
        3e-13_dp,'chi2/exponential')
    call close(ll(mixed,1,3.4_dp),gamma_censored_loglikelihood(mixed,1.7_dp,scale=2.0_dp), &
        3e-13_dp,'chi2/gamma')
    call close(ll(mixed,2,1.0_dp),halfnorm_censored_loglikelihood(mixed),3e-13_dp,'chi/half-normal')
    call close(ll(mixed,2,2.0_dp),rayleigh_censored_loglikelihood(mixed),3e-13_dp,'chi/Rayleigh')
    call close(ll(mixed,2,3.0_dp),maxwell_censored_loglikelihood(mixed),3e-13_dp,'chi/Maxwell')
    call close(ll(mixed,3,0.5_dp),halfnorm_censored_loglikelihood(mixed),3e-13_dp,'Nakagami/half-normal')
    call close(ll(mixed,3,1.0_dp),rayleigh_censored_loglikelihood(mixed,scale=sqrt(0.5_dp)), &
        3e-13_dp,'Nakagami/Rayleigh')
    call close(nakagami_logsf(1e200_dp,1e-300_dp)/1e100_dp,-1.0_dp, &
        3e-13_dp,'tiny-shape scalar overflow guard')
    one = make_censored_data(right=[1e200_dp])
    call close(ll(one,3,1e-300_dp)/1e100_dp,-1.0_dp,3e-13_dp,'tiny-shape far tail')
    one = make_censored_data(left=[nearest(0.0_dp,1.0_dp)])
    call close(ll(one,2,1.0_dp),log(nearest(0.0_dp,1.0_dp)) &
        +0.5_dp*log(2.0_dp/acos(-1.0_dp)),3e-13_dp,'subnormal half-normal CDF')
    one = make_censored_data(interval=reshape([0.0_dp,nearest(0.0_dp,1.0_dp)],[1,2]))
    call close(ll(one,2,1.0_dp),log(nearest(0.0_dp,1.0_dp)) &
        +0.5_dp*log(2.0_dp/acos(-1.0_dp)),3e-13_dp,'subnormal half-normal interval')
    write(*,'(a,i0,a,i0)') 'Gamma-power censored checks: ',checks,', failures: ',failures
    if (failures /= 0) error stop 1
contains
    pure function ll(data, family, shape, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        integer, intent(in) :: family !! chi2=1, chi=2, Nakagami=3
        real(dp), intent(in) :: shape !! positive family shape
        real(dp), intent(in), optional :: loc !! finite support origin, default zero
        real(dp), intent(in), optional :: scale !! positive scale, default one
        real(dp) :: value
        select case (family)
        case (1)
            value = chi2_censored_loglikelihood(data,shape,loc,scale)
        case (2)
            value = chi_censored_loglikelihood(data,shape,loc,scale)
        case (3)
            value = nakagami_censored_loglikelihood(data,shape,loc,scale)
        end select
    end function ll

    pure function scalar(x, family, shape, part) result(value)
        implicit none
        real(dp), intent(in) :: x !! scalar observation or limit
        integer, intent(in) :: family !! chi2=1, chi=2, Nakagami=3
        real(dp), intent(in) :: shape !! positive family shape
        integer, intent(in) :: part !! density=1, CDF=2, SF=3
        real(dp) :: value
        select case (family)
        case (1)
            select case (part)
            case (1)
                value = chi2_logpdf(x,shape)
            case (2)
                value = chi2_logcdf(x,shape)
            case (3)
                value = chi2_logsf(x,shape)
            end select
        case (2)
            select case (part)
            case (1)
                value = chi_logpdf(x,shape)
            case (2)
                value = chi_logcdf(x,shape)
            case (3)
                value = chi_logsf(x,shape)
            end select
        case (3)
            select case (part)
            case (1)
                value = nakagami_logpdf(x,shape)
            case (2)
                value = nakagami_logcdf(x,shape)
            case (3)
                value = nakagami_logsf(x,shape)
            end select
        end select
    end function scalar

    subroutine fit_family(data, family, lower, upper, fit)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        integer, intent(in) :: family !! chi2=1, chi=2, Nakagami=3
        real(dp), intent(in) :: lower(3) !! shape, location, scale lower bounds
        real(dp), intent(in) :: upper(3) !! upper bounds; equality fixes a parameter
        type(fit_result), intent(out) :: fit !! fit and diagnostics
        select case (family)
        case (1)
            call chi2_fit(data,lower,upper,fit)
        case (2)
            call chi_fit(data,lower,upper,fit)
        case (3)
            call nakagami_fit(data,lower,upper,fit)
        end select
    end subroutine fit_family

    subroutine check(ok, label)
        implicit none
        logical, intent(in) :: ok !! true when the assertion passes
        character(*), intent(in) :: label !! assertion description
        checks = checks+1
        if (ok) return
        failures = failures+1
        write(*,*) 'FAIL: ',label
    end subroutine check

    subroutine close(actual, expected, tol, label)
        implicit none
        real(dp), intent(in) :: actual !! observed result
        real(dp), intent(in) :: expected !! reference value
        real(dp), intent(in) :: tol !! relative tolerance with unit absolute floor
        character(*), intent(in) :: label !! assertion description
        logical :: ok
        ok = actual == expected
        if (ieee_is_finite(actual) .and. ieee_is_finite(expected)) &
            ok = abs(actual-expected) <= tol*max(1.0_dp,abs(expected))
        call check(ok,label)
        if (.not. ok) write(*,*) actual,expected
    end subroutine close
end program test_censored_radial_gamma
