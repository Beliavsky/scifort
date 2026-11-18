! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_censored_noncentral
    use scifort_stats, only : dp, censored_data, make_censored_data, fit_result, &
        fit_status_invalid_input, chi2_censored_loglikelihood, &
        rayleigh_censored_loglikelihood, f_censored_loglikelihood, &
        ncx2_censored_loglikelihood, ncx2_censored_nnlf, ncx2_fit, &
        ncx2_logpdf, ncx2_logcdf, ncx2_logsf, ncx2_ppf, ncx2_isf, &
        rice_censored_loglikelihood, rice_censored_nnlf, rice_fit, &
        rice_logpdf, rice_logcdf, rice_logsf, rice_ppf, rice_isf, &
        ncf_censored_loglikelihood, ncf_censored_nnlf, ncf_fit, &
        ncf_logpdf, ncf_logcdf, ncf_logsf, ncf_ppf, ncf_isf
    use scifort_math, only : positive_infinity, quiet_nan
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    implicit none
    type(censored_data) :: data, one, empty
    type(fit_result) :: fit, old
    real(dp) :: par(6), lower(6), upper(6), x, actual, expected, previous, p, inf, nan
    integer :: checks, failures, family, i, j, k, n
    checks = 0
    failures = 0
    inf = positive_infinity(0.0_dp)
    nan = quiet_nan(0.0_dp)
    empty = make_censored_data()
    data = make_censored_data(uncensored=[0.2_dp,0.4_dp,0.6_dp],left=[0.1_dp], &
        right=[0.8_dp],interval=reshape([0.3_dp,0.5_dp,0.35_dp,0.55_dp],[2,2]))
    do family = 1, 3
        n = 4
        if (family == 2) n = 3
        if (family == 3) n = 5
        do j = 1, 4
            call parameters(family,par)
            if (family == 1) par(1) = real(j,dp)/2.0_dp
            if (family == 2) par(1) = real(j,dp)/2.0_dp
            if (family == 3) par(1) = real(j,dp)/2.0_dp
            previous = -inf
            do i = 1, 25
                x = 0.05_dp+real(i,dp)/10.0_dp
                one = make_censored_data(uncensored=[x])
                call close(ll(one,family,par),scalar(x,family,par,1),3e-12_dp,'density')
                one = make_censored_data(left=[x])
                actual = ll(one,family,par)
                call close(actual,scalar(x,family,par,2),3e-12_dp,'CDF')
                call check(actual >= previous .and. actual <= 0.0_dp,'CDF monotonicity')
                previous = actual
                one = make_censored_data(right=[x])
                call close(ll(one,family,par),scalar(x,family,par,3),3e-12_dp,'SF')
                one = make_censored_data(interval=reshape([x,nearest(x,1.0_dp)],[1,2]))
                expected = scalar(x,family,par,1)+log(nearest(x,1.0_dp)-x)
                call close(ll(one,family,par),expected,3e-12_dp,'adjacent interval')
                one = make_censored_data(interval=reshape([x,x+0.1_dp],[1,2]))
                if (scalar(x,family,par,2) < log(0.5_dp)) then
                    expected = log(exp(scalar(x+0.1_dp,family,par,2))-exp(scalar(x,family,par,2)))
                else
                    expected = log(exp(scalar(x,family,par,3))-exp(scalar(x+0.1_dp,family,par,3)))
                end if
                call close(ll(one,family,par),expected,3e-11_dp,'ordinary interval')
            end do
            do i = 1, 9
                p = real(i,dp)/10.0_dp
                x = quantile(p,family,par,.false.)
                one = make_censored_data(left=[x])
                call close(ll(one,family,par),log(p),5e-11_dp,'CDF/PPF')
                x = quantile(p,family,par,.true.)
                one = make_censored_data(right=[x])
                call close(ll(one,family,par),log(p),5e-11_dp,'SF/ISF')
            end do
        end do
        call parameters(family,par)
        do k = -1, 1
            par(n) = 10.0_dp**(200*k)
            par(n-1) = -0.3_dp*par(n)
            x = 0.3_dp*par(n)
            one = make_censored_data(uncensored=[x])
            expected = scalar(0.6_dp,family,par,1)-log(par(n))
            call close(ll(one,family,par),expected,3e-12_dp,'affine density')
            one = make_censored_data(interval=reshape([x,nearest(x,1.0_dp)],[1,2]))
            expected = expected+log(nearest(x,1.0_dp)-x)
            call close(ll(one,family,par),expected,3e-12_dp,'affine adjacent interval')
        end do
        call parameters(family,par)
        call close(ll(empty,family,par),0.0_dp,0.0_dp,'valid empty')
        par(n) = 0.0_dp
        call check(ieee_is_nan(ll(empty,family,par)),'invalid empty scale')
        par(n) = inf
        call check(ieee_is_nan(ll(empty,family,par)),'infinite scale')
        call parameters(family,par)
        par(1) = nan
        call check(ieee_is_nan(ll(empty,family,par)),'NaN shape')
        call parameters(family,par)
        one = make_censored_data(uncensored=[inf])
        call check(ieee_is_nan(ll(one,family,par)),'infinite exact invalid')
        one = make_censored_data(left=[nan])
        call check(ieee_is_nan(ll(one,family,par)),'NaN data invalid')
        one = make_censored_data(interval=reshape([1.0_dp,0.9_dp],[1,2]))
        call check(ieee_is_nan(ll(one,family,par)),'reversed interval')
        one = make_censored_data(interval=reshape([-inf,inf],[1,2]))
        call check(ieee_is_nan(ll(one,family,par)),'whole infinite interval invalid')
        one = make_censored_data(interval=reshape([0.7_dp,0.7_dp],[1,2]))
        call close(ll(one,family,par),scalar(0.7_dp,family,par,1),3e-12_dp,'equal endpoints exact')
        one = make_censored_data(interval=reshape([-inf,0.7_dp],[1,2]))
        call close(ll(one,family,par),scalar(0.7_dp,family,par,2),3e-12_dp,'left interval')
        one = make_censored_data(interval=reshape([0.7_dp,inf],[1,2]))
        call close(ll(one,family,par),scalar(0.7_dp,family,par,3),3e-12_dp,'right interval')
        lower = par
        upper = par
        call fit_family(data,family,lower,upper,fit)
        call check(fit%success,'fixed fit converged')
        call close(fit%nllf,-ll(data,family,par),3e-13_dp,'fixed NNLF')
        call fit_family(empty,family,lower,upper,fit)
        call check(fit%status == fit_status_invalid_input,'empty fit rejected')
        one = make_censored_data(uncensored=[0.2_dp,0.4_dp,0.6_dp])
        call fit_family(one,family,lower,upper,fit)
        call fit_array([0.2_dp,0.4_dp,0.6_dp],family,lower,upper,old)
        call check(fit%success .and. old%success,'old generic retained')
        call close(fit%nllf,old%nllf,3e-12_dp,'array/exact fit agreement')
        lower(n) = 0.3_dp
        upper(n) = 3.0_dp
        call fit_family(data,family,lower,upper,fit)
        call check(fit%success,'one-parameter fit converges')
        call check(fit%nllf <= -ll(data,family,par),'one-parameter fit improves')
        select case (family)
        case (1)
            expected = 10.974656410900508_dp
        case (2)
            expected = 10.055647712310694_dp
        case (3)
            expected = 10.612152223011206_dp
        end select
        call close(fit%params(n),0.3000000072924005_dp,3e-5_dp,'independent scale fit')
        call close(fit%nllf,expected,2e-8_dp,'independent fitted NNLF')
        lower = par-0.1_dp
        upper = par+0.1_dp
        select case (family)
        case (1)
            call ncx2_fit(data,lower(:4),upper(:4),fit,guess=par(:4),max_iter=20)
        case (2)
            call rice_fit(data,lower(:3),upper(:3),fit,guess=par(:3),max_iter=20)
        case (3)
            call ncf_fit(data,lower(:5),upper(:5),fit,guess=par(:5),max_iter=20)
        end select
        call check(ieee_is_finite(fit%nllf),'all-free finite candidate')
        call check(fit%nllf <= -ll(data,family,par),'all-free improves initial point')
    end do
    call close(ncx2_censored_loglikelihood(data,2.0_dp,0.0_dp), &
        chi2_censored_loglikelihood(data,2.0_dp),3e-12_dp,'central chi-squared identity')
    call close(rice_censored_loglikelihood(data,0.0_dp), &
        rayleigh_censored_loglikelihood(data),3e-12_dp,'Rayleigh identity')
    call close(ncf_censored_loglikelihood(data,2.0_dp,3.0_dp,0.0_dp), &
        f_censored_loglikelihood(data,2.0_dp,3.0_dp),3e-12_dp,'central F identity')
    one = make_censored_data(left=[0.0_dp])
    call close(ncx2_censored_loglikelihood(one,2.0_dp,1.0_dp),-inf,0.0_dp,'support CDF')
    one = make_censored_data(right=[0.0_dp])
    call close(rice_censored_loglikelihood(one,1.0_dp),0.0_dp,0.0_dp,'support SF')
    one = make_censored_data(right=[1500.0_dp])
    call close(ncx2_censored_loglikelihood(one,2.0_dp,0.0_dp),-750.0_dp,3e-13_dp,'central tail')
    one = make_censored_data(left=[1400.0_dp])
    expected = -exp(-700.0_dp)
    call close(ncx2_censored_loglikelihood(one,2.0_dp,0.0_dp)/expected,1.0_dp, &
        3e-12_dp,'tiny log-CDF relative accuracy')
    empty = make_censored_data()
    call check(ieee_is_nan(ncx2_censored_loglikelihood(empty,0.01_dp,1.0_dp)),'df range')
    call check(ieee_is_nan(ncf_censored_loglikelihood(empty,2.0_dp,2.0_dp,101.0_dp)),'nc range')
    call check(ieee_is_nan(rice_censored_loglikelihood(empty,-1.0_dp)),'Rice range')
    one = make_censored_data(uncensored=[0.0_dp])
    call close(ncx2_censored_loglikelihood(one,2.0_dp,2.0_dp),-1.0_dp-log(2.0_dp), &
        3e-13_dp,'endpoint finite density')
    call close(ncx2_censored_loglikelihood(one,1.0_dp,2.0_dp),inf,0.0_dp,'endpoint infinite density')
    call close(ncx2_censored_loglikelihood(one,3.0_dp,2.0_dp),-inf,0.0_dp,'endpoint zero density')
    call close(rice_censored_loglikelihood(one,1.0_dp),-inf,0.0_dp,'Rice endpoint density')
    one = make_censored_data(right=[1e20_dp])
    call check(ieee_is_nan(ncx2_censored_loglikelihood(one,2.0_dp,1.0_dp)), &
        'unconverged series returns NaN')
    one = make_censored_data(left=[huge(1.0_dp)])
    expected = ncx2_logcdf(2.0_dp,2.0_dp,1.0_dp)
    call close(ncx2_censored_loglikelihood(one,2.0_dp,1.0_dp, &
        loc=-huge(1.0_dp),scale=huge(1.0_dp)),expected,3e-12_dp,'overflowed physical distance')
    data = make_censored_data(uncensored=[0.2_dp,0.4_dp,0.6_dp],left=[0.1_dp],right=[0.8_dp])
    call close(ncx2_censored_loglikelihood(data,2.0_dp,0.0_dp), &
        -ncx2_censored_nnlf(data,2.0_dp,0.0_dp),0.0_dp,'ncx2 NNLF')
    call close(rice_censored_loglikelihood(data,0.0_dp), &
        -rice_censored_nnlf(data,0.0_dp),0.0_dp,'Rice NNLF')
    call close(ncf_censored_loglikelihood(data,2.0_dp,3.0_dp,0.0_dp), &
        -ncf_censored_nnlf(data,2.0_dp,3.0_dp,0.0_dp),0.0_dp,'ncf NNLF')
    x = nearest(0.0_dp,1.0_dp)
    one = make_censored_data(interval=reshape([x,2.0_dp*x],[1,2]))
    call close(ncf_censored_loglikelihood(one,2.0_dp,3.0_dp,1.0_dp), &
        -0.5_dp+log(x),3e-13_dp,'subnormal beta interval')
    one = make_censored_data(left=[x])
    call close(ncf_censored_loglikelihood(one,2.0_dp,3.0_dp,1.0_dp), &
        -0.5_dp+log(x),3e-13_dp,'subnormal beta CDF')
    one = make_censored_data(uncensored=[x])
    call close(ncf_censored_loglikelihood(one,2.0_dp,3.0_dp,1.0_dp), &
        -0.5_dp,3e-13_dp,'subnormal beta density')
    one = make_censored_data(right=[23.407231640316144_dp])
    call close(rice_censored_loglikelihood(one,5.975890905275955_dp), &
        -155.0215267911996304_dp,3e-13_dp,'Rice tail multiprecision reference')
    one = make_censored_data(left=[0.0009567362136092416_dp])
    call close(rice_censored_loglikelihood(one,9.134584309754336_dp), &
        -56.3174188050420075_dp,3e-13_dp,'Rice lower multiprecision reference')
    one = make_censored_data(uncensored=[298394586479.9011_dp])
    call close(ncf_censored_loglikelihood(one,481.6793871540879_dp, &
        0.8139110183831583_dp,58.75326915349704_dp),-38.2727002183747198_dp, &
        3e-13_dp,'noncentral F density multiprecision reference')
    call check(ieee_is_nan(ncf_censored_loglikelihood(empty, &
        1.0_dp,1.0_dp+epsilon(1.0_dp),0.0_dp,scale=huge(1.0_dp))), &
        'transformed scale overflow rejected')
    call check(ieee_is_nan(ncf_censored_loglikelihood(empty, &
        2.0_dp,1.0_dp,0.0_dp,scale=tiny(1.0_dp))), &
        'subnormal transformed scale rejected')
    write(*,'(a,i0,a,i0)') 'noncentral censoring: ',checks,' checks; failures: ',failures
    if (failures /= 0) error stop 1
contains
    subroutine parameters(family,par)
        implicit none
        integer, intent(in) :: family !! ncx2=1, Rice=2, noncentral F=3
        real(dp), intent(out) :: par(6) !! shapes, loc, scale with unused entries zero
        par = 0.0_dp
        select case (family)
        case (1)
            par(:4) = [2.0_dp,1.0_dp,0.0_dp,1.0_dp]
        case (2)
            par(:3) = [1.0_dp,0.0_dp,1.0_dp]
        case (3)
            par(:5) = [2.0_dp,3.0_dp,1.0_dp,0.0_dp,1.0_dp]
        end select
    end subroutine parameters

    pure function ll(data,family,par) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        integer, intent(in) :: family !! ncx2=1, Rice=2, noncentral F=3
        real(dp), intent(in) :: par(6) !! shapes, loc, scale; unused tail ignored
        real(dp) :: value
        select case (family)
        case (1)
            value = ncx2_censored_loglikelihood(data,par(1),par(2),par(3),par(4))
        case (2)
            value = rice_censored_loglikelihood(data,par(1),par(2),par(3))
        case (3)
            value = ncf_censored_loglikelihood(data,par(1),par(2),par(3),par(4),par(5))
        end select
    end function ll

    pure function scalar(x,family,par,part) result(value)
        implicit none
        real(dp), intent(in) :: x !! standardized observation
        integer, intent(in) :: part !! density=1, CDF=2, SF=3
        integer, intent(in) :: family !! ncx2=1, Rice=2, noncentral F=3
        real(dp), intent(in) :: par(6) !! shapes, loc, scale; unused tail ignored
        real(dp) :: value
        select case (family)
        case (1)
            select case (part)
            case (1)
                value = ncx2_logpdf(x,par(1),par(2))
            case (2)
                value = ncx2_logcdf(x,par(1),par(2))
            case (3)
                value = ncx2_logsf(x,par(1),par(2))
            end select
        case (2)
            select case (part)
            case (1)
                value = rice_logpdf(x,par(1))
            case (2)
                value = rice_logcdf(x,par(1))
            case (3)
                value = rice_logsf(x,par(1))
            end select
        case (3)
            select case (part)
            case (1)
                value = ncf_logpdf(x,par(1),par(2),par(3))
            case (2)
                value = ncf_logcdf(x,par(1),par(2),par(3))
            case (3)
                value = ncf_logsf(x,par(1),par(2),par(3))
            end select
        end select
    end function scalar

    pure function quantile(p,family,par,upper) result(value)
        implicit none
        real(dp), intent(in) :: p !! probability in (0,1)
        logical, intent(in) :: upper !! select survival quantile
        integer, intent(in) :: family !! ncx2=1, Rice=2, noncentral F=3
        real(dp), intent(in) :: par(6) !! shapes, loc, scale; unused tail ignored
        real(dp) :: value
        select case (family)
        case (1)
            if (upper) then
                value = ncx2_isf(p,par(1),par(2))
            else
                value = ncx2_ppf(p,par(1),par(2))
            end if
        case (2)
            if (upper) then
                value = rice_isf(p,par(1))
            else
                value = rice_ppf(p,par(1))
            end if
        case (3)
            if (upper) then
                value = ncf_isf(p,par(1),par(2),par(3))
            else
                value = ncf_ppf(p,par(1),par(2),par(3))
            end if
        end select
    end function quantile

    subroutine fit_family(data,family,lower,upper,result)
        implicit none
        type(censored_data), intent(in) :: data !! observations
        real(dp), intent(in) :: lower(6), upper(6) !! bounds in shape, loc, scale order
        type(fit_result), intent(out) :: result !! fit and diagnostics
        integer, intent(in) :: family !! ncx2=1, Rice=2, noncentral F=3
        select case (family)
        case (1)
            call ncx2_fit(data,lower(:4),upper(:4),result)
        case (2)
            call rice_fit(data,lower(:3),upper(:3),result)
        case (3)
            call ncf_fit(data,lower(:5),upper(:5),result)
        end select
    end subroutine fit_family

    subroutine fit_array(data,family,lower,upper,result)
        implicit none
        real(dp), intent(in) :: data(:) !! exact observations
        real(dp), intent(in) :: lower(6), upper(6) !! bounds in shape, loc, scale order
        type(fit_result), intent(out) :: result !! fit and diagnostics
        integer, intent(in) :: family !! ncx2=1, Rice=2, noncentral F=3
        select case (family)
        case (1)
            call ncx2_fit(data,lower(:4),upper(:4),result)
        case (2)
            call rice_fit(data,lower(:3),upper(:3),result)
        case (3)
            call ncf_fit(data,lower(:5),upper(:5),result)
        end select
    end subroutine fit_array

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
        real(dp), intent(in) :: tol !! scaled tolerance with unit absolute floor
        character(*), intent(in) :: label !! assertion description
        logical :: ok
        ok = actual == expected
        if (ieee_is_finite(actual) .and. ieee_is_finite(expected)) &
            ok = abs(actual-expected) <= tol*max(1.0_dp,abs(expected))
        call check(ok,label)
        if (.not. ok) write(*,*) actual,expected
    end subroutine close
end program test_censored_noncentral
