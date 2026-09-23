! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Fisher's noncentral hypergeometric distribution matching
! scipy.stats.nchypergeom_fisher.  The finite support is normalized directly
! in log space, so no external biased-urn library is required.
module scifort_nchypergeom_fisher
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : log1p_safe, negative_infinity, quiet_nan
    implicit none
    private

    public :: nchypergeom_fisher_cdf, nchypergeom_fisher_isf, nchypergeom_fisher_logcdf
    public :: nchypergeom_fisher_logpmf, nchypergeom_fisher_logsf, nchypergeom_fisher_pmf
    public :: nchypergeom_fisher_ppf, nchypergeom_fisher_sf
    public :: nchypergeom_fisher_logpmf_derivative_odds

    interface nchypergeom_fisher_pmf
        module procedure nchypergeom_fisher_pmf_real
        module procedure nchypergeom_fisher_pmf_int
    end interface
    interface nchypergeom_fisher_logpmf
        module procedure nchypergeom_fisher_logpmf_real
        module procedure nchypergeom_fisher_logpmf_int
    end interface
    interface nchypergeom_fisher_cdf
        module procedure nchypergeom_fisher_cdf_real
        module procedure nchypergeom_fisher_cdf_int
    end interface
    interface nchypergeom_fisher_sf
        module procedure nchypergeom_fisher_sf_real
        module procedure nchypergeom_fisher_sf_int
    end interface
    interface nchypergeom_fisher_logcdf
        module procedure nchypergeom_fisher_logcdf_real
        module procedure nchypergeom_fisher_logcdf_int
    end interface
    interface nchypergeom_fisher_logsf
        module procedure nchypergeom_fisher_logsf_real
        module procedure nchypergeom_fisher_logsf_int
    end interface

contains

    pure elemental function nchypergeom_fisher_logpmf_real(k, m, n, draws, odds, loc) result(y)
        real(dp), intent(in) :: k !! observed Type-I count
        real(dp), intent(in) :: m !! positive integer population size
        real(dp), intent(in) :: n !! integer number of Type-I objects in [0,m]
        real(dp), intent(in) :: draws !! integer sample size in [0,m]
        real(dp), intent(in) :: odds !! strictly positive Type-I odds ratio
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y, count, lower, upper
        count = shifted_count(k, loc)
        if (.not. valid_parameters(m, n, draws, odds) .or. ieee_is_nan(count)) then
            y = quiet_nan(k)
            return
        end if
        lower = support_lower(m, n, draws); upper = support_upper(n, draws)
        if (.not. ieee_is_finite(count) .or. count /= aint(count) .or. count < lower .or. count > upper) then
            y = negative_infinity(k)
        else
            y = log_weight(int(count), m, n, draws, odds) - log_normalizer(m, n, draws, odds)
        end if
    end function nchypergeom_fisher_logpmf_real

    pure elemental function nchypergeom_fisher_logpmf_int(k, m, n, draws, odds, loc) result(y)
        integer, intent(in) :: k !! observed Type-I count
        integer, intent(in) :: m !! positive population size
        integer, intent(in) :: n !! number of Type-I objects
        integer, intent(in) :: draws !! sample size
        real(dp), intent(in) :: odds !! strictly positive Type-I odds ratio
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = nchypergeom_fisher_logpmf_real(real(k,dp), real(m,dp), real(n,dp), real(draws,dp), odds, loc)
    end function nchypergeom_fisher_logpmf_int

    pure elemental function nchypergeom_fisher_pmf_real(k, m, n, draws, odds, loc) result(y)
        real(dp), intent(in) :: k !! observed Type-I count
        real(dp), intent(in) :: m !! positive integer population size
        real(dp), intent(in) :: n !! integer number of Type-I objects
        real(dp), intent(in) :: draws !! integer sample size
        real(dp), intent(in) :: odds !! strictly positive Type-I odds ratio
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = exp(nchypergeom_fisher_logpmf_real(k,m,n,draws,odds,loc))
    end function nchypergeom_fisher_pmf_real

    pure elemental function nchypergeom_fisher_pmf_int(k, m, n, draws, odds, loc) result(y)
        integer, intent(in) :: k !! observed Type-I count
        integer, intent(in) :: m !! positive population size
        integer, intent(in) :: n !! number of Type-I objects
        integer, intent(in) :: draws !! sample size
        real(dp), intent(in) :: odds !! strictly positive Type-I odds ratio
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = nchypergeom_fisher_pmf_real(real(k,dp),real(m,dp),real(n,dp),real(draws,dp),odds,loc)
    end function nchypergeom_fisher_pmf_int

    pure elemental function nchypergeom_fisher_cdf_real(k, m, n, draws, odds, loc) result(y)
        real(dp), intent(in) :: k !! upper limit of P(X <= k)
        real(dp), intent(in) :: m !! positive integer population size
        real(dp), intent(in) :: n !! integer number of Type-I objects
        real(dp), intent(in) :: draws !! integer sample size
        real(dp), intent(in) :: odds !! strictly positive Type-I odds ratio
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y, sf, logcdf, logsf
        call tails(floor_count(k,loc),m,n,draws,odds,y,sf,logcdf,logsf)
    end function nchypergeom_fisher_cdf_real

    pure elemental function nchypergeom_fisher_cdf_int(k, m, n, draws, odds, loc) result(y)
        integer, intent(in) :: k !! upper integer limit
        integer, intent(in) :: m !! positive population size
        integer, intent(in) :: n !! number of Type-I objects
        integer, intent(in) :: draws !! sample size
        real(dp), intent(in) :: odds !! strictly positive Type-I odds ratio
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = nchypergeom_fisher_cdf_real(real(k,dp),real(m,dp),real(n,dp),real(draws,dp),odds,loc)
    end function nchypergeom_fisher_cdf_int

    pure elemental function nchypergeom_fisher_sf_real(k, m, n, draws, odds, loc) result(y)
        real(dp), intent(in) :: k !! lower limit of P(X > k)
        real(dp), intent(in) :: m !! positive integer population size
        real(dp), intent(in) :: n !! integer number of Type-I objects
        real(dp), intent(in) :: draws !! integer sample size
        real(dp), intent(in) :: odds !! strictly positive Type-I odds ratio
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y, cdf, logcdf, logsf
        call tails(floor_count(k,loc),m,n,draws,odds,cdf,y,logcdf,logsf)
    end function nchypergeom_fisher_sf_real

    pure elemental function nchypergeom_fisher_sf_int(k, m, n, draws, odds, loc) result(y)
        integer, intent(in) :: k !! lower integer limit
        integer, intent(in) :: m !! positive population size
        integer, intent(in) :: n !! number of Type-I objects
        integer, intent(in) :: draws !! sample size
        real(dp), intent(in) :: odds !! strictly positive Type-I odds ratio
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = nchypergeom_fisher_sf_real(real(k,dp),real(m,dp),real(n,dp),real(draws,dp),odds,loc)
    end function nchypergeom_fisher_sf_int

    pure elemental function nchypergeom_fisher_logcdf_real(k, m, n, draws, odds, loc) result(y)
        real(dp), intent(in) :: k !! upper limit of log(P(X <= k))
        real(dp), intent(in) :: m !! positive integer population size
        real(dp), intent(in) :: n !! integer number of Type-I objects
        real(dp), intent(in) :: draws !! integer sample size
        real(dp), intent(in) :: odds !! strictly positive Type-I odds ratio
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y, cdf, sf, logsf
        call tails(floor_count(k,loc),m,n,draws,odds,cdf,sf,y,logsf)
    end function nchypergeom_fisher_logcdf_real

    pure elemental function nchypergeom_fisher_logcdf_int(k, m, n, draws, odds, loc) result(y)
        integer, intent(in) :: k !! upper integer limit
        integer, intent(in) :: m !! positive population size
        integer, intent(in) :: n !! number of Type-I objects
        integer, intent(in) :: draws !! sample size
        real(dp), intent(in) :: odds !! strictly positive Type-I odds ratio
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y=nchypergeom_fisher_logcdf_real(real(k,dp),real(m,dp),real(n,dp),real(draws,dp),odds,loc)
    end function nchypergeom_fisher_logcdf_int

    pure elemental function nchypergeom_fisher_logsf_real(k, m, n, draws, odds, loc) result(y)
        real(dp), intent(in) :: k !! lower limit of log(P(X > k))
        real(dp), intent(in) :: m !! positive integer population size
        real(dp), intent(in) :: n !! integer number of Type-I objects
        real(dp), intent(in) :: draws !! integer sample size
        real(dp), intent(in) :: odds !! strictly positive Type-I odds ratio
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y, cdf, sf, logcdf
        call tails(floor_count(k,loc),m,n,draws,odds,cdf,sf,logcdf,y)
    end function nchypergeom_fisher_logsf_real

    pure elemental function nchypergeom_fisher_logsf_int(k, m, n, draws, odds, loc) result(y)
        integer, intent(in) :: k !! lower integer limit
        integer, intent(in) :: m !! positive population size
        integer, intent(in) :: n !! number of Type-I objects
        integer, intent(in) :: draws !! sample size
        real(dp), intent(in) :: odds !! strictly positive Type-I odds ratio
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y=nchypergeom_fisher_logsf_real(real(k,dp),real(m,dp),real(n,dp),real(draws,dp),odds,loc)
    end function nchypergeom_fisher_logsf_int

    pure elemental function nchypergeom_fisher_ppf(probability, m, n, draws, odds, loc) result(y)
        real(dp), intent(in) :: probability !! lower-tail probability in [0,1]
        real(dp), intent(in) :: m !! positive integer population size
        real(dp), intent(in) :: n !! integer number of Type-I objects
        real(dp), intent(in) :: draws !! integer sample size
        real(dp), intent(in) :: odds !! strictly positive Type-I odds ratio
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y, shift, lower, upper
        integer :: lo, hi, mid
        shift = optional_loc(loc)
        if (.not. valid_parameters(m,n,draws,odds) .or. .not. ieee_is_finite(shift) .or. &
            .not. valid_probability(probability)) then
            y = quiet_nan(probability); return
        end if
        lower=support_lower(m,n,draws); upper=support_upper(n,draws)
        if (probability <= 0.0_dp) then
            y=shift+lower-1.0_dp
        else if (probability >= 1.0_dp) then
            y=shift+upper
        else
            lo=int(lower); hi=int(upper)
            do while (lo < hi)
                mid=lo+(hi-lo)/2
                if (nchypergeom_fisher_cdf_real(shift+real(mid,dp),m,n,draws,odds,shift) >= probability) then
                    hi=mid
                else
                    lo=mid+1
                end if
            end do
            y=shift+real(lo,dp)
        end if
    end function nchypergeom_fisher_ppf

    pure elemental function nchypergeom_fisher_isf(probability, m, n, draws, odds, loc) result(y)
        real(dp), intent(in) :: probability !! upper-tail probability in [0,1]
        real(dp), intent(in) :: m !! positive integer population size
        real(dp), intent(in) :: n !! integer number of Type-I objects
        real(dp), intent(in) :: draws !! integer sample size
        real(dp), intent(in) :: odds !! strictly positive Type-I odds ratio
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y, shift, lower, upper
        integer :: lo, hi, mid
        shift=optional_loc(loc)
        if (.not. valid_parameters(m,n,draws,odds) .or. .not. ieee_is_finite(shift) .or. &
            .not. valid_probability(probability)) then
            y=quiet_nan(probability); return
        end if
        lower=support_lower(m,n,draws); upper=support_upper(n,draws)
        if (probability >= 1.0_dp) then
            y=shift+lower-1.0_dp
        else if (probability <= 0.0_dp) then
            y=shift+upper
        else
            lo=int(lower); hi=int(upper)
            do while (lo < hi)
                mid=lo+(hi-lo)/2
                if (nchypergeom_fisher_sf_real(shift+real(mid,dp),m,n,draws,odds,shift) <= probability) then
                    hi=mid
                else
                    lo=mid+1
                end if
            end do
            y=shift+real(lo,dp)
        end if
    end function nchypergeom_fisher_isf

    pure subroutine nchypergeom_fisher_logpmf_derivative_odds(count, m, n, draws, odds, logf, dodds)
        real(dp), intent(in) :: count !! unshifted integer observation
        real(dp), intent(in) :: m !! positive integer population size
        real(dp), intent(in) :: n !! integer number of Type-I objects
        real(dp), intent(in) :: draws !! integer sample size
        real(dp), intent(in) :: odds !! strictly positive Type-I odds ratio
        real(dp), intent(out) :: logf !! log probability mass
        real(dp), intent(out) :: dodds !! derivative of log mass with respect to odds
        logf = nchypergeom_fisher_logpmf_real(count,m,n,draws,odds)
        if (.not. ieee_is_finite(logf)) then
            dodds = quiet_nan(count)
        else
            dodds = (count - expected_count(m,n,draws,odds)) / odds
        end if
    end subroutine nchypergeom_fisher_logpmf_derivative_odds

    pure subroutine tails(k,m,n,draws,odds,cdf,sf,logcdf,logsf)
        real(dp), intent(in) :: k !! shifted integer tail split
        real(dp), intent(in) :: m !! positive integer population size
        real(dp), intent(in) :: n !! integer number of Type-I objects
        real(dp), intent(in) :: draws !! integer sample size
        real(dp), intent(in) :: odds !! strictly positive Type-I odds ratio
        real(dp), intent(out) :: cdf !! lower-tail probability
        real(dp), intent(out) :: sf !! upper-tail probability
        real(dp), intent(out) :: logcdf !! logarithm of lower tail
        real(dp), intent(out) :: logsf !! logarithm of upper tail
        real(dp) :: lower, upper, logz
        integer :: ik, ilo, ihi
        if (.not. valid_parameters(m,n,draws,odds) .or. ieee_is_nan(k)) then
            cdf=quiet_nan(k); sf=cdf; logcdf=cdf; logsf=cdf; return
        end if
        lower=support_lower(m,n,draws); upper=support_upper(n,draws)
        if (k < lower) then
            cdf=0.0_dp; sf=1.0_dp; logcdf=negative_infinity(k); logsf=0.0_dp; return
        end if
        if (.not. ieee_is_finite(k) .or. k >= upper) then
            cdf=1.0_dp; sf=0.0_dp; logcdf=0.0_dp; logsf=negative_infinity(k); return
        end if
        ik=int(k); ilo=int(lower); ihi=int(upper); logz=log_normalizer(m,n,draws,odds)
        logcdf=logsum_range(ilo,ik,m,n,draws,odds)-logz
        logsf=logsum_range(ik+1,ihi,m,n,draws,odds)-logz
        cdf=exp(logcdf); sf=exp(logsf)
    end subroutine tails

    pure function expected_count(m,n,draws,odds) result(mean_count)
        real(dp), intent(in) :: m !! positive integer population size
        real(dp), intent(in) :: n !! integer number of Type-I objects
        real(dp), intent(in) :: draws !! integer sample size
        real(dp), intent(in) :: odds !! strictly positive Type-I odds ratio
        real(dp) :: mean_count, logz, lt, maxlog, sumw, sumkw
        integer :: j, ilo, ihi
        ilo=int(support_lower(m,n,draws)); ihi=int(support_upper(n,draws)); logz=log_normalizer(m,n,draws,odds)
        maxlog=negative_infinity(1.0_dp)
        do j=ilo,ihi
            maxlog=max(maxlog,log_weight(j,m,n,draws,odds)-logz)
        end do
        sumw=0.0_dp; sumkw=0.0_dp
        do j=ilo,ihi
            lt=exp(log_weight(j,m,n,draws,odds)-logz-maxlog)
            sumw=sumw+lt; sumkw=sumkw+real(j,dp)*lt
        end do
        mean_count=sumkw/sumw
    end function expected_count

    pure function log_normalizer(m,n,draws,odds) result(value)
        real(dp), intent(in) :: m !! positive integer population size
        real(dp), intent(in) :: n !! integer number of Type-I objects
        real(dp), intent(in) :: draws !! integer sample size
        real(dp), intent(in) :: odds !! strictly positive Type-I odds ratio
        real(dp) :: value
        value=logsum_range(int(support_lower(m,n,draws)),int(support_upper(n,draws)),m,n,draws,odds)
    end function log_normalizer

    pure function logsum_range(lo,hi,m,n,draws,odds) result(value)
        integer, intent(in) :: lo !! inclusive lower support index
        integer, intent(in) :: hi !! inclusive upper support index
        real(dp), intent(in) :: m !! positive integer population size
        real(dp), intent(in) :: n !! integer number of Type-I objects
        real(dp), intent(in) :: draws !! integer sample size
        real(dp), intent(in) :: odds !! strictly positive Type-I odds ratio
        real(dp) :: value
        integer :: j
        value=negative_infinity(1.0_dp)
        do j=lo,hi
            value=logadd(value,log_weight(j,m,n,draws,odds))
        end do
    end function logsum_range

    pure elemental function log_weight(j,m,n,draws,odds) result(value)
        integer, intent(in) :: j !! support index
        real(dp), intent(in) :: m !! positive integer population size
        real(dp), intent(in) :: n !! integer number of Type-I objects
        real(dp), intent(in) :: draws !! integer sample size
        real(dp), intent(in) :: odds !! strictly positive Type-I odds ratio
        real(dp) :: value
        value=log_choose(n,real(j,dp))+log_choose(m-n,draws-real(j,dp))+real(j,dp)*log(odds)
    end function log_weight

    pure elemental function log_choose(total,chosen) result(y)
        real(dp), intent(in) :: total !! nonnegative integer upper argument
        real(dp), intent(in) :: chosen !! integer lower argument in [0,total]
        real(dp) :: y
        y=log_gamma(total+1.0_dp)-log_gamma(chosen+1.0_dp)-log_gamma(total-chosen+1.0_dp)
    end function log_choose

    pure elemental function logadd(a,b) result(y)
        real(dp), intent(in) :: a !! first logarithm
        real(dp), intent(in) :: b !! second logarithm
        real(dp) :: y,hi,lo
        hi=max(a,b); lo=min(a,b)
        if (.not. ieee_is_finite(hi)) then; y=hi; else; y=hi+log1p_safe(exp(lo-hi)); end if
    end function logadd

    pure elemental function support_lower(m,n,draws) result(y)
        real(dp), intent(in) :: m !! population size
        real(dp), intent(in) :: n !! Type-I count
        real(dp), intent(in) :: draws !! sample size
        real(dp) :: y
        y=max(0.0_dp,draws-(m-n))
    end function support_lower

    pure elemental function support_upper(n,draws) result(y)
        real(dp), intent(in) :: n !! Type-I count
        real(dp), intent(in) :: draws !! sample size
        real(dp) :: y
        y=min(n,draws)
    end function support_upper

    pure elemental function shifted_count(k,loc) result(count)
        real(dp), intent(in) :: k !! observation
        real(dp), intent(in), optional :: loc !! support shift
        real(dp) :: count
        count=k-optional_loc(loc)
    end function shifted_count

    pure elemental function floor_count(k,loc) result(count)
        real(dp), intent(in) :: k !! CDF/SF argument
        real(dp), intent(in), optional :: loc !! support shift
        real(dp) :: count
        count=k-optional_loc(loc)
        if (ieee_is_finite(count)) count=floor(count)
    end function floor_count

    pure elemental function optional_loc(loc) result(shift)
        real(dp), intent(in), optional :: loc !! support shift
        real(dp) :: shift
        shift=0.0_dp
        if (present(loc)) shift=loc
    end function optional_loc

    pure elemental logical function valid_parameters(m,n,draws,odds)
        real(dp), intent(in) :: m !! population size
        real(dp), intent(in) :: n !! Type-I count
        real(dp), intent(in) :: draws !! sample size
        real(dp), intent(in) :: odds !! Type-I odds ratio
        valid_parameters=ieee_is_finite(m).and.ieee_is_finite(n).and.ieee_is_finite(draws).and. &
            ieee_is_finite(odds).and.m>0.0_dp.and.n>=0.0_dp.and.draws>=0.0_dp.and.n<=m.and.draws<=m.and. &
            m==aint(m).and.n==aint(n).and.draws==aint(draws).and.m<=real(huge(0),dp).and.odds>0.0_dp
    end function valid_parameters

    pure elemental logical function valid_probability(p)
        real(dp), intent(in) :: p !! probability candidate
        valid_probability=ieee_is_finite(p).and.p>=0.0_dp.and.p<=1.0_dp
    end function valid_probability

end module scifort_nchypergeom_fisher
