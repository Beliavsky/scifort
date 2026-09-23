! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Wallenius' noncentral hypergeometric distribution matching
! scipy.stats.nchypergeom_wallenius.  The PMF is computed from the defining
! sequential biased-urn process by a finite-state dynamic program.
module scifort_nchypergeom_wallenius
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : log1p_safe, negative_infinity, quiet_nan
    implicit none
    private

    public :: nchypergeom_wallenius_cdf, nchypergeom_wallenius_isf
    public :: nchypergeom_wallenius_logcdf, nchypergeom_wallenius_logpmf
    public :: nchypergeom_wallenius_logsf, nchypergeom_wallenius_pmf
    public :: nchypergeom_wallenius_ppf, nchypergeom_wallenius_sf
    public :: nchypergeom_wallenius_logpmf_derivative_odds

    interface nchypergeom_wallenius_pmf
        module procedure nchypergeom_wallenius_pmf_real
        module procedure nchypergeom_wallenius_pmf_int
    end interface
    interface nchypergeom_wallenius_logpmf
        module procedure nchypergeom_wallenius_logpmf_real
        module procedure nchypergeom_wallenius_logpmf_int
    end interface
    interface nchypergeom_wallenius_cdf
        module procedure nchypergeom_wallenius_cdf_real
        module procedure nchypergeom_wallenius_cdf_int
    end interface
    interface nchypergeom_wallenius_sf
        module procedure nchypergeom_wallenius_sf_real
        module procedure nchypergeom_wallenius_sf_int
    end interface
    interface nchypergeom_wallenius_logcdf
        module procedure nchypergeom_wallenius_logcdf_real
        module procedure nchypergeom_wallenius_logcdf_int
    end interface
    interface nchypergeom_wallenius_logsf
        module procedure nchypergeom_wallenius_logsf_real
        module procedure nchypergeom_wallenius_logsf_int
    end interface

contains

    pure function nchypergeom_wallenius_logpmf_real(k, m, n, draws, odds, loc) result(y)
        real(dp), intent(in) :: k !! observed Type-I count
        real(dp), intent(in) :: m !! positive integer population size
        real(dp), intent(in) :: n !! integer number of Type-I objects in [0,m]
        real(dp), intent(in) :: draws !! integer sample size in [0,m]
        real(dp), intent(in) :: odds !! strictly positive Type-I odds ratio
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y, count
        real(dp), allocatable :: logp(:)
        integer :: lo, hi, ik
        count = shifted_count(k, loc)
        if (.not. valid_parameters(m,n,draws,odds) .or. ieee_is_nan(count)) then
            y = quiet_nan(k); return
        end if
        lo = int(support_lower(m,n,draws)); hi = int(support_upper(n,draws))
        if (.not. ieee_is_finite(count) .or. abs(count-aint(count)) > 0.0_dp .or. count < real(lo,dp) .or. count > real(hi,dp)) then
            y = negative_infinity(k); return
        end if
        call build_logpmf(int(m),int(n),int(draws),odds,logp)
        ik = int(count)
        y = logp(ik)
    end function nchypergeom_wallenius_logpmf_real

    pure function nchypergeom_wallenius_logpmf_int(k, m, n, draws, odds, loc) result(y)
        integer, intent(in) :: k !! observed Type-I count
        integer, intent(in) :: m !! positive population size
        integer, intent(in) :: n !! number of Type-I objects
        integer, intent(in) :: draws !! sample size
        real(dp), intent(in) :: odds !! strictly positive Type-I odds ratio
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = nchypergeom_wallenius_logpmf_real(real(k,dp),real(m,dp),real(n,dp),real(draws,dp),odds,loc)
    end function nchypergeom_wallenius_logpmf_int

    pure function nchypergeom_wallenius_pmf_real(k, m, n, draws, odds, loc) result(y)
        real(dp), intent(in) :: k !! observed Type-I count
        real(dp), intent(in) :: m !! positive integer population size
        real(dp), intent(in) :: n !! integer number of Type-I objects
        real(dp), intent(in) :: draws !! integer sample size
        real(dp), intent(in) :: odds !! strictly positive Type-I odds ratio
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = exp(nchypergeom_wallenius_logpmf_real(k,m,n,draws,odds,loc))
    end function nchypergeom_wallenius_pmf_real

    pure function nchypergeom_wallenius_pmf_int(k, m, n, draws, odds, loc) result(y)
        integer, intent(in) :: k !! observed Type-I count
        integer, intent(in) :: m !! positive population size
        integer, intent(in) :: n !! number of Type-I objects
        integer, intent(in) :: draws !! sample size
        real(dp), intent(in) :: odds !! strictly positive Type-I odds ratio
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = nchypergeom_wallenius_pmf_real(real(k,dp),real(m,dp),real(n,dp),real(draws,dp),odds,loc)
    end function nchypergeom_wallenius_pmf_int

    pure function nchypergeom_wallenius_cdf_real(k, m, n, draws, odds, loc) result(y)
        real(dp), intent(in) :: k !! upper limit of P(X <= k)
        real(dp), intent(in) :: m !! positive integer population size
        real(dp), intent(in) :: n !! integer number of Type-I objects
        real(dp), intent(in) :: draws !! integer sample size
        real(dp), intent(in) :: odds !! strictly positive Type-I odds ratio
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y, sf, logcdf, logsf
        call tails(floor_count(k,loc),m,n,draws,odds,y,sf,logcdf,logsf)
    end function nchypergeom_wallenius_cdf_real

    pure function nchypergeom_wallenius_cdf_int(k, m, n, draws, odds, loc) result(y)
        integer, intent(in) :: k !! upper integer limit
        integer, intent(in) :: m !! positive population size
        integer, intent(in) :: n !! number of Type-I objects
        integer, intent(in) :: draws !! sample size
        real(dp), intent(in) :: odds !! strictly positive Type-I odds ratio
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = nchypergeom_wallenius_cdf_real(real(k,dp),real(m,dp),real(n,dp),real(draws,dp),odds,loc)
    end function nchypergeom_wallenius_cdf_int

    pure function nchypergeom_wallenius_sf_real(k, m, n, draws, odds, loc) result(y)
        real(dp), intent(in) :: k !! lower limit of P(X > k)
        real(dp), intent(in) :: m !! positive integer population size
        real(dp), intent(in) :: n !! integer number of Type-I objects
        real(dp), intent(in) :: draws !! integer sample size
        real(dp), intent(in) :: odds !! strictly positive Type-I odds ratio
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y, cdf, logcdf, logsf
        call tails(floor_count(k,loc),m,n,draws,odds,cdf,y,logcdf,logsf)
    end function nchypergeom_wallenius_sf_real

    pure function nchypergeom_wallenius_sf_int(k, m, n, draws, odds, loc) result(y)
        integer, intent(in) :: k !! lower integer limit
        integer, intent(in) :: m !! positive population size
        integer, intent(in) :: n !! number of Type-I objects
        integer, intent(in) :: draws !! sample size
        real(dp), intent(in) :: odds !! strictly positive Type-I odds ratio
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = nchypergeom_wallenius_sf_real(real(k,dp),real(m,dp),real(n,dp),real(draws,dp),odds,loc)
    end function nchypergeom_wallenius_sf_int

    pure function nchypergeom_wallenius_logcdf_real(k, m, n, draws, odds, loc) result(y)
        real(dp), intent(in) :: k !! upper limit of log(P(X <= k))
        real(dp), intent(in) :: m !! positive integer population size
        real(dp), intent(in) :: n !! integer number of Type-I objects
        real(dp), intent(in) :: draws !! integer sample size
        real(dp), intent(in) :: odds !! strictly positive Type-I odds ratio
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y, cdf, sf, logsf
        call tails(floor_count(k,loc),m,n,draws,odds,cdf,sf,y,logsf)
    end function nchypergeom_wallenius_logcdf_real

    pure function nchypergeom_wallenius_logcdf_int(k, m, n, draws, odds, loc) result(y)
        integer, intent(in) :: k !! upper integer limit
        integer, intent(in) :: m !! positive population size
        integer, intent(in) :: n !! number of Type-I objects
        integer, intent(in) :: draws !! sample size
        real(dp), intent(in) :: odds !! strictly positive Type-I odds ratio
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = nchypergeom_wallenius_logcdf_real(real(k,dp),real(m,dp),real(n,dp),real(draws,dp),odds,loc)
    end function nchypergeom_wallenius_logcdf_int

    pure function nchypergeom_wallenius_logsf_real(k, m, n, draws, odds, loc) result(y)
        real(dp), intent(in) :: k !! lower limit of log(P(X > k))
        real(dp), intent(in) :: m !! positive integer population size
        real(dp), intent(in) :: n !! integer number of Type-I objects
        real(dp), intent(in) :: draws !! integer sample size
        real(dp), intent(in) :: odds !! strictly positive Type-I odds ratio
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y, cdf, sf, logcdf
        call tails(floor_count(k,loc),m,n,draws,odds,cdf,sf,logcdf,y)
    end function nchypergeom_wallenius_logsf_real

    pure function nchypergeom_wallenius_logsf_int(k, m, n, draws, odds, loc) result(y)
        integer, intent(in) :: k !! lower integer limit
        integer, intent(in) :: m !! positive population size
        integer, intent(in) :: n !! number of Type-I objects
        integer, intent(in) :: draws !! sample size
        real(dp), intent(in) :: odds !! strictly positive Type-I odds ratio
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = nchypergeom_wallenius_logsf_real(real(k,dp),real(m,dp),real(n,dp),real(draws,dp),odds,loc)
    end function nchypergeom_wallenius_logsf_int

    pure function nchypergeom_wallenius_ppf(probability, m, n, draws, odds, loc) result(y)
        real(dp), intent(in) :: probability !! lower-tail probability in [0,1]
        real(dp), intent(in) :: m !! positive integer population size
        real(dp), intent(in) :: n !! integer number of Type-I objects
        real(dp), intent(in) :: draws !! integer sample size
        real(dp), intent(in) :: odds !! strictly positive Type-I odds ratio
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y, shift
        integer :: lo, hi, mid
        shift = optional_loc(loc)
        if (.not. valid_parameters(m,n,draws,odds) .or. .not. ieee_is_finite(shift) .or. &
            .not. valid_probability(probability)) then
            y = quiet_nan(probability); return
        end if
        lo=int(support_lower(m,n,draws)); hi=int(support_upper(n,draws))
        if (probability <= 0.0_dp) then
            y=shift+real(lo-1,dp); return
        else if (probability >= 1.0_dp) then
            y=shift+real(hi,dp); return
        end if
        do while (lo < hi)
            mid=lo+(hi-lo)/2
            if (nchypergeom_wallenius_cdf_real(shift+real(mid,dp),m,n,draws,odds,shift) >= probability) then
                hi=mid
            else
                lo=mid+1
            end if
        end do
        y=shift+real(lo,dp)
    end function nchypergeom_wallenius_ppf

    pure function nchypergeom_wallenius_isf(probability, m, n, draws, odds, loc) result(y)
        real(dp), intent(in) :: probability !! upper-tail probability in [0,1]
        real(dp), intent(in) :: m !! positive integer population size
        real(dp), intent(in) :: n !! integer number of Type-I objects
        real(dp), intent(in) :: draws !! integer sample size
        real(dp), intent(in) :: odds !! strictly positive Type-I odds ratio
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y, shift
        integer :: lo, hi, mid
        shift = optional_loc(loc)
        if (.not. valid_parameters(m,n,draws,odds) .or. .not. ieee_is_finite(shift) .or. &
            .not. valid_probability(probability)) then
            y = quiet_nan(probability); return
        end if
        lo=int(support_lower(m,n,draws)); hi=int(support_upper(n,draws))
        if (probability >= 1.0_dp) then
            y=shift+real(lo-1,dp); return
        else if (probability <= 0.0_dp) then
            y=shift+real(hi,dp); return
        end if
        do while (lo < hi)
            mid=lo+(hi-lo)/2
            if (nchypergeom_wallenius_sf_real(shift+real(mid,dp),m,n,draws,odds,shift) <= probability) then
                hi=mid
            else
                lo=mid+1
            end if
        end do
        y=shift+real(lo,dp)
    end function nchypergeom_wallenius_isf

    pure subroutine nchypergeom_wallenius_logpmf_derivative_odds(count,m,n,draws,odds,logf,dodds)
        real(dp), intent(in) :: count !! unshifted integer observation
        real(dp), intent(in) :: m !! positive integer population size
        real(dp), intent(in) :: n !! integer number of Type-I objects
        real(dp), intent(in) :: draws !! integer sample size
        real(dp), intent(in) :: odds !! strictly positive Type-I odds ratio
        real(dp), intent(out) :: logf !! log probability mass
        real(dp), intent(out) :: dodds !! derivative of log mass with respect to odds
        real(dp), allocatable :: p(:), dpdw(:)
        integer :: ik, lo, hi
        if (.not. valid_parameters(m,n,draws,odds) .or. .not. ieee_is_finite(count) .or. abs(count-aint(count)) > 0.0_dp) then
            logf=quiet_nan(count); dodds=logf; return
        end if
        lo=int(support_lower(m,n,draws)); hi=int(support_upper(n,draws)); ik=int(count)
        if (ik < lo .or. ik > hi) then
            logf=negative_infinity(count); dodds=quiet_nan(count); return
        end if
        call build_pmf_derivative(int(m),int(n),int(draws),odds,p,dpdw)
        if (p(ik) <= 0.0_dp) then
            logf=negative_infinity(count); dodds=quiet_nan(count)
        else
            logf=log(p(ik)); dodds=dpdw(ik)/p(ik)
        end if
    end subroutine nchypergeom_wallenius_logpmf_derivative_odds

    pure subroutine tails(k,m,n,draws,odds,cdf,sf,logcdf,logsf)
        real(dp), intent(in) :: k !! unshifted integer tail split
        real(dp), intent(in) :: m !! positive integer population size
        real(dp), intent(in) :: n !! integer number of Type-I objects
        real(dp), intent(in) :: draws !! integer sample size
        real(dp), intent(in) :: odds !! strictly positive Type-I odds ratio
        real(dp), intent(out) :: cdf !! lower-tail probability
        real(dp), intent(out) :: sf !! upper-tail probability
        real(dp), intent(out) :: logcdf !! logarithm of lower tail
        real(dp), intent(out) :: logsf !! logarithm of upper tail
        real(dp), allocatable :: logp(:)
        real(dp) :: lc, ls
        integer :: lo, hi, ik, j
        if (.not. valid_parameters(m,n,draws,odds) .or. ieee_is_nan(k)) then
            cdf=quiet_nan(k); sf=cdf; logcdf=cdf; logsf=cdf; return
        end if
        lo=int(support_lower(m,n,draws)); hi=int(support_upper(n,draws))
        if (k < real(lo,dp)) then
            cdf=0.0_dp; sf=1.0_dp; logcdf=negative_infinity(k); logsf=0.0_dp; return
        end if
        if (.not. ieee_is_finite(k) .or. k >= real(hi,dp)) then
            cdf=1.0_dp; sf=0.0_dp; logcdf=0.0_dp; logsf=negative_infinity(k); return
        end if
        call build_logpmf(int(m),int(n),int(draws),odds,logp)
        ik=int(k); lc=negative_infinity(1.0_dp); ls=negative_infinity(1.0_dp)
        do j=lo,ik; lc=logadd(lc,logp(j)); end do
        do j=ik+1,hi; ls=logadd(ls,logp(j)); end do
        logcdf=lc; logsf=ls; cdf=exp(lc); sf=exp(ls)
    end subroutine tails

    pure subroutine build_logpmf(m,n,draws,odds,logp)
        integer, intent(in) :: m !! population size
        integer, intent(in) :: n !! Type-I population count
        integer, intent(in) :: draws !! number of sequential draws
        real(dp), intent(in) :: odds !! positive Type-I odds ratio
        real(dp), allocatable, intent(out) :: logp(:) !! log PMF indexed by Type-I count
        real(dp), allocatable :: cur(:), nxt(:)
        real(dp) :: den, loga, logb
        integer :: j, x, xmin, xmax, aleft, bleft
        allocate(cur(0:draws),nxt(0:draws),logp(0:draws))
        cur=negative_infinity(1.0_dp); cur(0)=0.0_dp
        do j=0,draws-1
            nxt=negative_infinity(1.0_dp)
            xmin=max(0,j-(m-n)); xmax=min(j,n)
            do x=xmin,xmax
                if (.not. ieee_is_finite(cur(x))) cycle
                aleft=n-x; bleft=(m-n)-(j-x); den=odds*real(aleft,dp)+real(bleft,dp)
                if (aleft > 0) then
                    loga=log(odds*real(aleft,dp)/den)
                    nxt(x+1)=logadd(nxt(x+1),cur(x)+loga)
                end if
                if (bleft > 0) then
                    logb=log(real(bleft,dp)/den)
                    nxt(x)=logadd(nxt(x),cur(x)+logb)
                end if
            end do
            cur=nxt
        end do
        logp=cur
    end subroutine build_logpmf

    pure subroutine build_pmf_derivative(m,n,draws,odds,p,dpdw)
        integer, intent(in) :: m !! population size
        integer, intent(in) :: n !! Type-I population count
        integer, intent(in) :: draws !! number of sequential draws
        real(dp), intent(in) :: odds !! positive Type-I odds ratio
        real(dp), allocatable, intent(out) :: p(:) !! PMF indexed by Type-I count
        real(dp), allocatable, intent(out) :: dpdw(:) !! odds derivative of PMF
        real(dp), allocatable :: cur(:), nxt(:), dcur(:), dnxt(:)
        real(dp) :: den, q, dq
        integer :: j, x, xmin, xmax, aleft, bleft
        allocate(cur(0:draws),nxt(0:draws),dcur(0:draws),dnxt(0:draws),p(0:draws),dpdw(0:draws))
        cur=0.0_dp; dcur=0.0_dp; cur(0)=1.0_dp
        do j=0,draws-1
            nxt=0.0_dp; dnxt=0.0_dp
            xmin=max(0,j-(m-n)); xmax=min(j,n)
            do x=xmin,xmax
                aleft=n-x; bleft=(m-n)-(j-x); den=odds*real(aleft,dp)+real(bleft,dp)
                if (aleft > 0) then
                    q=odds*real(aleft,dp)/den
                    dq=real(aleft,dp)*real(bleft,dp)/(den*den)
                    nxt(x+1)=nxt(x+1)+cur(x)*q
                    dnxt(x+1)=dnxt(x+1)+dcur(x)*q+cur(x)*dq
                end if
                if (bleft > 0) then
                    if (aleft > 0) then
                        q=odds*real(aleft,dp)/den
                        dq=real(aleft,dp)*real(bleft,dp)/(den*den)
                    else
                        q=0.0_dp; dq=0.0_dp
                    end if
                    nxt(x)=nxt(x)+cur(x)*(1.0_dp-q)
                    dnxt(x)=dnxt(x)+dcur(x)*(1.0_dp-q)-cur(x)*dq
                end if
            end do
            cur=nxt; dcur=dnxt
        end do
        p=cur; dpdw=dcur
    end subroutine build_pmf_derivative

    pure elemental function logadd(a,b) result(y)
        real(dp), intent(in) :: a !! first logarithm
        real(dp), intent(in) :: b !! second logarithm
        real(dp) :: y, hi, lo
        hi=max(a,b); lo=min(a,b)
        if (.not. ieee_is_finite(hi)) then
            y=hi
        else
            y=hi+log1p_safe(exp(lo-hi))
        end if
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
            abs(m-aint(m))<=0.0_dp.and.abs(n-aint(n))<=0.0_dp.and. &
            abs(draws-aint(draws))<=0.0_dp.and.m<=real(huge(0),dp).and.odds>0.0_dp
    end function valid_parameters

    pure elemental logical function valid_probability(p)
        real(dp), intent(in) :: p !! probability candidate
        valid_probability=ieee_is_finite(p).and.p>=0.0_dp.and.p<=1.0_dp
    end function valid_probability

end module scifort_nchypergeom_wallenius
