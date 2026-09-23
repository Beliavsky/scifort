! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Poisson-binomial distribution matching scipy.stats.poisson_binom.
! Probabilities are supplied as a rank-one shape vector p(:).
module scifort_poisson_binom
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : log1p_safe, negative_infinity, quiet_nan
    implicit none
    private

    public :: poisson_binom_cdf, poisson_binom_isf, poisson_binom_logcdf
    public :: poisson_binom_logpmf, poisson_binom_logsf, poisson_binom_pmf
    public :: poisson_binom_ppf, poisson_binom_sf
    public :: poisson_binom_logpmf_derivative_p

    interface poisson_binom_pmf
        module procedure poisson_binom_pmf_real
        module procedure poisson_binom_pmf_int
    end interface
    interface poisson_binom_logpmf
        module procedure poisson_binom_logpmf_real
        module procedure poisson_binom_logpmf_int
    end interface
    interface poisson_binom_cdf
        module procedure poisson_binom_cdf_real
        module procedure poisson_binom_cdf_int
    end interface
    interface poisson_binom_sf
        module procedure poisson_binom_sf_real
        module procedure poisson_binom_sf_int
    end interface
    interface poisson_binom_logcdf
        module procedure poisson_binom_logcdf_real
        module procedure poisson_binom_logcdf_int
    end interface
    interface poisson_binom_logsf
        module procedure poisson_binom_logsf_real
        module procedure poisson_binom_logsf_int
    end interface

contains

    pure function poisson_binom_logpmf_real(k,p,loc) result(y)
        real(dp), intent(in) :: k !! observed success count
        real(dp), intent(in) :: p(:) !! Bernoulli success probabilities in [0,1]
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y, count
        real(dp), allocatable :: logp(:)
        integer :: ik
        count=shifted_count(k,loc)
        if (.not. valid_parameters(p) .or. ieee_is_nan(count)) then
            y=quiet_nan(k); return
        end if
        if (.not. ieee_is_finite(count) .or. abs(count-aint(count))>0.0_dp .or. &
            count<0.0_dp .or. count>real(size(p),dp)) then
            y=negative_infinity(k); return
        end if
        call build_logpmf(p,logp); ik=int(count); y=logp(ik)
    end function poisson_binom_logpmf_real

    pure function poisson_binom_logpmf_int(k,p,loc) result(y)
        integer, intent(in) :: k !! observed success count
        real(dp), intent(in) :: p(:) !! Bernoulli success probabilities in [0,1]
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y=poisson_binom_logpmf_real(real(k,dp),p,loc)
    end function poisson_binom_logpmf_int

    pure function poisson_binom_pmf_real(k,p,loc) result(y)
        real(dp), intent(in) :: k !! observed success count
        real(dp), intent(in) :: p(:) !! Bernoulli success probabilities in [0,1]
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y=exp(poisson_binom_logpmf_real(k,p,loc))
    end function poisson_binom_pmf_real

    pure function poisson_binom_pmf_int(k,p,loc) result(y)
        integer, intent(in) :: k !! observed success count
        real(dp), intent(in) :: p(:) !! Bernoulli success probabilities in [0,1]
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y=poisson_binom_pmf_real(real(k,dp),p,loc)
    end function poisson_binom_pmf_int

    pure function poisson_binom_cdf_real(k,p,loc) result(y)
        real(dp), intent(in) :: k !! upper limit of P(X <= k)
        real(dp), intent(in) :: p(:) !! Bernoulli success probabilities in [0,1]
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y, sf, logcdf, logsf
        call tails(floor_count(k,loc),p,y,sf,logcdf,logsf)
    end function poisson_binom_cdf_real

    pure function poisson_binom_cdf_int(k,p,loc) result(y)
        integer, intent(in) :: k !! upper integer limit
        real(dp), intent(in) :: p(:) !! Bernoulli success probabilities in [0,1]
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y=poisson_binom_cdf_real(real(k,dp),p,loc)
    end function poisson_binom_cdf_int

    pure function poisson_binom_sf_real(k,p,loc) result(y)
        real(dp), intent(in) :: k !! lower limit of P(X > k)
        real(dp), intent(in) :: p(:) !! Bernoulli success probabilities in [0,1]
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y, cdf, logcdf, logsf
        call tails(floor_count(k,loc),p,cdf,y,logcdf,logsf)
    end function poisson_binom_sf_real

    pure function poisson_binom_sf_int(k,p,loc) result(y)
        integer, intent(in) :: k !! lower integer limit
        real(dp), intent(in) :: p(:) !! Bernoulli success probabilities in [0,1]
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y=poisson_binom_sf_real(real(k,dp),p,loc)
    end function poisson_binom_sf_int

    pure function poisson_binom_logcdf_real(k,p,loc) result(y)
        real(dp), intent(in) :: k !! upper limit of log(P(X <= k))
        real(dp), intent(in) :: p(:) !! Bernoulli success probabilities in [0,1]
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y, cdf, sf, logsf
        call tails(floor_count(k,loc),p,cdf,sf,y,logsf)
    end function poisson_binom_logcdf_real

    pure function poisson_binom_logcdf_int(k,p,loc) result(y)
        integer, intent(in) :: k !! upper integer limit
        real(dp), intent(in) :: p(:) !! Bernoulli success probabilities in [0,1]
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y=poisson_binom_logcdf_real(real(k,dp),p,loc)
    end function poisson_binom_logcdf_int

    pure function poisson_binom_logsf_real(k,p,loc) result(y)
        real(dp), intent(in) :: k !! lower limit of log(P(X > k))
        real(dp), intent(in) :: p(:) !! Bernoulli success probabilities in [0,1]
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y, cdf, sf, logcdf
        call tails(floor_count(k,loc),p,cdf,sf,logcdf,y)
    end function poisson_binom_logsf_real

    pure function poisson_binom_logsf_int(k,p,loc) result(y)
        integer, intent(in) :: k !! lower integer limit
        real(dp), intent(in) :: p(:) !! Bernoulli success probabilities in [0,1]
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y=poisson_binom_logsf_real(real(k,dp),p,loc)
    end function poisson_binom_logsf_int

    pure function poisson_binom_ppf(probability,p,loc) result(y)
        real(dp), intent(in) :: probability !! lower-tail probability in [0,1]
        real(dp), intent(in) :: p(:) !! Bernoulli success probabilities in [0,1]
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y, shift
        integer :: lo, hi, mid
        shift=optional_loc(loc)
        if (.not.valid_parameters(p) .or. .not.ieee_is_finite(shift) .or. .not.valid_probability(probability)) then
            y=quiet_nan(probability); return
        end if
        lo=0; hi=size(p)
        if (probability<=0.0_dp) then; y=shift-1.0_dp; return; end if
        if (probability>=1.0_dp) then; y=shift+real(hi,dp); return; end if
        do while (lo<hi)
            mid=lo+(hi-lo)/2
            if (poisson_binom_cdf_real(shift+real(mid,dp),p,shift)>=probability) then
                hi=mid
            else
                lo=mid+1
            end if
        end do
        y=shift+real(lo,dp)
    end function poisson_binom_ppf

    pure function poisson_binom_isf(probability,p,loc) result(y)
        real(dp), intent(in) :: probability !! upper-tail probability in [0,1]
        real(dp), intent(in) :: p(:) !! Bernoulli success probabilities in [0,1]
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y, shift
        integer :: lo, hi, mid
        shift=optional_loc(loc)
        if (.not.valid_parameters(p) .or. .not.ieee_is_finite(shift) .or. .not.valid_probability(probability)) then
            y=quiet_nan(probability); return
        end if
        lo=0; hi=size(p)
        if (probability>=1.0_dp) then; y=shift-1.0_dp; return; end if
        if (probability<=0.0_dp) then; y=shift+real(hi,dp); return; end if
        do while (lo<hi)
            mid=lo+(hi-lo)/2
            if (poisson_binom_sf_real(shift+real(mid,dp),p,shift)<=probability) then
                hi=mid
            else
                lo=mid+1
            end if
        end do
        y=shift+real(lo,dp)
    end function poisson_binom_isf

    pure subroutine poisson_binom_logpmf_derivative_p(count,p,logf,dpvec)
        real(dp), intent(in) :: count !! unshifted integer observation
        real(dp), intent(in) :: p(:) !! Bernoulli success probabilities in [0,1]
        real(dp), intent(out) :: logf !! log probability mass
        real(dp), intent(out) :: dpvec(:) !! derivatives of log mass with respect to p(:)
        real(dp), allocatable :: q(:), mass(:)
        real(dp) :: pk, deriv
        integer :: j, k, n
        n=size(p); dpvec=quiet_nan(count)
        if (size(dpvec)/=n .or. .not.valid_parameters(p) .or. .not.ieee_is_finite(count) .or. &
            abs(count-aint(count))>0.0_dp) then
            logf=quiet_nan(count); return
        end if
        k=int(count)
        if (k<0 .or. k>n) then; logf=negative_infinity(count); return; end if
        logf=poisson_binom_logpmf_real(count,p)
        if (.not.ieee_is_finite(logf)) return
        pk=exp(logf)
        allocate(q(max(0,n-1)))
        do j=1,n
            if (n==1) then
                if (k==0) then; deriv=-1.0_dp; else; deriv=1.0_dp; end if
            else
                if (j>1) q(1:j-1)=p(1:j-1)
                if (j<n) q(j:n-1)=p(j+1:n)
                call build_pmf(q,mass)
                deriv=0.0_dp
                if (k>=1) deriv=deriv+mass(k-1)
                if (k<=n-1) deriv=deriv-mass(k)
            end if
            dpvec(j)=deriv/pk
        end do
    end subroutine poisson_binom_logpmf_derivative_p

    pure subroutine tails(k,p,cdf,sf,logcdf,logsf)
        real(dp), intent(in) :: k !! unshifted integer tail split
        real(dp), intent(in) :: p(:) !! Bernoulli success probabilities in [0,1]
        real(dp), intent(out) :: cdf !! lower-tail probability
        real(dp), intent(out) :: sf !! upper-tail probability
        real(dp), intent(out) :: logcdf !! logarithm of lower tail
        real(dp), intent(out) :: logsf !! logarithm of upper tail
        real(dp), allocatable :: logp(:)
        real(dp) :: lc, ls
        integer :: ik, j, n
        n=size(p)
        if (.not.valid_parameters(p) .or. ieee_is_nan(k)) then
            cdf=quiet_nan(k); sf=cdf; logcdf=cdf; logsf=cdf; return
        end if
        if (k<0.0_dp) then
            cdf=0.0_dp; sf=1.0_dp; logcdf=negative_infinity(k); logsf=0.0_dp; return
        end if
        if (.not.ieee_is_finite(k) .or. k>=real(n,dp)) then
            cdf=1.0_dp; sf=0.0_dp; logcdf=0.0_dp; logsf=negative_infinity(k); return
        end if
        call build_logpmf(p,logp); ik=int(k)
        lc=negative_infinity(1.0_dp); ls=negative_infinity(1.0_dp)
        do j=0,ik; lc=logadd(lc,logp(j)); end do
        do j=ik+1,n; ls=logadd(ls,logp(j)); end do
        logcdf=lc; logsf=ls; cdf=exp(lc); sf=exp(ls)
    end subroutine tails

    pure subroutine build_logpmf(p,logp)
        real(dp), intent(in) :: p(:) !! Bernoulli success probabilities
        real(dp), allocatable, intent(out) :: logp(:) !! log PMF indexed 0:size(p)
        real(dp), allocatable :: cur(:), nxt(:)
        real(dp) :: lp, lq
        integer :: i, k, n
        n=size(p); allocate(cur(0:n),nxt(0:n),logp(0:n))
        cur=negative_infinity(1.0_dp); cur(0)=0.0_dp
        do i=1,n
            nxt=negative_infinity(1.0_dp)
            if (p(i)>0.0_dp) lp=log(p(i))
            if (p(i)<1.0_dp) lq=log1p_safe(-p(i))
            do k=0,i-1
                if (.not.ieee_is_finite(cur(k))) cycle
                if (p(i)<1.0_dp) nxt(k)=logadd(nxt(k),cur(k)+lq)
                if (p(i)>0.0_dp) nxt(k+1)=logadd(nxt(k+1),cur(k)+lp)
            end do
            cur=nxt
        end do
        logp=cur
    end subroutine build_logpmf

    pure subroutine build_pmf(p,mass)
        real(dp), intent(in) :: p(:) !! Bernoulli success probabilities
        real(dp), allocatable, intent(out) :: mass(:) !! PMF indexed 0:size(p)
        real(dp), allocatable :: cur(:), nxt(:)
        integer :: i, k, n
        n=size(p); allocate(cur(0:n),nxt(0:n),mass(0:n)); cur=0.0_dp; cur(0)=1.0_dp
        do i=1,n
            nxt=0.0_dp
            do k=0,i-1
                nxt(k)=nxt(k)+cur(k)*(1.0_dp-p(i))
                nxt(k+1)=nxt(k+1)+cur(k)*p(i)
            end do
            cur=nxt
        end do
        mass=cur
    end subroutine build_pmf

    pure elemental function logadd(a,b) result(y)
        real(dp), intent(in) :: a !! first logarithm
        real(dp), intent(in) :: b !! second logarithm
        real(dp) :: y, hi, lo
        hi=max(a,b); lo=min(a,b)
        if (.not.ieee_is_finite(hi)) then; y=hi; else; y=hi+log1p_safe(exp(lo-hi)); end if
    end function logadd

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
        shift=0.0_dp; if (present(loc)) shift=loc
    end function optional_loc

    pure logical function valid_parameters(p)
        real(dp), intent(in) :: p(:) !! probability vector
        integer :: i
        valid_parameters=size(p)>=1
        if (.not.valid_parameters) return
        do i=1,size(p)
            if (.not.ieee_is_finite(p(i)) .or. p(i)<0.0_dp .or. p(i)>1.0_dp) then
                valid_parameters=.false.; return
            end if
        end do
    end function valid_parameters

    pure elemental logical function valid_probability(q)
        real(dp), intent(in) :: q !! probability candidate
        valid_probability=ieee_is_finite(q).and.q>=0.0_dp.and.q<=1.0_dp
    end function valid_probability

end module scifort_poisson_binom
