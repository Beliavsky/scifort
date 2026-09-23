! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Skellam distribution matching scipy.stats.skellam.
module scifort_skellam
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : negative_infinity, positive_infinity, quiet_nan
    use scifort_ncx2, only : ncx2_cdf, ncx2_sf
    implicit none
    private

    public :: skellam_cdf, skellam_isf, skellam_logcdf, skellam_logpmf
    public :: skellam_logsf, skellam_pmf, skellam_ppf, skellam_sf
    public :: skellam_logpmf_derivatives

    interface skellam_pmf
        module procedure skellam_pmf_real
        module procedure skellam_pmf_int
    end interface
    interface skellam_logpmf
        module procedure skellam_logpmf_real
        module procedure skellam_logpmf_int
    end interface
    interface skellam_cdf
        module procedure skellam_cdf_real
        module procedure skellam_cdf_int
    end interface
    interface skellam_sf
        module procedure skellam_sf_real
        module procedure skellam_sf_int
    end interface
    interface skellam_logcdf
        module procedure skellam_logcdf_real
        module procedure skellam_logcdf_int
    end interface
    interface skellam_logsf
        module procedure skellam_logsf_real
        module procedure skellam_logsf_int
    end interface
contains
    pure elemental function skellam_logpmf_real(k, mu1, mu2, loc) result(y)
        real(dp), intent(in) :: k !! lattice observation
        real(dp), intent(in) :: mu1 !! positive first Poisson mean
        real(dp), intent(in) :: mu2 !! positive second Poisson mean
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y, count, x
        integer :: n
        count = k - optional_loc(loc)
        if (.not. valid_shapes(mu1,mu2) .or. ieee_is_nan(count)) then
            y = quiet_nan(k)
        else if (.not. ieee_is_finite(count) .or. count /= aint(count)) then
            y = negative_infinity(k)
        else
            n = int(abs(count))
            x = 2.0_dp * sqrt(mu1*mu2)
            y = -(mu1+mu2) + 0.5_dp*count*log(mu1/mu2) + log_besseli_integer(n,x)
        end if
    end function skellam_logpmf_real
    pure elemental function skellam_logpmf_int(k, mu1, mu2, loc) result(y)
        integer, intent(in) :: k !! integer observation
        real(dp), intent(in) :: mu1 !! positive first Poisson mean
        real(dp), intent(in) :: mu2 !! positive second Poisson mean
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = skellam_logpmf_real(real(k,dp),mu1,mu2,loc)
    end function skellam_logpmf_int
    pure elemental function skellam_pmf_real(k, mu1, mu2, loc) result(y)
        real(dp), intent(in) :: k !! lattice observation
        real(dp), intent(in) :: mu1 !! positive first Poisson mean
        real(dp), intent(in) :: mu2 !! positive second Poisson mean
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = exp(skellam_logpmf_real(k,mu1,mu2,loc))
    end function skellam_pmf_real
    pure elemental function skellam_pmf_int(k, mu1, mu2, loc) result(y)
        integer, intent(in) :: k !! integer observation
        real(dp), intent(in) :: mu1 !! positive first Poisson mean
        real(dp), intent(in) :: mu2 !! positive second Poisson mean
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y = skellam_pmf_real(real(k,dp),mu1,mu2,loc)
    end function skellam_pmf_int
    pure elemental function skellam_cdf_real(k, mu1, mu2, loc) result(y)
        real(dp), intent(in) :: k !! upper limit of P(X <= k); shifted floor is used
        real(dp), intent(in) :: mu1 !! positive first Poisson mean
        real(dp), intent(in) :: mu2 !! positive second Poisson mean
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y, count
        integer :: n
        count = floor(k-optional_loc(loc))
        if (.not. valid_shapes(mu1,mu2) .or. ieee_is_nan(count)) then
            y = quiet_nan(k)
        else if (.not. ieee_is_finite(count)) then
            if (count < 0.0_dp) then; y=0.0_dp; else; y=1.0_dp; end if
        else
            n = int(count)
            if (n >= 0) then
                y = ncx2_sf(2.0_dp*mu1, 2.0_dp*real(n+1,dp), 2.0_dp*mu2)
            else
                y = ncx2_cdf(2.0_dp*mu2, -2.0_dp*real(n,dp), 2.0_dp*mu1)
            end if
        end if
    end function skellam_cdf_real
    pure elemental function skellam_cdf_int(k, mu1, mu2, loc) result(y)
        integer, intent(in) :: k !! integer upper limit
        real(dp), intent(in) :: mu1 !! positive first Poisson mean
        real(dp), intent(in) :: mu2 !! positive second Poisson mean
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y=skellam_cdf_real(real(k,dp),mu1,mu2,loc)
    end function skellam_cdf_int
    pure elemental function skellam_sf_real(k, mu1, mu2, loc) result(y)
        real(dp), intent(in) :: k !! lower limit of P(X > k); shifted floor is used
        real(dp), intent(in) :: mu1 !! positive first Poisson mean
        real(dp), intent(in) :: mu2 !! positive second Poisson mean
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y, count
        integer :: n
        count=floor(k-optional_loc(loc))
        if (.not.valid_shapes(mu1,mu2) .or. ieee_is_nan(count)) then
            y=quiet_nan(k)
        else if (.not.ieee_is_finite(count)) then
            if (count<0.0_dp) then; y=1.0_dp; else; y=0.0_dp; end if
        else
            n=int(count)
            if (n>=0) then
                y=ncx2_cdf(2.0_dp*mu1,2.0_dp*real(n+1,dp),2.0_dp*mu2)
            else
                y=ncx2_sf(2.0_dp*mu2,-2.0_dp*real(n,dp),2.0_dp*mu1)
            end if
        end if
    end function skellam_sf_real
    pure elemental function skellam_sf_int(k, mu1, mu2, loc) result(y)
        integer, intent(in) :: k !! integer lower limit
        real(dp), intent(in) :: mu1 !! positive first Poisson mean
        real(dp), intent(in) :: mu2 !! positive second Poisson mean
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y=skellam_sf_real(real(k,dp),mu1,mu2,loc)
    end function skellam_sf_int
    pure elemental function skellam_logcdf_real(k, mu1, mu2, loc) result(y)
        real(dp), intent(in) :: k !! upper limit of log(P(X <= k))
        real(dp), intent(in) :: mu1 !! positive first Poisson mean
        real(dp), intent(in) :: mu2 !! positive second Poisson mean
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y,q
        q=skellam_cdf_real(k,mu1,mu2,loc)
        if (ieee_is_nan(q)) then; y=q; else if(q<=0.0_dp) then; y=negative_infinity(q); else; y=log(q); end if
    end function skellam_logcdf_real
    pure elemental function skellam_logcdf_int(k, mu1, mu2, loc) result(y)
        integer, intent(in) :: k !! integer upper limit
        real(dp), intent(in) :: mu1 !! positive first Poisson mean
        real(dp), intent(in) :: mu2 !! positive second Poisson mean
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y=skellam_logcdf_real(real(k,dp),mu1,mu2,loc)
    end function skellam_logcdf_int
    pure elemental function skellam_logsf_real(k, mu1, mu2, loc) result(y)
        real(dp), intent(in) :: k !! lower limit of log(P(X > k))
        real(dp), intent(in) :: mu1 !! positive first Poisson mean
        real(dp), intent(in) :: mu2 !! positive second Poisson mean
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y,q
        q=skellam_sf_real(k,mu1,mu2,loc)
        if (ieee_is_nan(q)) then; y=q; else if(q<=0.0_dp) then; y=negative_infinity(q); else; y=log(q); end if
    end function skellam_logsf_real
    pure elemental function skellam_logsf_int(k, mu1, mu2, loc) result(y)
        integer, intent(in) :: k !! integer lower limit
        real(dp), intent(in) :: mu1 !! positive first Poisson mean
        real(dp), intent(in) :: mu2 !! positive second Poisson mean
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y
        y=skellam_logsf_real(real(k,dp),mu1,mu2,loc)
    end function skellam_logsf_int
    pure elemental function skellam_ppf(q, mu1, mu2, loc) result(y)
        real(dp), intent(in) :: q !! lower-tail probability in [0,1]
        real(dp), intent(in) :: mu1 !! positive first Poisson mean
        real(dp), intent(in) :: mu2 !! positive second Poisson mean
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y,shift
        integer :: lo,hi,mid
        shift=optional_loc(loc)
        if (.not.valid_shapes(mu1,mu2) .or. .not.valid_probability(q) .or. .not.ieee_is_finite(shift)) then
            y=quiet_nan(q)
        else if(q<=0.0_dp) then; y=negative_infinity(q)
        else if(q>=1.0_dp) then; y=positive_infinity(q)
        else
            lo=int(floor(mu1-mu2-12.0_dp*sqrt(mu1+mu2)-20.0_dp))
            hi=int(ceiling(mu1-mu2+12.0_dp*sqrt(mu1+mu2)+20.0_dp))
            do while(skellam_cdf_real(shift+real(lo,dp),mu1,mu2,shift)>=q); lo=lo*2-1; end do
            do while(skellam_cdf_real(shift+real(hi,dp),mu1,mu2,shift)<q); hi=hi*2+1; end do
            do while(lo<hi)
                mid=lo+(hi-lo)/2
                if(skellam_cdf_real(shift+real(mid,dp),mu1,mu2,shift)>=q) then; hi=mid; else; lo=mid+1; end if
            end do
            y=shift+real(lo,dp)
        end if
    end function skellam_ppf
    pure elemental function skellam_isf(q, mu1, mu2, loc) result(y)
        real(dp), intent(in) :: q !! upper-tail probability in [0,1]
        real(dp), intent(in) :: mu1 !! positive first Poisson mean
        real(dp), intent(in) :: mu2 !! positive second Poisson mean
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: y,shift
        integer :: lo,hi,mid
        shift=optional_loc(loc)
        if (.not.valid_shapes(mu1,mu2) .or. .not.valid_probability(q) .or. .not.ieee_is_finite(shift)) then
            y=quiet_nan(q)
        else if(q>=1.0_dp) then; y=negative_infinity(q)
        else if(q<=0.0_dp) then; y=positive_infinity(q)
        else
            lo=int(floor(mu1-mu2-12.0_dp*sqrt(mu1+mu2)-20.0_dp))
            hi=int(ceiling(mu1-mu2+12.0_dp*sqrt(mu1+mu2)+20.0_dp))
            do while(skellam_sf_real(shift+real(lo,dp),mu1,mu2,shift)<=q); lo=lo*2-1; end do
            do while(skellam_sf_real(shift+real(hi,dp),mu1,mu2,shift)>q); hi=hi*2+1; end do
            do while(lo<hi)
                mid=lo+(hi-lo)/2
                if(skellam_sf_real(shift+real(mid,dp),mu1,mu2,shift)<=q) then; hi=mid; else; lo=mid+1; end if
            end do
            y=shift+real(lo,dp)
        end if
    end function skellam_isf

    pure subroutine skellam_logpmf_derivatives(k, mu1, mu2, logf, dmu1, dmu2)
        real(dp), intent(in) :: k !! integer-valued unshifted observation
        real(dp), intent(in) :: mu1 !! positive first Poisson mean
        real(dp), intent(in) :: mu2 !! positive second Poisson mean
        real(dp), intent(out) :: logf !! log probability mass
        real(dp), intent(out) :: dmu1 !! derivative with respect to mu1
        real(dp), intent(out) :: dmu2 !! derivative with respect to mu2
        real(dp) :: x, log_i0, log_im, log_ip, dlogi_dx
        integer :: n
        if (.not.valid_shapes(mu1,mu2) .or. k/=aint(k)) then
            logf=quiet_nan(k); dmu1=logf; dmu2=logf; return
        end if
        logf=skellam_logpmf_real(k,mu1,mu2,0.0_dp)
        n=int(abs(k))
        x=2.0_dp*sqrt(mu1*mu2)
        log_i0=log_besseli_integer(n,x)
        if (n == 0) then
            log_ip=log_besseli_integer(1,x)
            dlogi_dx=exp(log_ip-log_i0)
        else
            log_im=log_besseli_integer(n-1,x)
            log_ip=log_besseli_integer(n+1,x)
            dlogi_dx=0.5_dp*(exp(log_im-log_i0)+exp(log_ip-log_i0))
        end if
        dmu1=-1.0_dp+0.5_dp*k/mu1+dlogi_dx*sqrt(mu2/mu1)
        dmu2=-1.0_dp-0.5_dp*k/mu2+dlogi_dx*sqrt(mu1/mu2)
    end subroutine skellam_logpmf_derivatives

    pure elemental function log_besseli_integer(n,x) result(y)
        integer, intent(in) :: n !! nonnegative integer order
        real(dp), intent(in) :: x !! nonnegative real argument
        real(dp) :: y,u,logt,term,sumv
        integer :: j0,j
        if(n<0 .or. x<0.0_dp .or. ieee_is_nan(x)) then; y=quiet_nan(x); return; end if
        if(x==0.0_dp) then
            if(n==0) then; y=0.0_dp; else; y=negative_infinity(x); end if
            return
        end if
        u=0.25_dp*x*x
        j0=max(0,int(floor(0.5_dp*(sqrt(real(n*n,dp)+x*x)-real(n,dp)-2.0_dp))))
        logt=(2.0_dp*real(j0,dp)+real(n,dp))*log(0.5_dp*x)-log_gamma(real(j0+1,dp))-log_gamma(real(j0+n+1,dp))
        sumv=1.0_dp
        term=1.0_dp
        do j=j0-1,0,-1
            term=term*real((j+1)*(j+n+1),dp)/u
            sumv=sumv+term
            if(term<epsilon(sumv)*sumv) exit
        end do
        term=1.0_dp
        do j=j0+1,j0+100000
            term=term*u/(real(j,dp)*real(j+n,dp))
            sumv=sumv+term
            if(term<epsilon(sumv)*sumv .and. j>j0+8) exit
        end do
        y=logt+log(sumv)
    end function log_besseli_integer

    pure elemental function optional_loc(loc) result(value)
        real(dp), intent(in), optional :: loc !! support shift
        real(dp) :: value
        value=0.0_dp; if(present(loc)) value=loc
    end function optional_loc
    pure elemental logical function valid_shapes(mu1,mu2)
        real(dp), intent(in) :: mu1 !! first Poisson mean
        real(dp), intent(in) :: mu2 !! second Poisson mean
        valid_shapes=ieee_is_finite(mu1).and.ieee_is_finite(mu2).and.mu1>0.0_dp.and.mu2>0.0_dp
    end function valid_shapes
    pure elemental logical function valid_probability(q)
        real(dp), intent(in) :: q !! probability to validate
        valid_probability=q>=0.0_dp.and.q<=1.0_dp
    end function valid_probability
end module scifort_skellam
