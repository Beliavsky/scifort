! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
! Independent constrained derivative-jump smoothing; see CODE_PROVENANCE.md.
module scifort_interpolate_fixed_smoothing
    use scifort_kinds, only: dp
    use scifort_interpolate_bspline, only: bspline,bspline_array_interpolator_type
    use scifort_linalg_svd, only: svd
    use scifort_linalg_qr, only: lstsq
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: fixed_knots_valid,fixed_smoothing
contains
    logical function fixed_knots_valid(t,k,m,bounds,periodic) result(ok)
        implicit none
        real(dp), intent(in) :: t(:) !! Full finite knots, simple on the base interval.
        integer, intent(in) :: k !! Degree, one through five.
        integer, intent(in) :: m !! Number of sites, including any closing site.
        real(dp), intent(in) :: bounds(2) !! Required base interval endpoints.
        logical, intent(in) :: periodic !! Whether periodic boundary constraints apply.
        integer :: n,maximum
        ok=.false.
        n=size(t)
        maximum=m+k+1
        if (periodic) maximum=m+2*k
        if (n<2*(k+1) .or. n>maximum) return
        if (.not.all(ieee_is_finite(t))) return
        if (any(t(2:)<t(:n-1))) return
        if (t(k+1)/=bounds(1) .or. t(n-k)/=bounds(2)) return
        if (any(t(k+2:n-k)<=t(k+1:n-k-1))) return
        ok=.true.
    end function fixed_knots_valid

    subroutine null_basis(a,z,code)
        implicit none
        real(dp), intent(in) :: a(:,:) !! Finite constraints, at most as many rows as columns.
        real(dp), allocatable, intent(out) :: z(:,:) !! Orthonormal basis of numerical null space.
        integer, intent(out) :: code !! Zero success; three numerical or allocation failure.
        real(dp), allocatable :: left(:,:),right(:,:),sigma(:)
        integer :: m,n,r,info
        code=3
        m=size(a,1)
        n=size(a,2)
        allocate(left(m,m),right(n,n),sigma(m),stat=info)
        if (info/=0) return
        call svd(a,left,sigma,right,info,full_matrices=.true.)
        if (info/=0) return
        r=count(sigma>real(max(m,n),dp)*epsilon(1.0_dp)*maxval(sigma))
        z=transpose(right(r+1:,:))
        code=0
    end subroutine null_basis

    subroutine fixed_smoothing(u,x,w,k,s,t,bounds,periodic,coeff,status)
        implicit none
        real(dp), intent(in) :: u(:) !! Finite strictly increasing physical parameter sites.
        real(dp), intent(in) :: x(:,:) !! Finite coordinates (sites,channels).
        real(dp), intent(in), optional :: w(:) !! Positive residual multipliers, default one.
        integer, intent(in) :: k !! Degree, one through five.
        real(dp), intent(in) :: s !! Positive shared weighted squared-residual target.
        real(dp), intent(in) :: t(:) !! Validated full knots, retained by the caller.
        real(dp), intent(in) :: bounds(2) !! Physical base interval endpoints.
        logical, intent(in) :: periodic !! Enforce matching endpoint derivatives through k-1.
        real(dp), allocatable, intent(out) :: coeff(:,:) !! Unpadded physical coefficients.
        integer, intent(out) :: status !! Zero success; three range, size or numerical failure.
        type(bspline_array_interpolator_type) :: basis,physical
        real(dp), allocatable :: tt(:),uu(:),weights(:),identity(:,:),z(:,:),constraints(:,:)
        real(dp), allocatable :: design(:,:),rhs(:,:),jumps(:,:),a(:,:),d(:,:),v(:,:),c(:,:)
        real(dp), allocatable :: polynomial(:,:),aug(:,:),brhs(:,:),trial(:,:),values(:,:)
        real(dp) :: span,xscale,wscale,logtarget,target,fp,limit,mu,lo,hi,scale
        integer :: m,n,ndim,active,nj,i,j,r,info,iteration,terms
        real(dp) :: elementary(0:5),choose
        logical :: limited
        status=3
        m=size(u)
        n=size(t)-k-1
        ndim=size(x,2)
        ! Dense fallback for moderate problems; automatic banded solvers are unchanged.
        if (n>256 .or. m>10000 .or. ndim>256) return
        if (maxval(abs(t))>huge(1.0_dp)/real(4*(k+1),dp)) return
        span=bounds(2)-bounds(1)
        if (span<tiny(1.0_dp)) return
        allocate(tt(size(t)),uu(m),weights(m),identity(n,n),stat=info)
        if (info/=0) return
        tt=t-bounds(1)
        if (span<1.0_dp) then
            if (maxval(abs(tt))>huge(1.0_dp)*span) return
        end if
        tt=tt/span
        if (maxval(abs(tt))>1e4_dp) return
        uu=(u-bounds(1))/span
        limit=tiny(1.0_dp)**(1.0_dp/real(4*(k+1),dp))
        do i=2,size(tt)
            if (t(i)>t(i-1)) then
                if (t(i)-t(i-1)<tiny(1.0_dp)) return
                if (tt(i)-tt(i-1)<limit) return
            end if
        end do
        if (minval(uu(2:)-uu(:m-1))<limit) return
        active=m
        if (periodic) active=m-1
        weights=1.0_dp
        if (present(w)) weights(:active)=w(:active)
        wscale=maxval(weights(:active))
        weights(:active)=weights(:active)/wscale
        if (minval(weights(:active))<tiny(1.0_dp)**0.125_dp) return
        xscale=maxval(abs(x(:active,:)))
        if (xscale==0.0_dp) xscale=1.0_dp
        if (any(x(:active,:)/=0.0_dp .and. x(:active,:)/xscale==0.0_dp)) return
        logtarget=log(s)-2*log(xscale)-2*log(wscale)
        if (logtarget<log(tiny(1.0_dp))) return
        target=exp(min(logtarget,log(real(active*ndim,dp)+1.0_dp)))
        identity=0.0_dp
        do i=1,n
            identity(i,i)=1.0_dp
        end do
        basis=bspline(tt,identity,k,status=info)
        if (info/=0) return
        if (periodic) then
            allocate(constraints(k,n),stat=info)
            if (info/=0) return
            do i=1,k
                constraints(i,:)=basis%derivative(0.0_dp,i-1)-basis%derivative(1.0_dp,i-1)
                scale=maxval(abs(constraints(i,:)))
                if (scale>0.0_dp) constraints(i,:)=constraints(i,:)/scale
            end do
            call null_basis(constraints,z,info)
            if (info/=0) return
        else
            z=identity
        end if
        design=basis%evaluate(uu(:active))
        do i=1,active
            design(i,:)=design(i,:)*weights(i)
        end do
        rhs=x(:active,:)/xscale*spread(weights(:active),2,ndim)
        nj=size(t)-2*k-2
        allocate(jumps(nj,n),stat=info)
        if (info/=0) return
        do i=1,nj
            j=k+1+i
            jumps(i,:)=basis%derivative((tt(j)+tt(j+1))/2,k)- &
                basis%derivative((tt(j-1)+tt(j))/2,k)
        end do
        a=matmul(design,z)
        d=matmul(jumps,z)
        if (nj>0) then
            scale=maxval(abs(d))
            if (scale>0.0_dp) d=d/scale
        end if
        allocate(v(size(z,2),ndim),stat=info)
        if (info/=0) return
        ! Marsden's identity gives polynomial coefficients without a numerical null space.
        terms=k+1
        if (periodic) terms=1
        allocate(polynomial(n,terms),trial(terms,ndim),stat=info)
        if (info/=0) return
        polynomial(:,1)=1.0_dp
        if (.not.periodic) then
            do i=1,n
                elementary=0.0_dp
                elementary(0)=1.0_dp
                do j=1,k
                    do r=j,1,-1
                        elementary(r)=elementary(r)+tt(i+j)*elementary(r-1)
                    end do
                end do
                choose=1.0_dp
                do r=1,k
                    choose=choose*real(k-r+1,dp)/real(r,dp)
                    polynomial(i,r+1)=elementary(r)/choose
                end do
            end do
        end if
        call lstsq(matmul(design,polynomial),rhs,trial,info)
        if (info/=0) return
        c=matmul(polynomial,trial)
        fp=sum((matmul(design,c)-rhs)**2)
        limited=fp<=target
        deallocate(trial)
        if (limited) goto 800
        call lstsq(a,rhs,v,info)
        if (info/=0) return
        fp=sum((matmul(a,v)-rhs)**2)
        if (fp>1.001_dp*target) return
        ! A feasible positive penalty is required even when the data matrix is singular.
        allocate(aug(active+nj,size(z,2)),brhs(active+nj,ndim),stat=info)
        if (info/=0) return
        brhs=0.0_dp
        brhs(:active,:)=rhs
        lo=0.0_dp
        hi=1.0_dp
        do iteration=1,64
            call penalized(hi)
            if (info/=0) return
            if (fp>=target) exit
            lo=hi
            hi=hi*16.0_dp
        end do
        if (fp<target) return
        do iteration=1,160
            mu=(lo+hi)/2
            call penalized(mu)
            if (info/=0) return
            if (abs(fp-target)<=0.001_dp*target) exit
            if (fp>target) then
                hi=mu
            else
                lo=mu
            end if
        end do
        if (abs(fp-target)>0.001_dp*target) return
800     if (.not.limited) c=matmul(z,v)
        if (.not.all(ieee_is_finite(c))) return
        if (xscale>1.0_dp) then
            if (maxval(abs(c))>huge(1.0_dp)/xscale) return
        end if
        physical=bspline(t,c*xscale,k,status=info)
        if (info/=0) return
        values=physical%evaluate(u(:active))/xscale
        if (.not.all(ieee_is_finite(values))) return
        fp=sum(((values-x(:active,:)/xscale)*spread(weights(:active),2,ndim))**2)
        if (fp>1.0011_dp*target) return
        coeff=c*xscale
        status=0
    contains
        subroutine penalized(multiplier)
            implicit none
            real(dp), intent(in) :: multiplier !! Square root of the shared penalty multiplier.
            aug(:active,:)=a
            aug(active+1:,:)=multiplier*d
            call lstsq(aug,brhs,v,info)
            if (info/=0) return
            fp=sum((matmul(a,v)-rhs)**2)
        end subroutine penalized
    end subroutine fixed_smoothing
end module scifort_interpolate_fixed_smoothing
