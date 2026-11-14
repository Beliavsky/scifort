! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
! Independent cyclic B-spline least squares; see CODE_PROVENANCE.md.
module scifort_interpolate_periodic_lsq
    use scifort_kinds, only: dp
    use scifort_interpolate_bspline, only: bspline,bspline_array_interpolator_type
    use scifort_linalg_qr, only: lstsq
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: periodic_lsq
contains
    subroutine periodic_lsq(u,x,w,k,interior,knots,coeff,status)
        implicit none
        real(dp), intent(in) :: u(:) !! Finite increasing sites, including the period endpoint.
        real(dp), intent(in) :: x(:,:) !! Sample coordinates; final row is ignored.
        real(dp), intent(in), optional :: w(:) !! Positive active residual weights; last ignored.
        integer, intent(in) :: k !! Spline degree, one through five.
        real(dp), intent(in) :: interior(:) !! Nonempty increasing knots strictly within the period.
        real(dp), allocatable, intent(out) :: knots(:) !! Canonical full periodic knots on success.
        real(dp), allocatable, intent(out) :: coeff(:,:) !! Unpadded cyclic coefficients by channel.
        integer, intent(out) :: status !! Zero success; one invalid; three numerical or size failure.
        type(bspline_array_interpolator_type) :: basis
        real(dp), allocatable :: tt(:),t(:),uu(:),weights(:),wrap(:,:),a(:,:),rhs(:,:),c(:,:)
        real(dp), allocatable :: solution(:,:)
        real(dp) :: span,xscale,wscale,gap
        integer :: m,ndim,ni,n,nc,i,info
        status=1
        m=size(u)
        ndim=size(x,2)
        if (k<1 .or. k>5 .or. m<=k .or. ndim<1) return
        if (size(x,1)/=m) return
        if (size(interior)<1 .or. size(interior)>m-2) return
        if (.not.all(ieee_is_finite(u)) .or. .not.all(ieee_is_finite(interior))) return
        if (.not.all(ieee_is_finite(x(:m-1,:)))) return
        if (any(u(2:)<=u(:m-1))) return
        if (any(interior<=u(1)) .or. any(interior>=u(m))) return
        if (any(interior(2:)<=interior(:size(interior)-1))) return
        if (present(w)) then
            if (size(w)/=m) return
            if (.not.all(ieee_is_finite(w(:m-1)))) return
            if (any(w(:m-1)<=0.0_dp)) return
        end if
        status=3
        ni=size(interior)+1
        if (ni>256 .or. m>10000 .or. ndim>256) return
        if (maxval(abs(u))>huge(1.0_dp)/real(4*(k+1),dp)) return
        span=u(m)-u(1)
        if (span<tiny(1.0_dp)) return
        n=ni+2*k+1
        nc=ni+k
        allocate(t(n),tt(n),uu(m),weights(m-1),wrap(nc,ni),solution(ni,ndim), &
            c(nc,ndim),stat=info)
        if (info/=0) return
        t(k+1:n-k)=[u(1),interior,u(m)]
        ! Sequential extension also handles cycles shorter than the degree.
        do i=1,k
            t(k+1-i)=t(k+1-i+ni)-span
            t(n-k+i)=t(k+1+i)+span
        end do
        if (.not.all(ieee_is_finite(t))) return
        if (minval(t(2:)-t(:n-1))<tiny(1.0_dp)) return
        tt=(t-u(1))/span
        uu=(u-u(1))/span
        gap=tiny(1.0_dp)**(1.0_dp/real(4*(k+1),dp))
        if (minval(tt(2:)-tt(:n-1))<gap) return
        if (minval(uu(2:)-uu(:m-1))<gap) return
        weights=1.0_dp
        if (present(w)) weights=w(:m-1)
        wscale=maxval(weights)
        weights=weights/wscale
        if (minval(weights)<tiny(1.0_dp)**0.125_dp) return
        xscale=maxval(abs(x(:m-1,:)))
        if (xscale==0.0_dp) xscale=1.0_dp
        if (any(x(:m-1,:)/=0.0_dp .and. x(:m-1,:)/xscale==0.0_dp)) return
        wrap=0.0_dp
        do i=1,nc
            wrap(i,modulo(i-1,ni)+1)=1.0_dp
        end do
        basis=bspline(tt,wrap,k,status=info)
        if (info/=0) return
        a=basis%evaluate(uu(:m-1))
        do i=1,m-1
            a(i,:)=a(i,:)*weights(i)
        end do
        rhs=x(:m-1,:)/xscale*spread(weights,2,ndim)
        call lstsq(a,rhs,solution,info)
        if (info/=0 .or. .not.all(ieee_is_finite(solution))) return
        if (xscale>1.0_dp) then
            if (maxval(abs(solution))>huge(1.0_dp)/xscale) return
        end if
        do i=1,nc
            c(i,:)=solution(modulo(i-1,ni)+1,:)*xscale
            if (any(solution(modulo(i-1,ni)+1,:)/=0.0_dp .and. c(i,:)==0.0_dp)) return
        end do
        call move_alloc(t,knots)
        call move_alloc(c,coeff)
        status=0
    end subroutine periodic_lsq
end module scifort_interpolate_periodic_lsq
