! SPDX-License-Identifier: MIT AND BSD-3-Clause
! Copyright (c) 2001-2002 Enthought, Inc. 2003, SciPy Developers.
! Copyright (c) 2026 SciFort contributors.
! Adapted selection logic from SciPy v1.18.0 _fitpack_repro.py,
! _bsplines.py and src/__fitpack.cc; see CODE_PROVENANCE.md and
! LICENSES/scipy-generate-knots-BSD-3-Clause.txt.
module scifort_interpolate_generate_knots
    use scifort_kinds, only: dp
    use scifort_interpolate_bspline, only: bspline_design_matrix
    use scifort_linalg_qr, only: lstsq
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: generate_knots, knot_generator_type

    type :: knot_generator_type
        private
        real(dp), allocatable :: x(:),u(:),y(:,:),w(:),t(:)
        real(dp) :: origin=0.0_dp,span=1.0_dp,target=0.0_dp,fp=0.0_dp
        integer :: k=3,nest=0,nmax=0,nplus=1
        logical :: periodic=.false.,ready=.false.,first=.true.,terminal=.false.
    contains
        procedure :: next => next_knots
    end type knot_generator_type

    interface generate_knots
        module procedure generate_scalar,generate_channels
    end interface generate_knots
contains
    function generate_scalar(x,y,w,xb,xe,k,s,nest,bc_type,status) result(gen)
        implicit none
        real(dp), intent(in) :: x(:) !! Finite strictly increasing sample sites, at most 256.
        real(dp), intent(in) :: y(:) !! Finite observations, one per sample site.
        real(dp), intent(in), optional :: w(:) !! Positive residual weights, default one.
        real(dp), intent(in), optional :: xb,xe !! Enclosing bounds, default first/last site.
        integer, intent(in), optional :: k !! Degree one through five, default three.
        real(dp), intent(in), optional :: s !! Nonnegative shared residual target, default zero.
        integer, intent(in), optional :: nest !! Soft stopping count, at least 2*(k+1).
        character(*), intent(in), optional :: bc_type !! 'not-a-knot' (default) or 'periodic'.
        integer, intent(out), optional :: status !! 0 success, 1 invalid, 3 range/size/numerical.
        type(knot_generator_type) :: gen
        gen=generate_channels(x,reshape(y,[size(y),1]),w,xb,xe,k,s,nest,bc_type,status)
    end function generate_scalar

    function generate_channels(x,y,w,xb,xe,k,s,nest,bc_type,status) result(gen)
        implicit none
        real(dp), intent(in) :: x(:) !! Finite strictly increasing sites, at most 256.
        real(dp), intent(in) :: y(:,:) !! Finite observations (sites, channels), at most 256 channels.
        real(dp), intent(in), optional :: w(:) !! Positive residual weights, default one.
        real(dp), intent(in), optional :: xb,xe !! Enclosing bounds; periodic bounds equal endpoints.
        integer, intent(in), optional :: k !! Degree one through five, default three.
        real(dp), intent(in), optional :: s !! Nonnegative sum of weighted squared errors, default zero.
        integer, intent(in), optional :: nest !! Soft knot count; the initial refinement may exceed it.
        character(*), intent(in), optional :: bc_type !! 'not-a-knot' (default) or 'periodic'.
        integer, intent(out), optional :: status !! 0 success, 1 invalid, 3 range/size/numerical.
        type(knot_generator_type) :: gen
        real(dp), allocatable :: residuals(:)
        real(dp) :: target,left,right,yscale,wscale,logtarget,gap,scale
        integer :: m,ndim,code,j,active
        code=1
        m=size(x)
        ndim=size(y,2)
        if (present(k)) gen%k=k
        if (gen%k<1 .or. gen%k>5 .or. m<=gen%k .or. ndim<1) goto 900
        if (size(y,1)/=m) goto 900
        if (.not.all(ieee_is_finite(x)) .or. .not.all(ieee_is_finite(y))) goto 900
        if (any(x(2:)<=x(:m-1))) goto 900
        if (present(bc_type)) then
            select case(trim(bc_type))
            case('periodic')
                gen%periodic=.true.
            case('not-a-knot')
            case default
                goto 900
            end select
        end if
        target=0.0_dp
        if (present(s)) target=s
        if (.not.ieee_is_finite(target)) goto 900
        if (target<0.0_dp) goto 900
        if (target==0.0_dp .and. (present(w) .or. present(nest))) goto 900
        left=x(1)
        right=x(m)
        if (present(xb)) left=xb
        if (present(xe)) right=xe
        if (.not.ieee_is_finite(left) .or. .not.ieee_is_finite(right)) goto 900
        if (left>x(1) .or. right<x(m)) goto 900
        if (gen%periodic) then
            if (left/=x(1) .or. right/=x(m)) goto 900
            do j=1,ndim
                scale=max(1.0_dp,abs(y(1,j)),abs(y(m,j)))
                if (abs(y(1,j)/scale-y(m,j)/scale)>1.0e-12_dp) goto 900
            end do
        end if
        if (present(w)) then
            if (size(w)/=m) goto 900
            if (.not.all(ieee_is_finite(w))) goto 900
            if (any(w<=0.0_dp)) goto 900
        end if
        gen%nmax=m+gen%k+1
        if (gen%periodic) gen%nmax=m+2*gen%k
        gen%nest=max(gen%nmax,2*gen%k+3)
        if (present(nest)) then
            if (nest<2*(gen%k+1)) goto 900
            gen%nest=nest
        end if
        code=3
        if (m>256 .or. ndim>256) goto 900
        if (max(abs(left),abs(right))>huge(1.0_dp)/real(4*(gen%k+1),dp)) goto 900
        gen%origin=left
        gen%span=right-left
        if (gen%span<tiny(1.0_dp)) goto 900
        gen%x=x
        gen%u=(x-left)/gen%span
        gap=tiny(1.0_dp)**(1.0_dp/real(4*(gen%k+1),dp))
        if (minval(gen%u(2:)-gen%u(:m-1))<gap) goto 900
        yscale=maxval(abs(y))
        if (yscale==0.0_dp) yscale=1.0_dp
        gen%y=y/yscale
        if (any(y/=0.0_dp .and. gen%y==0.0_dp)) goto 900
        allocate(gen%w(m))
        gen%w=1.0_dp
        if (present(w)) gen%w=w
        active=m
        if (gen%periodic) then
            active=m-1
            gen%w(m)=0.0_dp
        end if
        wscale=maxval(gen%w)
        gen%w=gen%w/wscale
        if (minval(gen%w(:active))<tiny(1.0_dp)**0.125_dp) goto 900
        if (target>0.0_dp) then
            logtarget=log(target)-2.0_dp*(log(yscale)+log(wscale))
            if (logtarget<log(tiny(1.0_dp))) goto 900
            gen%target=exp(min(logtarget,log(huge(1.0_dp)/16.0_dp)))
        end if
        if (target==0.0_dp) then
            call interpolation_knots(gen)
            gen%terminal=.true.
        else
            allocate(gen%t(2*(gen%k+1)))
            gen%t(:gen%k+1)=left
            gen%t(gen%k+2:)=right
            if (.not.gen%periodic .and. gen%nmax==size(gen%t)) gen%terminal=.true.
            if (gen%periodic) then
                call wrap_knots(gen)
                call fit_residuals(gen,residuals,gen%fp,code)
                if (code/=0) goto 900
                if (gen%fp-gen%target<0.001_dp*gen%target .or. &
                    gen%nmax==size(gen%t)) then
                    gen%terminal=.true.
                else
                    gen%t=[gen%t(:gen%k+1),x((m+1)/2),gen%t(gen%k+2:)]
                    call wrap_knots(gen)
                end if
            end if
        end if
        gen%ready=.true.
        code=0
900     if (present(status)) status=code
    end function generate_channels

    subroutine next_knots(self,t,done,status)
        implicit none
        class(knot_generator_type), intent(inout) :: self !! Initialized generator; advances one step.
        real(dp), allocatable, intent(out) :: t(:) !! Full knots; empty at completion or failure.
        logical, intent(out) :: done !! True when no vector is returned, including failure.
        integer, intent(out) :: status !! 0 vector/completion, 1 uninitialized, 3 numerical failure.
        real(dp), allocatable :: residuals(:)
        real(dp) :: oldfp,delta,excess,ratio,fp
        integer :: code,j,n,proposed
        allocate(t(0))
        done=.true.
        status=1
        if (.not.self%ready) return
        status=0
        if (self%first) then
            self%first=.false.
        else
            if (self%terminal) return
            oldfp=self%fp
            call fit_residuals(self,residuals,fp,code)
            if (code/=0) goto 800
            self%fp=fp
            excess=fp-self%target
            if (excess<0.001_dp*self%target) then
                self%terminal=.true.
                return
            end if
            n=size(self%t)
            if (n==2*(self%k+1)) then
                self%nplus=1
            else
                delta=oldfp-fp
                proposed=2*self%nplus
                if (delta>0.001_dp*self%target) then
                    ! Saturate before integer conversion; only this bounded range matters.
                    if (excess<2.0_dp*delta) then
                        ratio=excess/delta
                        proposed=int(real(self%nplus,dp)*ratio)
                    end if
                end if
                self%nplus=min(2*self%nplus,max(proposed,self%nplus/2,1))
            end if
            do j=1,self%nplus
                call insert_residual_knot(self,residuals,code)
                if (code/=0) goto 800
                if (size(self%t)>=self%nmax) then
                    call interpolation_knots(self)
                    self%terminal=.true.
                    exit
                end if
                if (size(self%t)>=self%nest) then
                    self%terminal=.true.
                    exit
                end if
                if (j<self%nplus) then
                    call fit_residuals(self,residuals,fp,code)
                    if (code/=0) goto 800
                end if
            end do
        end if
        n=size(self%t)
        if (.not.all(ieee_is_finite(self%t))) goto 800
        if (any(self%t(self%k+2:n-self%k)<=self%t(self%k+1:n-self%k-1))) goto 800
        if (self%periodic .and. self%target>0.0_dp) then
            if (any(self%t(2:)<=self%t(:n-1))) goto 800
        end if
        t=self%t
        done=.false.
        return
800     status=3
        self%ready=.false.
    end subroutine next_knots

    subroutine fit_residuals(gen,residuals,fp,status)
        implicit none
        class(knot_generator_type), intent(in) :: gen !! Normalized samples and physical full knots.
        real(dp), allocatable, intent(out) :: residuals(:) !! Weighted squared errors by sample.
        real(dp), intent(out) :: fp !! Sum over samples and channels; periodic last row ignored.
        integer, intent(out) :: status !! 0 success; 3 numerical failure.
        real(dp), allocatable :: a(:,:),b(:,:),c(:,:),basis(:,:),rhs(:,:),error(:,:)
        integer :: nc,ni,m,j,info
        status=3
        fp=0.0_dp
        m=size(gen%x)
        if (gen%periodic) m=m-1
        nc=size(gen%t)-gen%k-1
        ni=nc
        if (gen%periodic) ni=nc-gen%k
        basis=bspline_design_matrix(gen%u(:m),(gen%t-gen%origin)/gen%span,gen%k,info)
        if (info/=0) return
        allocate(a(m,ni),c(ni,size(gen%y,2)))
        a=0.0_dp
        do j=1,nc
            a(:,modulo(j-1,ni)+1)=a(:,modulo(j-1,ni)+1)+basis(:,j)
        end do
        b=a*spread(gen%w(:m),2,ni)
        rhs=gen%y(:m,:)*spread(gen%w(:m),2,size(gen%y,2))
        call lstsq(b,rhs,c,info)
        if (info/=0 .or. .not.all(ieee_is_finite(c))) return
        if (maxval(abs(c))>1.0e100_dp) return
        error=matmul(b,c)-rhs
        allocate(residuals(size(gen%x)))
        residuals=0.0_dp
        residuals(:m)=sum(error**2,dim=2)
        fp=sum(residuals)
        if (.not.ieee_is_finite(fp)) return
        status=0
    end subroutine fit_residuals

    subroutine insert_residual_knot(gen,residuals,status)
        implicit none
        class(knot_generator_type), intent(inout) :: gen !! Current knots; inserts one data-site knot.
        real(dp), intent(in) :: residuals(:) !! Weighted squared sample residuals.
        integer, intent(out) :: status !! 0 success; 3 if no splittable interval remains.
        real(dp) :: parts(size(gen%x)),total,best,new
        integer :: ends(size(gen%x)+1),interval,nc,i,nruns,chosen,index,pos
        interval=gen%k+2
        nc=size(gen%t)-gen%k-1
        nruns=0
        ends(1)=1
        total=0.0_dp
        do i=1,size(gen%x)
            total=total+residuals(i)
            if (gen%x(i)>=gen%t(interval) .and. interval<=nc) then
                nruns=nruns+1
                parts(nruns)=total-0.5_dp*residuals(i)
                total=0.5_dp*residuals(i)
                interval=interval+1
                ends(nruns+1)=i
            end if
        end do
        nruns=nruns+1
        ends(nruns+1)=size(gen%x)
        parts(nruns)=total
        chosen=0
        best=-1.0_dp
        do i=1,nruns
            if (ends(i+1)-ends(i)<=1) cycle
            if (parts(i)<=best) cycle
            best=parts(i)
            chosen=i
        end do
        status=3
        if (chosen==0) return
        index=(ends(chosen)+ends(chosen+1)+1)/2
        new=gen%x(index)
        if (any(gen%t==new)) return
        pos=count(gen%t<new)
        gen%t=[gen%t(:pos),new,gen%t(pos+1:)]
        if (gen%periodic) call wrap_knots(gen)
        status=0
    end subroutine insert_residual_knot

    subroutine wrap_knots(gen)
        implicit none
        class(knot_generator_type), intent(inout) :: gen !! Knots with base endpoints and interiors set.
        integer :: i,n,ni,k
        k=gen%k
        n=size(gen%t)
        ni=n-2*k-1
        do i=1,k
            gen%t(k+1-i)=gen%t(k+1-i+ni)-gen%span
            gen%t(n-k+i)=gen%t(k+1+i)+gen%span
        end do
    end subroutine wrap_knots

    subroutine interpolation_knots(gen)
        implicit none
        class(knot_generator_type), intent(inout) :: gen !! Replaces knots with interpolation convention.
        real(dp), allocatable :: core(:)
        integer :: m,k,i,j
        m=size(gen%x)
        k=gen%k
        if (allocated(gen%t)) deallocate(gen%t)
        if (gen%periodic .and. (k>1 .or. gen%target>0.0_dp)) then
            core=gen%x
            if (modulo(k,2)==0) core(2:m-1)= &
                gen%x(2:m-1)-0.5_dp*(gen%x(2:m-1)-gen%x(:m-2))
            allocate(gen%t(m+2*k))
            gen%t(k+1:m+k)=core
            call wrap_knots(gen)
        else
            allocate(gen%t(m+k+1))
            gen%t(:k+1)=gen%x(1)
            gen%t(m+1:)=gen%x(m)
            do i=1,m-k-1
                if (modulo(k,2)==1) then
                    gen%t(k+1+i)=gen%x(i+(k+1)/2)
                else
                    j=i+k/2
                    gen%t(k+1+i)=gen%x(j)+0.5_dp*(gen%x(j+1)-gen%x(j))
                end if
            end do
        end if
    end subroutine interpolation_knots
end module scifort_interpolate_generate_knots
