! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
! Reference driver for the separately downloaded pinned perazz/FITPACK core.
program upstream_grid
 use fitpack_core, only: spgrid
 use scifort_kinds, only: dp
 use scifort_interpolate_fitpack_bivariate, only: fitpack_bivariate_tck_type,bisplev
 use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
 implicit none
 type(fitpack_bivariate_tck_type) :: tck
 real(dp),allocatable :: u(:),v(:),r(:,:),qt(:),qp(:),tu(:),tv(:),c(:),wrk(:),data(:)
 integer,allocatable :: iwrk(:)
 real(dp) :: poles(2),s,fp
 integer :: mu,mv,nq,i,flags(6),iopt(3),ider(4),nu,nv,nuest,nvest,lwrk,kwrk,ier
 read(*,*) mu,mv,nq
 allocate(u(mu),v(mv),r(mu,mv),qt(nq),qp(nq))
 read(*,*) u
 read(*,*) v
 do i=1,mu
  read(*,*) r(i,:)
 end do
 read(*,*) flags
 read(*,*) poles
 do i=1,nq
  read(*,*) qt(i),qp(i)
 end do
 read(*,*) s
 iopt=[0,flags(1:2)]
 ider=[flags(3),flags(5),flags(4),flags(6)]
 do i=1,2
  if(ieee_is_nan(poles(i))) then
   ider(2*i-1)=-1
   poles(i)=0
  end if
 end do
 nuest=mu+8
 nvest=mv+7
 lwrk=12+nuest*(mv+nvest+3)+24*nvest+4*mu+8*mv+max(nuest,mv+nvest)
 kwrk=5+mu+mv+nuest+nvest
 allocate(tu(nuest),tv(nvest),c((nuest-4)*(nvest-4)),wrk(lwrk),iwrk(kwrk))
 tu=0
 tv=0
 c=0
 wrk=0
 iwrk=0
 nu=0
 nv=0
 fp=0
 data=reshape(transpose(r),[mu*mv])
 call spgrid(iopt,ider,mu,u,mv,v,data,poles(1),poles(2),s,nuest,nvest, &
  nu,tu,nv,tv,c,fp,wrk,lwrk,iwrk,kwrk,ier)
 print *, max(0,ier),-1,0,fp
 if(ier>0) stop
 tck%tx=tu(:nu)
 tck%ty=tv(:nv)
 tck%c=c(:(nu-4)*(nv-4))
 tck%kx=3
 tck%ky=3
 print '(*(es26.17e3,1x))',tck%tx
 print '(*(es26.17e3,1x))',tck%ty
 print '(*(es26.17e3,1x))',tck%c
 do i=1,nq
  print '(6(es26.17e3,1x))',bisplev(qt(i),qp(i),tck), &
   bisplev(qt(i),qp(i),tck,dx=1),bisplev(qt(i),qp(i),tck,dy=1), &
   bisplev(qt(i),qp(i),tck,dx=2),bisplev(qt(i),qp(i),tck,dy=2), &
   bisplev(qt(i),qp(i),tck,dx=1,dy=1)
 end do
end program
