! SPDX-License-Identifier: BSD-3-Clause
! Copyright (c) 2022, Federico Perini. All rights reserved.
! Copyright (c) 2026 SciFort contributors
!
! Adapted from perazz/fitpack src/fitpack_core.F90, commit
! 1907c7bad9f5c4557970e906bdd4909eb52b2d4c, by Federico Perini,
! based on Paul Dierckx's FITPACK. Retained license:
! LICENSES/perazz-fitpack-BSD-3-Clause.txt.
! Only the s=0 grid path and its numerical helpers are included.
module scifort_interpolate_sphere_grid_core
    use scifort_kinds, only : dp
    use scifort_constants, only : pi => scifort_pi
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    implicit none
    private
    public :: fitpack_sphere_grid
    integer, parameter :: max_order = 19, degree_3 = 3, degree_4 = 4, degree_5 = 5
    real(dp), parameter :: zero = 0.0_dp, one = 1.0_dp, two = 2.0_dp
    real(dp), parameter :: three = 3.0_dp, half = 0.5_dp, pi2 = 2.0_dp*pi
    real(dp), parameter :: third = 1.0_dp/3.0_dp
    logical, parameter :: fp_true = .true., fp_false = .false.
contains
    subroutine fitpack_sphere_grid(u,v,r,continuity,poles,given,exact,flat,tu,tv,c,status)
        implicit none
        real(dp), intent(in) :: u(:) !! Increasing interior colatitude grid, at least four sites.
        real(dp), intent(in) :: v(:) !! Increasing periodic longitude grid, at least four sites.
        real(dp), intent(in) :: r(:,:) !! Finite scaled samples, shape (size(u),size(v)).
        logical, intent(in) :: continuity(2) !! First derivative continuity at north/south poles.
        real(dp), intent(in) :: poles(2) !! Scaled supplied pole values, ignored unless given.
        logical, intent(in) :: given(2) !! Whether each pole value was supplied.
        logical, intent(in) :: exact(2) !! Exactness of supplied pole values.
        logical, intent(in) :: flat(2) !! Zero first derivatives at north/south poles.
        real(dp), allocatable, intent(out) :: tu(:) !! Full active colatitude knots.
        real(dp), allocatable, intent(out) :: tv(:) !! Full active longitude knots.
        real(dp), allocatable, intent(out) :: c(:) !! Longitude-fastest tensor coefficients.
        integer, intent(out) :: status !! 0 success, 3 arithmetic failure; caller validates inputs.
        real(dp), allocatable :: xknots(:),yknots(:),coef(:),data(:),wrk(:),fpu(:),fpv(:)
        integer, allocatable :: nru(:),nrv(:)
        integer :: mu,mv,mr,nuest,nvest,nc,nu,nv,numax,lwrk,ifsu,ifsv,ifbu,ifbv
        integer :: iopt(3),ider(4),idd(4),i,ju,nn
        real(dp) :: dr(6),step(2),r0,r1,fp,p,range
        allocate(tu(0),tv(0),c(0))
        status = 3
        mu = size(u)
        mv = size(v)
        if (mu < 4 .or. mv < 4) return
        if (size(r,1) /= mu .or. size(r,2) /= mv) return
        iopt = [0,merge(1,0,continuity(1)),merge(1,0,continuity(2))]
        ider = [-1,merge(1,0,flat(1)),-1,merge(1,0,flat(2))]
        if (given(1)) ider(1) = merge(1,0,exact(1))
        if (given(2)) ider(3) = merge(1,0,exact(2))
        r0 = sum(r(1,:))/real(mv,dp)
        r1 = sum(r(mu,:))/real(mv,dp)
        if (given(1)) r0 = poles(1)
        if (given(2)) r1 = poles(2)
        dr = [r0,zero,zero,r1,zero,zero]
        range = max(r0,r1,maxval(r))-min(r0,r1,minval(r))
        step = range
        mr = mu*mv
        numax = mu+6+iopt(2)+iopt(3)
        nuest = max(8,numax)
        nvest = mv+7
        nc = (nuest-4)*(nvest-4)
        lwrk = 12+nuest*(mv+nvest+3)+24*nvest+4*mu+8*mv+max(nuest,mv+nvest)
        allocate(xknots(nuest),yknots(nvest),coef(nc),data(mr), &
            wrk(lwrk),fpu(nuest),fpv(nvest),nru(mu),nrv(mv))
        xknots = zero
        yknots = zero
        coef = zero
        wrk = zero
        fpu = zero
        fpv = zero
        nru = 0
        nrv = 0
        data = reshape(transpose(r),[mr])
        nv = nvest
        yknots = [v(mv-2:mv)-pi2,v,v(1)+pi2,v(2:4)+pi2]
        ifsu = 0
        ifsv = 0
        ifbu = 0
        ifbv = 0
        p = -one
        fp = zero
        idd = ider
        if (idd(1) == 0) idd(1) = 1
        if (idd(3) == 0) idd(3) = 1
        if (ider(1) < 0 .or. (iopt(2) /= 0 .and. ider(2) == 0) .or. &
            ider(3) < 0 .or. (iopt(3) /= 0 .and. ider(4) == 0)) then
            nu = 4
            do i = 3, mu-2, 2
                nu = nu+1
                xknots(nu) = u(i)
            end do
            xknots(nu+1:nu+4) = pi
            nu = nu+4
            call fpopsp(ifsu,ifsv,ifbu,ifbv,u,mu,v,mv,data,mr,r0,r1,dr,iopt,idd, &
                xknots,nu,yknots,nv,nuest,nvest,p,step,coef,nc,fp,fpu,fpv,nru,nrv,wrk,lwrk)
            ifsu = 0
        end if
        nu = numax
        nn = nu-8
        xknots = zero
        ju = 2-iopt(2)
        xknots(5:nn+4) = u(ju:ju+nn-1)
        xknots(nu-3:nu) = pi
        call fpopsp(ifsu,ifsv,ifbu,ifbv,u,mu,v,mv,data,mr,r0,r1,dr,iopt,idd, &
            xknots,nu,yknots,nv,nuest,nvest,p,step,coef,nc,fp,fpu,fpv,nru,nrv,wrk,lwrk)
        nc = (nu-4)*(nv-4)
        if (.not. all(ieee_is_finite(coef(:nc))) .or. .not. ieee_is_finite(fp)) return
        tu = xknots(:nu)
        tv = yknots(:nv)
        c = coef(:nc)
        status = 0
    end subroutine fitpack_sphere_grid

    elemental logical function equal(a,b) result(ok)
        implicit none
        real(dp), intent(in) :: a !! First value in the spacing-based comparison.
        real(dp), intent(in) :: b !! Second value in the spacing-based comparison.
       ok = abs(a-b)<spacing(merge(a,b,abs(a)<abs(b)))
    end function equal

      pure function fp_knot_interval(t, x, l_start, l_max) result(l)
        implicit none
        real(dp),    intent(in) :: t(:) !! Nondecreasing B-spline knot vector.
        real(dp),    intent(in) :: x !! Query coordinate in the active knot span.
        integer, intent(in) :: l_start !! Initial interval index for the monotone search.
        integer, intent(in) :: l_max !! Last permitted knot interval index.
          integer :: l

          integer, parameter :: linear_threshold = 8
          integer :: lo, hi, mid

          ! for small ranges, use linear search (also optimal for warm starts)
          if (l_max - l_start <= linear_threshold) then
              l = l_start
              do while (x >= t(l+1) .and. l < l_max)
                  l = l + 1
              end do
              return
          end if

          ! handle boundary cases
          if (x < t(l_start+1)) then
              l = l_start
              return
          end if
          if (x >= t(l_max)) then
              l = l_max
              return
          end if

          ! binary search: find largest l such that t(l+1) <= x, then l is the interval
          lo = l_start
          hi = l_max
          do while (hi > lo)
              mid = (lo + hi) / 2
              if (t(mid+1) <= x) then
                  lo = mid + 1
              else
                  hi = mid
              end if
          end do
          l = lo
      end function fp_knot_interval

      pure subroutine fp_rotate_2mat_stride(h1, band1, h2, band2, a1, a2, na, &
                                             right, c, ldc, nrhs, j1_start, j1_end)
        implicit none
        integer, intent(in) :: band1 !! Band width of the first periodic factor block.
        integer, intent(in) :: band2 !! Band width of the wraparound factor block.
        integer, intent(in) :: na !! Leading dimension of the band-factor arrays.
        integer, intent(in) :: ldc !! Leading dimension of the right-hand-side storage.
        integer, intent(in) :: nrhs !! Number of right-hand-side channels.
        integer, intent(in) :: j1_start !! First row of the main periodic factor block.
        integer, intent(in) :: j1_end !! Last row of the main periodic factor block.
        real(dp),    intent(inout) :: h1(band1) !! Active row for the first periodic factor block.
        real(dp),    intent(inout) :: h2(band2) !! Active row for the wraparound factor block.
        real(dp),    intent(inout) :: a1(na,*) !! Main band block updated by plane rotations.
        real(dp),    intent(inout) :: a2(na,*) !! Wraparound band block updated by plane rotations.
        real(dp),    intent(inout) :: &
            right(nrhs) !! Right-hand-side row transformed during rotations.
        real(dp),    intent(inout) :: &
            c(ldc,*) !! Strided right-hand-side storage updated by rotations.

          real(dp) :: cos, sin, piv
          integer :: i2, ij, j

          ! phase 1: rotate h1 through a1 with shifting, cross-rotating h2/a2
          do j = j1_start, j1_end
              piv = h1(1)
              i2 = min(j1_end - j, band1 - 1)
              if (not_equal(piv, zero)) then
                  call fpgivs(piv, a1(j,1), cos, sin)
                  call fprota(cos, sin, right(1:nrhs), c(j, 1:nrhs))
                  call fprota(cos, sin, h2(1:band2), a2(j, 1:band2))
                  if (i2 > 0) call fprota(cos, sin, h1(2:i2+1), a1(j, 2:i2+1))
              endif
              h1(1:i2+1) = [h1(2:i2+1), zero]
          end do

          ! phase 2: rotate h2 through a2 (diagonal at a2(ij,j))
          do j = 1, band2
              ij = j1_end + j
              if (ij <= 0) cycle
              piv = h2(j)
              if (equal(piv, zero)) cycle
              call fpgivs(piv, a2(ij,j), cos, sin)
              call fprota(cos, sin, right(1:nrhs), c(ij, 1:nrhs))
              if (j < band2) call fprota(cos, sin, h2(j+1:band2), a2(ij, j+1:band2))
          end do
      end subroutine fp_rotate_2mat_stride

      pure subroutine fp_rotate_row_block(h, band, a, na, right, q, nrhs, irot_start)
        implicit none
        integer, intent(in) :: band !! Number of active diagonals in the band factor.
        integer, intent(in) :: na !! Leading dimension of the band-factor arrays.
        integer, intent(in) :: nrhs !! Number of right-hand-side channels.
        integer, intent(in) :: irot_start !! Number of already eliminated rows.
        real(dp),    intent(inout) :: &
            h(band) !! Active row rotated into the banded factor.
        real(dp),    intent(inout) :: a(na,*) !! Banded upper-triangular factor updated in place.
        real(dp),    intent(inout) :: &
            right(nrhs) !! Right-hand-side row transformed during rotations.
        real(dp),    intent(inout) :: q(nrhs,*) !! Transformed observation block workspace.

          real(dp) :: cos, sin, piv
          integer :: i, irot

          irot = irot_start
          do i = 1, band
              irot = irot + 1
              piv = h(i)
              if (equal(piv, zero)) cycle
              call fpgivs(piv, a(irot,1), cos, sin)
              call fprota(cos, sin, right(1:nrhs), q(1:nrhs, irot))
              if (i < band) call fprota(cos, sin, h(i+1:band), a(irot, 2:band-i+1))
          end do
      end subroutine fp_rotate_row_block

      pure subroutine fp_rotate_row_stride(h, band, a, na, right, c, ldc, nrhs, irot_start)
        implicit none
        integer, intent(in) :: band !! Number of active diagonals in the band factor.
        integer, intent(in) :: na !! Leading dimension of the band-factor arrays.
        integer, intent(in) :: ldc !! Leading dimension of the right-hand-side storage.
        integer, intent(in) :: nrhs !! Number of right-hand-side channels.
        integer, intent(in) :: irot_start !! Number of already eliminated rows.
        real(dp),    intent(inout) :: &
            h(band) !! Active row rotated into the banded factor.
        real(dp),    intent(inout) :: a(na,*) !! Banded upper-triangular factor updated in place.
        real(dp),    intent(inout) :: &
            right(nrhs) !! Right-hand-side row transformed during rotations.
        real(dp),    intent(inout) :: &
            c(ldc,*) !! Strided right-hand-side storage updated by rotations.

          real(dp) :: cos, sin, piv
          integer :: i, irot

          irot = irot_start
          do i = 1, band
              irot = irot + 1
              piv = h(i)
              if (equal(piv, zero)) cycle
              call fpgivs(piv, a(irot,1), cos, sin)
              call fprota(cos, sin, right(1:nrhs), c(irot, 1:nrhs))
              if (i < band) call fprota(cos, sin, h(i+1:band), a(irot, 2:band-i+1))
          end do
      end subroutine fp_rotate_row_stride

      pure function fpback(a,z,n,k,nest) result(c)
        implicit none

      !  ..scalar arguments..
        integer, intent(in) :: n !! Order of the upper-triangular system.
        integer, intent(in) :: k !! Number of diagonals stored in the upper band.
        integer, intent(in) :: nest !! Leading storage dimension of the triangular factor.
      !  ..array arguments..
        real(dp), intent(in) :: a(nest,k) !! Upper-triangular band factor, diagonal in column one.
        real(dp), intent(in) :: z(n) !! Right-hand side of the triangular solve.
      real(dp)              :: c(n)
      !  ..local scalars..
      real(dp) :: store
      integer :: i,i1,j,k1,l,m,jm1
      !  ..
      k1 = k-1

      if (equal(a(n,1),zero)) then
         c(n) = zero
      else
         c(n) = z(n)/a(n,1)
      end if

      i = n-1
      if (i>0) then

         jm1 = 1
         rows: do j=2,n
           store = z(i)
           i1 = merge(jm1,k1,j<=k1)
           m = i
           do l=1,i1
             m = m+1
             store = store-c(m)*a(i,l+1)
           end do
           if (equal(a(i,1),zero)) then
              c(i) = zero
           else
              c(i) = store/a(i,1)
           end if
           i = i-1
           jm1 = j
         end do rows

      end if

      end function fpback

      pure function fpbacp(a,b,z,n,k,k1,nest) result(c)
        implicit none

      !  ..scalar arguments..
        integer, intent(in) :: n !! Order of the periodic upper-triangular system.
        integer, intent(in) :: k !! Band width of the main factor block.
        integer, intent(in) :: k1 !! Band width of the wraparound factor block.
        integer, intent(in) :: nest !! Leading storage dimension of the triangular factor.
      !  ..array arguments..
        real(dp), intent(in) :: a(nest,k1) !! Main band block of the periodic triangular factor.
        real(dp), intent(in) :: b(nest,k) !! Wraparound block of the periodic triangular factor.
        real(dp), intent(in) :: z(n) !! Right-hand side of the triangular solve.
      real(dp) :: c(n)
      !  ..local scalars..
      integer :: i,i1,j,l,l0,l1,n2
      real(dp) :: store
      !  ..
      n2 = n-k
      l  = n
      wrapped_rows: do i=1,k
         store = z(l)
         j = k+2-i
         if (i/=1) then
             l0 = l
             do l1=j,k
               l0 = l0+1
               store = store-c(l0)*b(l,l1)
             end do
         endif
         if (equal(b(l,j-1),zero)) then
            c(l) = zero
         else
            c(l) = store/b(l,j-1)
         end if
         l = l-1
         if (l==0) exit wrapped_rows
      end do wrapped_rows

      banded_rows: if (l>0) then

         do i=1,n2
            store = z(i)
            l = n2
            do j=1,k
              l = l+1
              store = store-c(l)*b(i,j)
            end do
            c(i) = store
         end do

         i = n2
         if (equal(a(i,1),zero)) then
            c(i) = zero
         else
            c(i) = c(i)/a(i,1)
         end if

         if (i>1) then
            do j=2,n2
               i = i-1
               store = c(i)
               i1 = k
               if (j<=k) i1=j-1
               l = i
               do l0=1,i1
                 l = l+1
                 store = store-c(l)*a(i,l0+1)
               end do
               if (equal(a(i,1),zero)) then
                  c(i) = zero
               else
                  c(i) = store/a(i,1)
               end if
            end do
         end if

      end if banded_rows

      end function fpbacp

      pure function fpbspl(t,n,k,x,l) result(h)
        implicit none
        integer, intent(in) :: n !! Length of the nondecreasing knot vector.
        integer, intent(in) :: k !! B-spline polynomial degree, at most max_order.
        integer, intent(in) :: l !! Knot interval index containing the query coordinate.
        real(dp), intent(in) :: x !! Query coordinate in the active knot span.
        real(dp), intent(in) :: t(n) !! Nondecreasing B-spline knot vector.
         real(dp) :: h(max_order+1)

         ! local variables
         real(dp) :: f,hh(max_order+1)
         integer :: i,j,li,lj

         h(1) = one
         do j=1,k
           hh(1:j) = h(1:j)
           h(1) = zero
           do i=1,j
             li = l+i
             lj = li-j
             if (not_equal(t(li),t(lj))) then
                f = hh(i)/(t(li)-t(lj))
                h(i)   = h(i)+f*(t(li)-x)
                h(i+1) = f*(x-t(lj))
             else
                h(i+1) = zero
             endif
           end do
         end do

      end function fpbspl

      pure subroutine fpcyt1(a,n,nn)
        implicit none

          !  ..scalar arguments..
        integer, intent(in) :: n !! Order of the cyclic tridiagonal system.
        integer, intent(in) :: nn !! Leading storage dimension of the cyclic factor.
          !  ..array arguments..
        real(dp), intent(inout) :: a(nn,6) !! Cyclic tridiagonal entries, overwritten by factors.
          !  ..local scalars..
          real(dp) :: aa,beta,gamma,sum,teta,v
          integer :: i,n1,n2
          !  ..
          n2   = n-2
          beta = one/a(1,2)
          gamma = a(n,3)
          teta = a(1,1)*beta
          a(1,4) = beta
          a(1,5) = gamma
          a(1,6) = teta
          sum = gamma*teta

          internal_rows: do i=2,n2
             v = a(i-1,3)*beta
             aa = a(i,1)
             beta = one/(a(i,2)-aa*v)
             gamma = -gamma*v
             teta = -teta*aa*beta
             a(i,4) = beta
             a(i,5) = gamma
             a(i,6) = teta
             sum = sum+gamma*teta
          end do internal_rows

          n1 = n-1
          v = a(n2,3)*beta
          aa = a(n1,1)
          beta = one/(a(n1,2)-aa*v)
          gamma = a(n,1)-gamma*v
          teta = (a(n1,3)-teta*aa)*beta
          a(n1,4) = beta
          a(n1,5) = gamma
          a(n1,6) = teta
          a(n,4) = one/(a(n,2)-(sum+gamma*teta))

      return
      end subroutine fpcyt1

      pure subroutine fpcyt2(a,n,b,c,nn)
        implicit none
          !  ..scalar arguments..
        integer, intent(in) :: n !! Order of the factored cyclic system.
        integer, intent(in) :: nn !! Leading storage dimension of the cyclic factor.
          !  ..array arguments..
        real(dp), intent(in) :: a(nn,6) !! Factored cyclic system produced by fpcyt1.
        real(dp), intent(in) :: b(n) !! Right-hand side of the cyclic system.
        real(dp), intent(out) :: c(n) !! Solution of the factored cyclic system.
          !  ..local scalars..
          real(dp) :: cc,sum
          integer :: i,j,j1,n1
          !  ..
          c(1) = b(1)*a(1,4)
          sum  = c(1)*a(1,5)
          n1   = n-1
          do i=2,n1
             c(i) = (b(i)-a(i,1)*c(i-1))*a(i,4)
             sum = sum+c(i)*a(i,5)
          end do
          cc    = (b(n)-sum)*a(n,4)
          c(n)  = cc
          c(n1) = c(n1)-cc*a(n1,6)
          j = n1
          do i=3,n
             j1 = j-1
             c(j1) = c(j1)-c(j)*a(j1,3)*a(j1,4)-cc*a(j1,6)
             j = j1
          end do
          return

      end subroutine fpcyt2

      pure subroutine fpdisc(t,n,k2,b,nest)
        implicit none

      !  ..scalar arguments..
        integer, intent(in) :: n !! Length of the nondecreasing knot vector.
        integer, intent(in) :: k2 !! Polynomial degree plus two.
        integer, intent(in) :: nest !! Leading storage dimension of the triangular factor.
      !  ..array arguments..
        real(dp), intent(in) :: t(n) !! Nondecreasing B-spline knot vector.
        real(dp), intent(inout) :: b(nest,k2) !! Derivative-jump values at the interior knots.
      !  ..local scalars..
      real(dp) :: an,fac,prod
      integer :: i,ik,j,jk,k,k1,l,lj,lk,lmk,lp,nk1,nrint
      !  ..local array..
      real(dp) :: h(2*(degree_5+1))
      !  ..
      k1    = k2-1
      k     = k1-1
      nk1   = n-k1
      nrint = nk1-k
      an    = nrint
      fac   = an/(t(nk1+1)-t(k1))
      do l=k2,nk1
         lmk = l-k1
         do j=1,k1
            ik = j+k1
            lj = l+j
            lk = lj-k2
            h(j)  = t(l)-t(lk)
            h(ik) = t(l)-t(lj)
         end do
         lp = lmk
         do j=1,k2
           jk = j
           prod = h(j)
           do i=1,k
             jk = jk+1
             prod = prod*h(jk)*fac
           end do
           lk = lp+k1
           b(lmk,j) = (t(lk)-t(lp))/prod
           lp = lp+1
         end do
      end do
      return
      end subroutine fpdisc

      elemental subroutine fpgivs(piv,ww,cos,sin)
        implicit none
        real(dp), intent(in) :: piv !! Pivot to eliminate with a Givens rotation.
        real(dp), intent(inout) :: ww !! Diagonal entry updated by the Givens rotation.
        real(dp), intent(out) :: cos !! Cosine of the plane rotation.
        real(dp), intent(out) :: sin !! Sine of the plane rotation.
          !  ..local scalars..
          real(dp) :: dd,store

          store = abs(piv)
          !  evaluate only the taken branch: an if (not merge) avoids a spurious
          !  divide-by-zero in the unselected expression (ww/piv or piv/ww) when
          !  the other operand is zero, which would raise ieee_invalid under fpe traps.
          if (store>=ww) then
             dd = store*sqrt(one+(ww/piv)**2)
          else
             dd = ww   *sqrt(one+(piv/ww)**2)
          end if
          cos = ww/dd
          sin = piv/dd
          ww  = dd
          return
      end subroutine fpgivs

      pure subroutine fpgrsp(ifsu,ifsv,ifbu,ifbv,lback,u,mu,v,                         &
                             mv,r,mr,dr,iop0,iop1,tu,nu,tv,nv,p,c,nc,sq,fp,fpu,fpv,mm, &
                             mvnu,spu,spv,right,q,au,av1,av2,bu,bv,a0,a1,b0,b1,c0,c1,  &
                             cosi,nru,nrv)
        implicit none
      !  ..scalar arguments..
        real(dp), intent(in) :: p !! Smoothing parameter; negative selects unpenalized fitting.
        real(dp), intent(inout) :: fp !! Sum of squared grid residuals.
        real(dp), intent(inout) :: &
            sq !! Squared residual objective including eliminated constraints.
        logical, intent(in) :: lback !! Reuse cached factors and back-substitute only.
        integer, intent(in) :: mu !! Number of colatitude sites.
        integer, intent(in) :: mv !! Number of longitude sites.
        integer, intent(in) :: mr !! Length of the longitude-fastest observation vector.
        integer, intent(in) :: iop0 !! North-pole first-derivative continuity flag.
        integer, intent(in) :: iop1 !! South-pole first-derivative continuity flag.
        integer, intent(in) :: nu !! Active colatitude knot count.
        integer, intent(in) :: nv !! Active longitude knot count.
        integer, intent(in) :: nc !! Available tensor coefficient storage.
        integer, intent(in) :: mm !! Capacity of the row/right-hand-side workspace.
        integer, intent(in) :: mvnu !! Capacity of the transformed observation workspace.
        integer, intent(inout) :: ifsu !! Colatitude observation-factor cache flag.
        integer, intent(inout) :: ifsv !! Longitude observation-factor cache flag.
        integer, intent(inout) :: ifbu !! Colatitude penalty-factor cache flag.
        integer, intent(inout) :: ifbv !! Longitude penalty-factor cache flag.
      !  ..array arguments..
        real(dp), intent(in) :: u(mu) !! Interior colatitude sites in radians.
        real(dp), intent(in) :: v(mv) !! Periodic longitude sites in radians.
        real(dp), intent(in) :: r(mr) !! Longitude-fastest scaled grid observations.
        real(dp), intent(in) :: dr(6) !! Pole values and two first-derivative components per pole.
        real(dp), intent(in) :: tu(nu) !! Full colatitude knot vector.
        real(dp), intent(in) :: tv(nv) !! Full periodic longitude knot vector.

        real(dp), intent(inout) :: av1(nv,6) !! First band block of the periodic longitude factor.
        real(dp), intent(inout) :: av2(nv,4) !! Wraparound block of the periodic longitude factor.
        real(dp), intent(inout) :: cosi(nv,2) !! Periodic spline coefficients for cosine and sine.
        real(dp), intent(inout) :: spu(mu,4) !! Four active colatitude basis values per site.
        real(dp), intent(inout) :: spv(mv,4) !! Four active longitude basis values per site.
        real(dp), intent(inout) :: right(mm) !! Right-hand-side row transformed during rotations.
        real(dp), intent(inout) :: bu(nu,5) !! Colatitude derivative-jump penalty rows.
        real(dp), intent(inout) :: bv(nv,5) !! Longitude derivative-jump penalty rows.
        real(dp), intent(inout) :: a0(2,mv) !! North-pole observation constraint contributions.
        real(dp), intent(inout) :: b0(2,nv) !! North-pole derivative-jump constraint contributions.
        real(dp), intent(inout) :: c0(nv) !! North-pole adjacent coefficient row.
        real(dp), intent(inout) :: c1(nv) !! South-pole adjacent coefficient row.
        real(dp), intent(inout) :: a1(2,mv) !! South-pole observation constraint contributions.
        real(dp), intent(inout) :: b1(2,nv) !! South-pole derivative-jump constraint contributions.
        real(dp), intent(inout) :: q(mvnu) !! Transformed observation block workspace.
        real(dp), intent(inout) :: au(nu,5) !! Banded upper-triangular colatitude factor.
        real(dp), intent(inout) :: c(nc) !! Longitude-fastest output tensor B-spline coefficients.
        real(dp), intent(inout) :: fpu(nu) !! Residual contributions by colatitude interval.
        real(dp), intent(inout) :: fpv(nv) !! Residual contributions by longitude interval.
        integer, intent(inout) :: nru(mu) !! Colatitude interval indices, zero based.
        integer, intent(inout) :: nrv(mv) !! Longitude interval indices, zero based.
      !  ..local scalars..
      real(dp) :: arg,dr01,dr11,fac,fac0,fac1,pinv,term
      integer :: i,ii,ij,ik,irot,it,ir,i0,i1,i2,j,jj,jk,jper, &
                 j0,j1,k,k1,l,l0,l1,mvv,ncof,nrold,nroldu,nroldv,number, &
                 numu,numu1,numv,numv1,nuu,nu4,nu7,nu8,nu9,nv11,nv4,nv7,nv8,n1
      !  ..local arrays..
      real(dp) :: h(max_order+1),h1(degree_3+2),h2(degree_3+1)

      !  let
      !               |     (spu)      |            |     (spv)      |
      !        (au) = | -------------- |     (av) = | -------------- |
      !               | sqrt(1/p) (bu) |            | sqrt(1/p) (bv) |
      !
      !                                | r  ' 0 |
      !                            q = | ------ |
      !                                | 0  ' 0 |
      !
      !  with c      : the (nu-4) x (nv-4) matrix which contains the b-spline
      !                coefficients.
      !       r      : the mu x mv matrix which contains the function values.
      !       spu,spv: the mu x (nu-4), resp. mv x (nv-4) observation matrices
      !                according to the least-squares problems in the u-,resp.
      !                v-direction.
      !       bu,bv  : the (nu-7) x (nu-4),resp. (nv-7) x (nv-4) matrices
      !                containing the discontinuity jumps of the derivatives
      !                of the b-splines in the u-,resp.v-variable at the knots
      !  the b-spline coefficients of the smoothing spline are then calculated
      !  as the least-squares solution of the following over-determined linear
      !  system of equations
      !
      !  (1)  (av) c (au)' = q
      !
      !  subject to the constraints
      !
      !  (2)  c(i,nv-3+j) = c(i,j), j=1,2,3 ; i=1,2,...,nu-4
      !
      !  (3)  if iop0 = 0  c(1,j) = dr(1)
      !          iop0 = 1  c(1,j) = dr(1)
      !                    c(2,j) = dr(1)+(dr(2)*cosi(1,j)+dr(3)*cosi(2,j))*
      !                            tu(5)/3. = c0(j) , j=1,2,...nv-4
      !
      !  (4)  if iop1 = 0  c(nu-4,j) = dr(4)
      !          iop1 = 1  c(nu-4,j) = dr(4)
      !                    c(nu-5,j) = dr(4)+(dr(5)*cosi(1,j)+dr(6)*cosi(2,j))
      !                                *(tu(nu-4)-tu(nu-3))/3. = c1(j)
      !

      !  initialization
      nu4 = nu-4
      nu7 = nu-7
      nu8 = nu-8
      nu9 = nu-9
      nv4 = nv-4
      nv7 = nv-7
      nv8 = nv-8
      nv11 = nv-11
      nuu  = nu4-iop0-iop1-2
      pinv = merge(one/p,one,p>zero)

      !  it depends on the value of the flags ifsu,ifsv,ifbu,ifbv,iop0,iop1
      !  and on the value of p whether the matrices (spu), (spv), (bu), (bv),
      !  (cosi) still must be determined.
      if (ifsu==0) then
          !  calculate the non-zero elements of the matrix (spu) which is the observation matrix
          !  according to the least-squares spline approximation problem in the u-direction.
          l  = 4
          l1 = 5
          number = 0
          do it=1,mu
            arg = u(it)
            l0 = l
            l = fp_knot_interval(tu, arg, l, nu4)
            number = number + (l - l0)
            l1 = l + 1
            h = fpbspl(tu,nu,degree_3,arg,l)
            spu(it,1:4) = h(1:4)
            nru(it) = number
          end do
          ifsu = 1
      endif

      !  calculate the non-zero elements of the matrix (spv) which is the ob-
      !  servation matrix according to the least-squares spline approximation
      !  problem in the v-direction.
      if (ifsv==0) then

          l  = 4
          l1 = 5
          number = 0
          do it=1,mv
            arg = v(it)
            l0 = l
            l = fp_knot_interval(tv, arg, l, nv4)
            number = number + (l - l0)
            l1 = l + 1
            h = fpbspl(tv,nv,degree_3,arg,l)
            spv(it,1:4) = h(1:4)
            nrv(it) = number
          end do
          ifsv = 1

          !  calculate the coefficients of the interpolating splines for cos(v) and sin(v).
          if (iop0/=0 .or. iop1/=0) then

              cosi(1:nv4,:) = zero
              if (nv7>=4) then
                  do i=1,nv7
                     l = i+3
                     arg = tv(l)
                     h = fpbspl(tv,nv,degree_3,arg,l)
                     av1(i,1:3) = h(1:3)
                     cosi(i,1) = cos(arg)
                     cosi(i,2) = sin(arg)
                  end do
                  call fpcyt1(av1,nv7,nv)
                  do j=1,2
                     call fpcyt2(av1,nv7,cosi(:nv7,j),right,nv)
                     cosi(:nv4,j) = [right(nv7),right(1:nv7),right(1:2)]
                  end do
              endif
          endif

      endif

      if (p>zero) then
          !  calculate the non-zero elements of the matrix (bu).
          if (ifbu==0 .and. nu8/=0) then
             call fpdisc(tu,nu,degree_5,bu,nu)
             ifbu = 1
          endif

          !  calculate the non-zero elements of the matrix (bv).
          if (ifbv==0 .and. nv8/=0) then
             call fpdisc(tv,nv,degree_5,bv,nv)
             ifbv = 1
          endif
      endif
      !  substituting (2),(3) and (4) into (1), we obtain the overdetermined system
      !         (5)  (avv) (cc) (auu)' = (qq)
      !  from which the nuu*nv7 remaining coefficients
      !         c(i,j) , i=2+iop0,3+iop0,...,nu-5-iop1,j=1,2,...,nv-7.
      !  the elements of (cc), are then determined in the least-squares sense.
      !  simultaneously, we compute the resulting sum of squared residuals sq.
      dr01 = dr(1)
      dr11 = dr(4)
      a0(1,1:mv) = dr01
      a1(1,1:mv) = dr11
      if (nv8/=0 .and. p>zero) then
         b0(1,1:nv8) = zero
         b1(1,1:nv8) = zero
      endif
      mvv = mv
      if (iop0/=0) then
          fac       = (tu(5)-tu(4))*third
          c0(1:nv4) = dr01+fac*matmul(cosi(:nv4,:),dr(2:3))

          do i = 1, mv
              a0(2,i) = dot_product(spv(i,1:4),c0(nrv(i)+1:nrv(i)+4))
          end do

          if (nv8/=0 .and. p>zero) then
              do i = 1, nv8
                  b0(2,i) = pinv*dot_product(bv(i,1:5),c0(i:i+4))
              end do
              mvv = mv+nv8
          endif
      endif
      if (iop1/=0) then
          fac       = (tu(nu4)-tu(nu4+1))*third
          c1(1:nv4) = dr11 + fac*matmul(cosi(:nv4,:),dr(5:6))

          do i = 1, mv
              a1(2,i) = dot_product(c1(nrv(i)+1:nrv(i)+4),spv(i,1:4))
          end do

          if (nv8/=0 .and. p>zero) then
              do i = 1, nv8
                  b1(2,i) = pinv*dot_product(bv(i,1:5),c1(i:i+4))
              end do
              mvv = mv+nv8
          endif
      endif

      !  we first determine the matrices (auu) and (qq). then we reduce the matrix (auu) to an unit
      !  upper triangular form (ru) using givens rotations without square roots. we apply the same
      !  transformations to the rows of matrix qq to obtain the mv x nuu matrix g.
      !  we store matrix (ru) into au and g into q.
      l = mvv*nuu
      !  initialization.
      sq = zero

      if (l/=0) then
         q(1:l) = zero
         au(1:nuu,1:5) = zero
         l = 0
      endif
      nrold = 0
      n1 = nrold+1
      auu_iterations: do it=1,mu
         number = nru(it)

         auu_inner: do

         ! find the appropriate column of q.
         right(1:mvv) = zero

         if (nrold==number) then
            !  fetch a new row of matrix (spu).
            h(1:4) = spu(it,1:4)

            !  find the appropriate column of q.
            right(1:mv) = r(l+1:l+mv)
            l  = l+mv
            i0 = 1
            i1 = 4
         else

            if (p<=zero) then
               nrold = n1
               n1 = n1+1
               cycle auu_inner
            end if
            !  fetch a new row of matrix (bu).
            h(1:5) = bu(n1,1:5)*pinv
            i0 = 1
            i1 = 5
         end if

         j0 = n1
         j1 = nu7-number

         ! take into account that we eliminate the constraints (3)
         do while (j0-1<=iop0)
             fac0 = h(i0)
             right(1:mv) = right(1:mv)-fac0*a0(j0,1:mv)
             if (mv/=mvv) then
                 j = mv
                 do jj=1,nv8
                    j = j+1
                    right(j) = right(j)-fac0*b0(j0,jj)
                 end do
             endif
             j0 = j0+1
             i0 = i0+1
         end do

         ! take into account that we eliminate the constraints (4)
         do while (j1-1<=iop1)
            fac1 = h(i1)
            right(1:mv) = right(1:mv)-fac1*a1(j1,1:mv)
            if (mv/=mvv) then
                j = mv
                do jj=1,nv8
                   j = j+1
                   right(j) = right(j)-fac1*b1(j1,jj)
                end do
            endif
            j1 = j1+1
            i1 = i1-1
         end do

         irot = max(0,nrold-iop0-1)

         ! rotate the new row of matrix (auu) into triangle.
         if (i0<=i1) call fp_rotate_row_block(h(i0), i1-i0+1, au, nu, right, q, mvv, irot)

         ! we update the sum of squared residuals.
         sq = sq+sum(right(:mvv)**2)
         if (nrold==number) exit auu_inner
         nrold = n1
         n1 = n1+1
         end do auu_inner
      end do auu_iterations

      if (nuu/=0) then
! we determine the matrix (avv) and then we reduce her to an unit upper triangular form (rv)
          !  using givens rotations without square roots. we apply the same transformations to the
! columns of matrix g to obtain the (nv-7) x (nu-6-iop0-iop1) matrix h. we store matrix (rv)
! into av1 and av2, h into c. the nv7 x nv7 triangular unit upper matrix (rv) has the form
          !              | av1 '     |
          !       (rv) = |     ' av2 |
          !              |  0  '     |
          !  with (av2) a nv7 x 4 matrix and (av1) a nv11 x nv11 unit upper triangular matrix of
          !  bandwidth 5.
          ncof = nuu*nv7
          !  initialization.
          c(:ncof) = zero
          av1(1:nv4,1:5) = zero
          av2(1:nv4,1:4) = zero
          jper = 0
          nrold = 0
          avv_iterations: do it=1,mv
             number = nrv(it)

             avv_inner: do

             if (nrold==number) then

               ! fetch a new row of matrix (spv)
               h(1:5) = [spv(it,1:4),zero]

               ! find the appropriate row of g.
               right(1:nuu) = q(it:it+(nuu-1)*mvv:mvv)

            else

               if (p<=zero) then
                  nrold = nrold+1
                  cycle avv_inner
               endif
               ! fetch a new row of matrix (bv).
               n1 = nrold+1
               h(1:5) = bv(n1,1:5)*pinv

               ! find the appropriate row of g.
               if(mv/=mvv) then
                  l = mv+n1
                  right(1:nuu) = q(l:l+(nuu-1)*mvv:mvv)
               else
                  right(1:nuu) = zero
               endif

            endif

            ! test whether there are non-zero values in the new row of (avv)
            ! corresponding to the b-splines n(j;v),j=nv7+1,...,nv4.
            if (nrold<nv11) then

               ! rotation into triangle of the new row of (avv), in case the elements
               ! corresponding to the b-splines n(j;v),j=nv7+1,...,nv4 are all zero.
               call fp_rotate_row_stride(h, 5, av1, nv, right, c, nv7, nuu, nrold)

            else

               ! initialize the matrix (av2).
               if (jper==0) then

                  jk = nv11+1
                  do i=1,4
                     ik = jk
                     do j=1,5
                         if (ik<=0) exit
                         av2(ik,i) = av1(ik,j)
                         ik = ik-1
                     end do
                     jk = jk+1
                  end do
                  jper = 1

               endif

               ! if one of the non-zero elements of the new row corresponds to one of the b-splines
               ! n(j;v),j=nv7+1,...,nv4, we take account of condition (2) for setting up this row of
               ! (avv). the row is stored in h1(the part with respect to av1) and h2 (the part with
               ! respect to av2).
               h1 = zero
               h2 = zero
               j = nrold-nv11
               do i=1,5
                   j = j+1
                   l0 = j
                   l1 = l0-4
                   do while (l1>max(0,nv11))
                      l0 = l1-nv11
                      l1 = l0-4
                   end do
                   if (l1<=0) then
                      h2(l0) = h2(l0) + h(i)
                   else ! l1<=nv11
                      h1(l1) = h(i)
                   end if
               end do

               ! rotate the new row of (avv) into triangle.
               call fp_rotate_2mat_stride(h1, 5, h2, 4, av1, av2, nv, &
                                          right, c, nv7, nuu, 1, nv11)

            end if

            ! we update the sum of squared residuals.
            sq = sq+sum(right(:nuu)**2)

            if (nrold==number) exit avv_inner

            nrold = nrold+1
            end do avv_inner

          end do avv_iterations

          !  test whether the b-spline coefficients must be determined.
          if (lback) return
          !  backward substitution to obtain the b-spline coefficients as the
          !  solution of the linear system    (rv) (cr) (ru)' = h.
          !  first step: solve the system  (rv) (c1) = h.
          k = 1
          do i=1,nuu
             c(k:k+nv7-1) = fpbacp(av1,av2,c(k),nv7,degree_4,degree_5,nv)
             k = k+nv7
          end do

          !  second step: solve the system  (cr) (ru)' = (c1).
          do k=1,nv7
            right(:nuu) = c(k:k+(nuu-1)*nv7:nv7)
            c(k:k+(nuu-1)*nv7:nv7) = fpback(au,right(:nuu),nuu,degree_5,nu)
          end do

      endif

      !  calculate from the conditions (2)-(3)-(4), the remaining b-spline
      !  coefficients.
      ncof = nu4*nv4

      q(1:nv4) = dr01
      q(ncof-nv4+1:ncof) = dr11

      i = nv4
      j = 0
      if (iop0/=0) then
         q(i+1:i+nv4) = c0(1:nv4)
         i = i+nv4
      endif
      if(nuu/=0) then
         do l=1,nuu
            ii = i
            do k=1,nv7
               i = i+1
               j = j+1
               q(i) = c(j)
            end do
            do k=1,3
               ii = ii+1
               i  = i+1
               q(i) = q(ii)
            end do
         end do
      endif
      if(iop1/=0) then
          do l=1,nv4
             i = i+1
             q(i) = c1(l)
          end do
      endif
      c(1:ncof) = q(1:ncof)
      !  calculate the quantities
      !    res(i,j) = (r(i,j) - s(u(i),v(j)))**2 , i=1,2,..,mu;j=1,2,..,mv
      !    fp = sumi=1,mu(sumj=1,mv(res(i,j)))
      !    fpu(r) = sum''i(sumj=1,mv(res(i,j))) , r=1,2,...,nu-7
      !                  tu(r+3) <= u(i) <= tu(r+4)
      !    fpv(r) = sumi=1,mu(sum''j(res(i,j))) , r=1,2,...,nv-7
      !                  tv(r+3) <= v(j) <= tv(r+4)
      fp  = zero
      fpu = zero
      fpv = zero
      ir = 0
      nroldu = 0

      !  main loop for the different grid points.
      grid_u: do i1=1,mu
         numu = nru(i1)
         numu1 = numu+1
         nroldv = 0
         grid_v: do i2=1,mv
            numv = nrv(i2)
            numv1 = numv+1
            ir = ir+1
            ! evaluate s(u,v) at the current grid point by making the sum of the
            ! cross products of the non-zero b-splines at (u,v), multiplied with
            ! the appropriate b-spline coefficients.
            term = zero
            k1 = numu*nv4+numv
            do l1 = 1,4
               fac  = spu(i1,l1)
               term = term+fac*dot_product(spv(i2,1:4),c(k1+1:k1+4))
               k1   = k1+nv4
            end do
            ! calculate the squared residual at the current grid point.
            term = (r(ir)-term)**2
            ! adjust the different parameters.
            fp = fp+term
            fpu(numu1) = fpu(numu1)+term
            fpv(numv1) = fpv(numv1)+term
            fac = term*half
            if (numv/=nroldv) then
               fpv(numv1) = fpv(numv1)-fac
               fpv(numv) = fpv(numv)+fac
            endif
            nroldv = numv
            if (numu/=nroldu) then
               fpu(numu1) = fpu(numu1)-fac
               fpu(numu) = fpu(numu)+fac
            endif
         end do grid_v
         nroldu = numu
      end do grid_u
      return
      end subroutine fpgrsp

      pure subroutine fpopsp(ifsu,ifsv,ifbu,ifbv,u,mu,v,mv,r, &
                             mr,r0,r1,dr,iopt,ider,tu,nu,tv,nv,nuest,nvest,p,step,c,nc, &
                             fp,fpu,fpv,nru,nrv,wrk,lwrk)
        implicit none
      !
      !  ..scalar arguments..
        integer, intent(in) :: mu !! Number of colatitude sites.
        integer, intent(in) :: mv !! Number of longitude sites.
        integer, intent(in) :: mr !! Length of the longitude-fastest observation vector.
        integer, intent(in) :: nu !! Active colatitude knot count.
        integer, intent(in) :: nv !! Active longitude knot count.
        integer, intent(in) :: nuest !! Allocated colatitude knot capacity.
        integer, intent(in) :: nvest !! Allocated longitude knot capacity.
        integer, intent(in) :: nc !! Available tensor coefficient storage.
        integer, intent(in) :: lwrk !! Available work entries.
        integer, intent(inout) :: ifsu !! Colatitude observation-factor cache flag.
        integer, intent(inout) :: ifsv !! Longitude observation-factor cache flag.
        integer, intent(inout) :: ifbu !! Colatitude penalty-factor cache flag.
        integer, intent(inout) :: ifbv !! Longitude penalty-factor cache flag.
        real(dp), intent(in) :: r0 !! Supplied scaled north-pole value.
        real(dp), intent(in) :: r1 !! Supplied scaled south-pole value.
        real(dp), intent(inout) :: fp !! Sum of squared grid residuals.
        real(dp), intent(inout) :: p !! Smoothing parameter; negative selects unpenalized fitting.
      !  ..array arguments..
        integer, intent(in) :: ider(4) !! Pole value exactness and zero-slope flags.
        integer, intent(in) :: iopt(3) !! Fit mode and north/south continuity flags.
        integer, intent(inout) :: nru(mu) !! Colatitude interval indices, zero based.
        integer, intent(inout) :: nrv(mv) !! Longitude interval indices, zero based.
        real(dp), intent(in) :: u(mu) !! Interior colatitude sites in radians.
        real(dp), intent(in) :: v(mv) !! Periodic longitude sites in radians.
        real(dp), intent(in) :: r(mr) !! Longitude-fastest scaled grid observations.
        real(dp), intent(in) :: tu(nu) !! Full colatitude knot vector.
        real(dp), intent(in) :: tv(nv) !! Full periodic longitude knot vector.
        real(dp), intent(inout) :: &
            dr(6) !! Pole values and two first-derivative components per pole.
        real(dp), intent(inout) :: wrk(lwrk) !! Caller-owned factorization and constraint workspace.
        real(dp), intent(inout) :: step(2) !! Finite-difference optimization increments per pole.
        real(dp), intent(inout) :: c(nc) !! Longitude-fastest output tensor B-spline coefficients.
        real(dp), intent(inout) :: fpu(nu) !! Residual contributions by colatitude interval.
        real(dp), intent(inout) :: fpv(nv) !! Residual contributions by longitude interval.
      !  ..local scalars..
      real(dp) :: sq,sqq,sq0,sq1,step1,step2
      integer :: i,id0,iop0,iop1,i1,j,l,lau,lav1,lav2,la0,la1,lbu,lbv,lb0, &
       lb1,lc0,lc1,lcs,lq,lri,lsu,lsv,l1,l2,mm,mvnu,number, id1
      !  ..local arrays..
      integer :: nr(6)
      logical :: zeroed
      real(dp) :: delta(6),drr(6),sum(6),a(6,6),g(6)

      !  we partition the working space
      lsu  = 1
      lsv  = lsu+4*mu
      lri  = lsv+4*mv
      mm   = max(nuest,mv+nvest)
      lq   = lri+mm
      mvnu = nuest*(mv+nvest-8)
      lau  = lq+mvnu
      lav1 = lau+5*nuest
      lav2 = lav1+6*nvest
      lbu  = lav2+4*nvest
      lbv  = lbu+5*nuest
      la0  = lbv+5*nvest
      la1  = la0+2*mv
      lb0  = la1+2*mv
      lb1  = lb0+2*nvest
      lc0  = lb1+2*nvest
      lc1  = lc0+nvest
      lcs  = lc1+nvest
      !  we calculate the smoothing spline sp(u,v) according to the input
      !  values dr(i),i=1,...,6.
      iop0 = iopt(2)
      iop1 = iopt(3)
      id0  = ider(1)
      id1  = ider(3)
      sq   = zero
      sum  = zero

      call fpgrsp(ifsu,ifsv,ifbu,ifbv,fp_false,u,mu,v,mv,r,mr,dr,           &
                  iop0,iop1,tu,nu,tv,nv,p,c,nc,sq,fp,fpu,fpv,mm,mvnu,       &
                  wrk(lsu),wrk(lsv),wrk(lri),wrk(lq),wrk(lau),wrk(lav1),   &
                  wrk(lav2),wrk(lbu),wrk(lbv),wrk(la0),wrk(la1),wrk(lb0),  &
                  wrk(lb1),wrk(lc0),wrk(lc1),wrk(lcs),nru,nrv)

      sq0 = merge((r0-dr(1))**2,zero,id0==0)
      sq1 = merge((r1-dr(4))**2,zero,id1==0)
      sq  = sq+sq0+sq1

      ! in case all derivative values dr(i) are given (step<=0) or in case
      ! we have spline interpolation, we accept this spline as a solution.
      if (sq<=zero .or. all(step(1:2)<=zero)) return
      drr    = dr
      sqq    = sq
      sum(1) = sq

      ! number denotes the number of derivative values dr(i) that still must
      ! be optimized. let us denote these parameters by g(j),j=1,...,number.
      number = 0
      if (id0<=0) then
         number   = 1
         nr   (1) = 1
         delta(1) = step(1)
      endif
      if (iop0/=0 .and. ider(2)==0) then
         step2 = step(1)*three/(tu(5)-tu(4))
         nr   (number+1:number+2) = [2,3]
         delta(number+1:number+2) = step2
         number = number+2
      endif
      if (id1<=0) then
         number = number+1
         nr   (number) = 4
         delta(number) = step(2)
      end if
      if (iop1/=0 .and. ider(4)==0) then
         step2 = step(2)*three/(tu(nu)-tu(nu-4))
         nr   (number+1:number+2) = [5,6]
         delta(number+1:number+2) = step2
         number = number+2
      endif

      if (number==0) return

      ! the sum of squared residulas sq is a quadratic polynomial in the parameters g(j).
! we determine the unknown coefficients of this polymomial by calculating
! (number+1)*(number+2)/2
      ! different splines according to specific values for g(j).
      zeroed = .false.
      do i=1,number
         l      = nr(i)
         step1  = delta(i)
         drr(l) = dr(l)+step1
         call fpgrsp(ifsu,ifsv,ifbu,ifbv,fp_true,u,mu,v,mv,r,mr,drr, &
                     iop0,iop1,tu,nu,tv,nv,p,c,nc,sum(i),fp,fpu,fpv,mm,mvnu, &
                     wrk(lsu),wrk(lsv),wrk(lri),wrk(lq),wrk(lau),wrk(lav1), &
                     wrk(lav2),wrk(lbu),wrk(lbv),wrk(la0),wrk(la1),wrk(lb0), &
                     wrk(lb1),wrk(lc0),wrk(lc1),wrk(lcs),nru,nrv)

         if (id0==0) sq0 = (r0-drr(1))**2
         if (id1==0) sq1 = (r1-drr(4))**2
         sum(i) = sum(i)+sq0+sq1
         drr(l) = dr(l)-step1

         call fpgrsp(ifsu,ifsv,ifbu,ifbv,fp_true,u,mu,v,mv,r,mr,drr, &
                     iop0,iop1,tu,nu,tv,nv,p,c,nc,sqq,fp,fpu,fpv,mm,mvnu, &
                     wrk(lsu),wrk(lsv),wrk(lri),wrk(lq),wrk(lau),wrk(lav1), &
                     wrk(lav2),wrk(lbu),wrk(lbv),wrk(la0),wrk(la1),wrk(lb0), &
                     wrk(lb1),wrk(lc0),wrk(lc1),wrk(lcs),nru,nrv)

         if (id0==0) sq0 = (r0-drr(1))**2
         if (id1==0) sq1 = (r1-drr(4))**2

         sqq    = sqq+sq0+sq1
         drr(l) = dr(l)
         a(i,i) = (sum(i)+sqq-sq-sq)/step1**2
         if (a(i,i)<=zero) then
             zeroed = .true.
             exit
         end if
         g(i)   = (sqq-sum(i))/(2*step1)
      end do

      if (.not.zeroed) then

          if (number>1) then
             do i=2,number
                l1      = nr(i)
                step1   = delta(i)
                drr(l1) = dr(l1)+step1
                i1 = i-1
                do j=1,i1
                    l2      = nr(j)
                    step2   = delta(j)
                    drr(l2) = dr(l2)+step2
                    call fpgrsp(ifsu,ifsv,ifbu,ifbv,fp_true,u,mu,v,mv,r,mr,drr, &
                                iop0,iop1,tu,nu,tv,nv,p,c,nc,sqq,fp,fpu,fpv,mm,mvnu, &
                                wrk(lsu),wrk(lsv),wrk(lri),wrk(lq),wrk(lau),wrk(lav1), &
                                wrk(lav2),wrk(lbu),wrk(lbv),wrk(la0),wrk(la1),wrk(lb0), &
                                wrk(lb1),wrk(lc0),wrk(lc1),wrk(lcs),nru,nrv)
                    if(id0==0) sq0 = (r0-drr(1))**2
                    if(id1==0) sq1 = (r1-drr(4))**2
                    sqq = sqq+sq0+sq1
                    a(i,j) = (sq+sqq-sum(i)-sum(j))/(step1*step2)
                    drr(l2) = dr(l2)
                end do
                drr(l1) = dr(l1)
             end do
          end if

          ! the optimal values g(j) are found as the solution of the system
          ! d (sq) / d (g(j)) = 0 , j=1,...,number.
          call fpsysy(a,number,g)
          do i = 1, number
              dr(nr(i))=dr(nr(i))+g(i)
          end do

      endif

      ! we determine the spline sp(u,v) according to the optimal values g(j).
      call fpgrsp(ifsu,ifsv,ifbu,ifbv,fp_false,u,mu,v,mv,r,mr,dr, &
                  iop0,iop1,tu,nu,tv,nv,p,c,nc,sq,fp,fpu,fpv,mm,mvnu, &
                  wrk(lsu),wrk(lsv),wrk(lri),wrk(lq),wrk(lau),wrk(lav1), &
                  wrk(lav2),wrk(lbu),wrk(lbv),wrk(la0),wrk(la1),wrk(lb0), &
                  wrk(lb1),wrk(lc0),wrk(lc1),wrk(lcs),nru,nrv)
      if(id0==0) sq0 = (r0-dr(1))**2
      if(id1==0) sq1 = (r1-dr(4))**2
      sq = sq+sq0+sq1
      return
      end subroutine fpopsp

      elemental subroutine fprota(cos,sin,a,b)
        implicit none

          !  ..scalar arguments..
        real(dp), intent(in) :: cos !! Cosine of the plane rotation.
        real(dp), intent(in) :: sin !! Sine of the plane rotation.
        real(dp), intent(inout) :: a !! First component, overwritten by its rotation.
        real(dp), intent(inout) :: b !! Second component, overwritten by its rotation.

          ! ..local scalars..
          real(dp) :: stor1,stor2

          !  ..
          stor1 = a
          stor2 = b
          b = cos*stor2+sin*stor1
          a = cos*stor1-sin*stor2
          return

      end subroutine fprota

      pure subroutine fpsysy(a,n,g)
        implicit none

      !  ..scalar arguments..
        integer, intent(in) :: n !! Active order of the symmetric system, at most six.
      !  ..array arguments..
        real(dp), intent(inout) :: a(6,6) !! Symmetric matrix, overwritten by LDL factors.
        real(dp), intent(inout) :: &
            g(6) !! Right-hand side, overwritten by the symmetric-system solution.

      !  ..local scalars..
      real(dp) :: fac
      integer :: i,i1,j,k
      !  ..
      g(1) = g(1)/a(1,1)
      if (n==1) return
! decomposition of the symmetric matrix (a) = (l) * (d) *(l)'  with (l) a unit lower triangular
      !  matrix and (d) a diagonal matrix
      a(2:n,1) = a(2:n,1)/a(1,1)

      do i=2,n
         i1 = i-1
         do k=i,n
            fac = a(k,i)
            do j=1,i1
               fac = fac-a(j,j)*a(k,j)*a(i,j)
            end do
            a(k,i) = fac
            if (k>i) a(k,i) = fac/a(i,i)
         end do
      end do
      !  solve the system (l)*(d)*(l)'*(b) = (g).

      !  first step : solve (l)*(d)*(c) = (g).
      do i=2,n
         i1 = i-1
         fac = g(i)
         do j=1,i1
            fac = fac-g(j)*a(j,j)*a(i,j)
         end do
         g(i) = fac/a(i,i)
      end do

      !  second step : solve (l)'*(b) = (c)
      i = n
      do j=2,n
         i1 = i
         i = i-1
         fac = g(i)
         do k=i1,n
            fac = fac-g(k)*a(k,i)
         end do
         g(i) = fac
      end do
      return
      end subroutine fpsysy

    elemental logical function not_equal(a,b) result(ok)
        implicit none
        real(dp), intent(in) :: a !! First value in the spacing-based comparison.
        real(dp), intent(in) :: b !! Second value in the spacing-based comparison.
       ok = .not.equal(a,b)
    end function not_equal
end module scifort_interpolate_sphere_grid_core
