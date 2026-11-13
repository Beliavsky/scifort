! SPDX-License-Identifier: BSD-3-Clause
! Copyright (c) 2022, Federico Perini. All rights reserved.
! Copyright (c) 2026 SciFort contributors
!
! Adapted from perazz/fitpack src/fitpack_core.F90, commit
! 1907c7bad9f5c4557970e906bdd4909eb52b2d4c, by Federico Perini,
! based on Paul Dierckx's FITPACK. Retained license:
! LICENSES/perazz-fitpack-BSD-3-Clause.txt.
! Rectangular grid interpolation and adaptive smoothing paths are included.
module scifort_interpolate_sphere_grid_core
    use scifort_kinds, only : dp
    use scifort_constants, only : pi => scifort_pi
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    implicit none
    private
    public :: fitpack_sphere_grid
    ! Internal numerical helpers shared with the scattered spherical fitter.
    public :: equal,not_equal,fp_knot_interval,fpback,fpbspl,fpdisc,fpgivs,fprota
    public :: root_finding_iterate
    integer, parameter :: max_order = 19, degree_3 = 3, degree_4 = 4, degree_5 = 5
    real(dp), parameter :: zero = 0.0_dp, one = 1.0_dp, two = 2.0_dp
    real(dp), parameter :: three = 3.0_dp, half = 0.5_dp, pi2 = 2.0_dp*pi
    real(dp), parameter :: third = 1.0_dp/3.0_dp
    logical, parameter :: fp_true = .true., fp_false = .false.
    integer, parameter :: fitpack_ok = 0, fitpack_interpolating_ok = -1
    integer, parameter :: fitpack_leastsquares_ok = -2, fitpack_insufficient_storage = 1
    integer, parameter :: fitpack_s_too_small = 2, fitpack_maxit = 3
    integer, parameter :: knot_dim_none = 0, knot_dim_1 = 1, knot_dim_2 = 2, ione = 1
contains
    subroutine fitpack_sphere_grid(u,v,r,continuity,poles,given,exact,flat,s,tu,tv,c,fp,status)
        implicit none
        real(dp), intent(in) :: u(:) !! Increasing interior colatitude grid, at least four sites.
        real(dp), intent(in) :: v(:) !! Increasing periodic longitude grid, at least four sites.
        real(dp), intent(in) :: r(:,:) !! Finite scaled samples, shape (size(u),size(v)).
        logical, intent(in) :: continuity(2) !! First derivative continuity at north/south poles.
        real(dp), intent(in) :: poles(2) !! Scaled supplied pole values, ignored unless given.
        logical, intent(in) :: given(2) !! Whether each pole value was supplied.
        logical, intent(in) :: exact(2) !! Exactness of supplied pole values.
        logical, intent(in) :: flat(2) !! Zero first derivatives at north/south poles.
        real(dp), intent(in) :: s !! Nonnegative squared-residual target in scaled data units.
        real(dp), allocatable, intent(out) :: tu(:) !! Full active colatitude knots.
        real(dp), allocatable, intent(out) :: tv(:) !! Full active longitude knots.
        real(dp), allocatable, intent(out) :: c(:) !! Longitude-fastest tensor coefficients.
        real(dp), intent(out) :: fp !! Squared grid residual in scaled data units.
        integer, intent(out) :: status !! 0 success, 3 arithmetic failure; caller validates inputs.
        real(dp), allocatable :: xknots(:),yknots(:),coef(:),data(:),wrk(:),fpu(:),fpv(:)
        integer, allocatable :: nru(:),nrv(:),nrdatu(:),nrdatv(:)
        integer :: mu,mv,mr,nuest,nvest,nc,nu,nv,numax,lwrk,ifsu,ifsv,ifbu,ifbv
        integer :: iopt(3),ider(4),idd(4),i,ju,nn,ier,lastdi,nplusu,nplusv,lastu0,lastu1
        real(dp) :: dr(6),step(2),r0,r1,p,range,fp0,fpold,reducu,reducv
        allocate(tu(0),tv(0),c(0))
        status = 3
        fp = zero
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
        if (s > zero) then
            allocate(nrdatu(nuest),nrdatv(nvest))
            nrdatu = 0
            nrdatv = 0
            nu = 0
            nv = 0
            lastdi = 0
            nplusu = 0
            nplusv = 0
            lastu0 = 0
            lastu1 = 0
            fp0 = zero
            fpold = zero
            reducu = zero
            reducv = zero
            call fpspgr(iopt,ider,u,mu,v,mv,data,mr,r0,r1,s,nuest,nvest, &
                0.001_dp,20,nc,nu,xknots,nv,yknots,coef,fp,fp0,fpold,reducu,reducv, &
                fpu,fpv,dr,step,lastdi,nplusu,nplusv,lastu0,lastu1,nru,nrv,nrdatu, &
                nrdatv,wrk,lwrk,ier)
            if (ier > 0) return
        else
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
            fp = zero
        end if
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
      pure subroutine fpspgr(iopt,ider,u,mu,v,mv,r,mr,r0,r1,s, &
                             nuest,nvest,tol,maxit,nc,nu,tu,nv,tv,c,fp,fp0,fpold,reducu, &
                             reducv,fpintu,fpintv,dr,step,lastdi,nplusu,nplusv,lastu0, &
                             lastu1,nru,nrv,nrdatu,nrdatv,wrk,lwrk,ier)
        implicit none

      !  ..
      !  ..scalar arguments..
        integer, intent(in) :: mu !! Number of interior colatitude sites.
        integer, intent(in) :: mv !! Number of periodic longitude sites.
        integer, intent(in) :: mr !! Observation vector length.
        integer, intent(in) :: nuest !! Allocated colatitude knot capacity.
        integer, intent(in) :: nvest !! Allocated longitude knot capacity.
        integer, intent(in) :: maxit !! Maximum smoothing-parameter iterations.
        integer, intent(in) :: nc !! Allocated tensor coefficient capacity.
        integer, intent(in) :: lwrk !! Available numerical work entries.
        integer, intent(inout) :: nu !! Active colatitude knot count.
        integer, intent(inout) :: nv !! Active longitude knot count.
        integer, intent(inout) :: &
            lastdi !! Last refined axis: zero none, one latitude, two longitude.
        integer, intent(inout) :: nplusu !! Last number of added latitude knots.
        integer, intent(inout) :: nplusv !! Last number of added longitude knots.
        integer, intent(inout) :: lastu0 !! Cached north-pole range window size.
        integer, intent(inout) :: lastu1 !! Cached south-pole range window size.
        integer, intent(out) :: ier !! FITPACK termination code; nonpositive denotes success.
        real(dp), intent(in) :: r0 !! Supplied or estimated scaled north-pole observation.
        real(dp), intent(in) :: r1 !! Supplied or estimated scaled south-pole observation.
        real(dp), intent(in) :: s !! Nonnegative scaled squared-residual target.
        real(dp), intent(in) :: tol !! Relative residual target tolerance.
        real(dp), intent(inout) :: fp !! Squared grid residual of the current fit.
        real(dp), intent(inout) :: fp0 !! Squared residual of the polynomial fit.
        real(dp), intent(inout) :: fpold !! Previous knot-set squared residual.
        real(dp), intent(inout) :: reducu !! Last residual reduction from latitude refinement.
        real(dp), intent(inout) :: reducv !! Last residual reduction from longitude refinement.
      !  ..array arguments..
        integer, intent(in) :: iopt(3) !! Fit mode and north/south continuity flags.
        integer, intent(in) :: ider(4) !! Pole observation exactness and zero-slope flags.
        integer, intent(inout) :: nru(mu) !! Cached latitude interval for each site.
        integer, intent(inout) :: nrv(mv) !! Cached longitude interval for each site.
        integer, intent(inout) :: nrdatu(nuest) !! Eligible latitude data counts by knot interval.
        integer, intent(inout) :: nrdatv(nvest) !! Eligible longitude data counts by knot interval.
        real(dp), intent(in) :: u(mu) !! Interior colatitude sites in radians.
        real(dp), intent(in) :: v(mv) !! Periodic longitude sites in radians.
        real(dp), intent(in) :: r(mr) !! Longitude-fastest scaled grid observations.
        real(dp), intent(inout) :: &
            fpintu(nuest) !! Latitude-interval squared residual contributions.
        real(dp), intent(inout) :: &
            fpintv(nvest) !! Longitude-interval squared residual contributions.
        real(dp), intent(inout) :: dr(6) !! Pole values and two derivative components per pole.
        real(dp), intent(inout) :: tu(nuest) !! Full latitude knot vector and unused capacity.
        real(dp), intent(inout) :: tv(nvest) !! Full longitude knot vector and unused capacity.
        real(dp), intent(inout) :: c(nc) !! Longitude-fastest tensor coefficients.
        real(dp), intent(inout) :: wrk(lwrk) !! Caller-owned numerical workspace.
        real(dp), intent(inout) :: step(2) !! Finite-difference optimization increments per pole.
      !  ..local scalars..
      real(dp) :: acc,fpms,f1,f2,f3,p,p1,p2,p3,vb,ve,rmax,rmin,rn
      integer :: i,ifbu,ifbv,ifsu,ifsv,istart,iter,i1,i2,j,ju,ktu,l,l1,l2,l3,l4,mpm,mumin, &
                 mu0,mu1,nn,nplu,nplv,npl1,nrintu,nrintv,nue,numax,nve,nvmax
      !  ..local arrays..
      logical :: check1,check3,success
      integer :: idd(4)
      real(dp) :: drr(6)

      !   set constants
      real(dp), parameter :: period = pi2

      !   initialization
      ier   = fitpack_ok
      ifsu  = 0
      ifsv  = 0
      ifbu  = 0
      ifbv  = 0
      nplu  = 0
      nplv  = 0
      p     = -one
      mumin = 4 - merge(1,0,ider(1)>=0) - merge(1,0,iopt(2)==1 .and. ider(2)==1) &
                - merge(1,0,ider(3)>=0) - merge(1,0,iopt(3)==1 .and. ider(4)==1)
      if (mumin==0) mumin = 1

      vb = v(1)
      ve = vb+period

      ! acc denotes the absolute tolerance for the root of f(p)=s.
      acc = tol*s

      ! numax and nvmax denote the number of knots needed for interpolation.
      numax = mu+6+iopt(2)+iopt(3)
      nvmax = mv+7
      nue   = min(numax,nuest)
      nve   = min(nvmax,nvest)

      ! initialize wrk(1:4)
      if (iopt(1)==0) then
        fp0       = zero
        fpold     = zero
        reducu    = zero
        reducv    = zero
      endif

      ! *****
      ! part 1: determination of the number of knots and their position.
      ! *****
! given a set of knots we compute the least-squares spline sinf(u,v) and the corresponding sum
! of
      !  squared residuals fp = f(p=inf).
      !  if iopt(1)=-1  sinf(u,v) is the requested approximation.
      !  if iopt(1)>=0  we check whether we can accept the knots:
      !    if fp <= s we will continue with the current set of knots.
! if fp >  s we will increase the number of knots and compute the corresponding least-squares
! spline
      !               until finally fp <= s.
      !    the initial choice of knots depends on the value of s and iopt.
! if s=0 we have spline interpolation; in that case the number of knots in the u-direction
! equals
      !     nu=numax=mu+6+iopt(2)+iopt(3) and in the v-direction nv=nvmax=mv+7.
      !    if s>0 and
! iopt(1)=0 we first compute the least-squares polynomial,i.e. a spline without interior knots:
      !                nu=8 ; nv=8.
! iopt(1)=1 we start with the set of knots found at the last call of the routine, except for the
      !                case that s > fp0; then we compute the least-squares polynomial directly.
      ! *****

      if (iopt(1)>=0) then

          if (s<=zero) then

              !  if s = 0, s(u,v) is an interpolating spline.
              nu = numax
              nv = nvmax

              !  test whether the required storage space exceeds the available one.
              if (nu>nuest .or. nv>nvest) then
                 ier = fitpack_insufficient_storage
                 return
              end if

              !  find the position of the knots in the v-direction.
              tv(1:mv+7) = [v(mv-2:mv)-period,v(1:mv),ve,v(2:4)+period]

              !  if not all the derivative values g(i,j) are given, we will first
              !  estimate these values by computing a least-squares spline
              idd = ider
              if(idd(1)==0) idd(1) = 1
              if(idd(1)>0)  dr(1) = r0
              if(idd(3)==0) idd(3) = 1
              if(idd(3)>0)  dr(4) = r1
              if(ider(1)<0 .or. (iopt(2)/=0 .and. ider(2)==0) .or. &
                 ider(3)<0 .or. (iopt(3)/=0 .and. ider(4)==0)) then

                  ! we set up the knots in the u-direction for computing the least-squares spline.
                  i1 = 3
                  i2 = mu-2
                  nu = 4
                  do i=1,mu
                     if (i1>i2) exit
                     nu = nu+1
                     tu(nu) = u(i1)
                     i1 = i1+2
                  end do
                  do i=1,4
                     tu(i) = zero
                     nu = nu+1
                     tu(nu) = pi
                  end do

                  ! we compute the least-squares spline for estimating the derivatives.
                  call fpopsp(ifsu,ifsv,ifbu,ifbv,u,mu,v,mv,r,mr,r0,r1,dr,iopt,idd,         &
                              tu,nu,tv,nv,nuest,nvest,p,step,c,nc,fp,fpintu,fpintv,nru,nrv, &
                              wrk,lwrk)
                  ifsu = 0

              end if

              ! if all the derivatives at the origin are known, we compute the
              ! interpolating spline.
              ! we set up the knots in the u-direction, needed for interpolation.
              nn = numax-8
              if (nn/=0) then
                  ju = 2-iopt(2)
                  do l=1,nn
                    tu(l+4) = u(ju)
                    ju = ju+1
                  end do
                  nu = numax
                  l  = nu
                  do i=1,4
                     tu(i) = zero
                     tu(l) = pi
                     l = l-1
                  end do
              endif

              ! we compute the interpolating spline.
              call fpopsp(ifsu,ifsv,ifbu,ifbv,u,mu,v,mv,r,mr,r0,r1,dr,iopt,idd,         &
                          tu,nu,tv,nv,nuest,nvest,p,step,c,nc,fp,fpintu,fpintv,nru,nrv, &
                          wrk,lwrk)

              ier = fitpack_interpolating_ok
              fp  = zero
              return

          else

              !  if s>0 our initial choice of knots depends on the value of iopt(1).
              ier = fitpack_ok

              if (iopt(1)/=0) step = -step

              if (fp0<=s .or. iopt(1)==0) then

                  !  if iopt(1)=0 or iopt(1)=1 and fp0 <= s,we start computing the least-
                  !  squares polynomial (which is a spline without interior knots).
                  ier       = fitpack_leastsquares_ok
                  idd       = [ider(1),ione,ider(3),ione]
                  nu        = 8
                  nv        = 8
                  nrdatu(1) = mu-2+iopt(2)+iopt(3)
                  nrdatv(1) = mv-1
                  lastdi    = knot_dim_none
                  nplusu    = 0
                  nplusv    = 0
                  fp0       = zero
                  fpold     = zero
                  reducu    = zero
                  reducv    = zero
                  lastu0    = 0     ! force step(1)/step(2) (rmax-rmin per pole) to be computed on
                  lastu1    = 0     ! the first pass; step and the last* caches are inout persistent
! state, uninitialized on a fresh fit (f77 relied on zeroed statics)
                  step      = zero

              else

                  !  if iopt(1)=1 and fp0 > s we start computing the least-squares spline
                  !  according to the set of knots found at the last call of the routine.
                  !  we determine the number of grid coordinates u(i) inside each knot
                  !  interval (tu(l),tu(l+1)).
                  l = 5
                  j = 1
                  nrdatu(1) = 0
                  mu0 = 2-iopt(2)
                  mu1 = mu-1+iopt(3)
                  do i=mu0,mu1
                    nrdatu(j) = nrdatu(j)+1
                    if (u(i)>=tu(l)) then
                        nrdatu(j) = nrdatu(j)-1
                        l = l+1
                        j = j+1
                        nrdatu(j) = 0
                    endif
                  end do
                  !  we determine the number of grid coordinates v(i) inside each knot
                  !  interval (tv(l),tv(l+1)).
                  l = 5
                  j = 1
                  nrdatv(1) = 0
                  do i=2,mv
                    nrdatv(j) = nrdatv(j)+1
                    if (v(i)>=tv(l)) then
                       nrdatv(j) = nrdatv(j)-1
                       l = l+1
                       j = j+1
                       nrdatv(j) = 0
                    endif
                  end do
                  idd = ider

              endif

          endif

      endif

      !  main loop for the different sets of knots.mpm=mu+mv is a safe upper
      !  bound for the number of iterations
      mpm = mu+mv
      iterations: do iter=1,mpm
         ! find nrintu (nrintv) which is the number of knot intervals in the
         ! u-direction (v-direction).
         nrintu = nu-7
         nrintv = nv-7

         ! find the position of the additional knots which are needed for the
         ! b-spline representation of s(u,v).
         i = nu
         do j=1,4
            tu(j) = zero
            tu(i) = pi
            i = i-1
         end do
         l1 = 4
         l2 = l1
         l3 = nv-3
         l4 = l3
         tv(l2) = vb
         tv(l3) = ve
         do j=1,3
            l1 = l1+1
            l2 = l2-1
            l3 = l3+1
            l4 = l4-1
            tv(l2) = tv(l4)-pi2
            tv(l3) = tv(l1)+pi2
         end do

! find an estimate of the range of possible values for the optimal derivatives at the origin.
         ktu = max(mumin,nrdatu(1)+2-iopt(2))
         if (ktu/=lastu0) then
             l       = min(mr,mv*ktu)
             rmin    = min(r0,minval(r(1:l)))
             rmax    = max(r0,maxval(r(1:l)))
             step(1) = rmax-rmin
             lastu0  = ktu
         endif

         ktu = max(mumin,nrdatu(nrintu)+2-iopt(3))
         if (ktu/=lastu1) then
            l    = min(mr,mv*ktu)
            rmin = min(r1,minval(r(mr-l+1:mr)))
            rmax = max(r1,maxval(r(mr-l+1:mr)))
            step(2) = rmax-rmin
            lastu1 = ktu
         endif

         ! find the least-squares spline sinf(u,v).
         call fpopsp(ifsu,ifsv,ifbu,ifbv,u,mu,v,mv,r,mr,r0,r1,dr,iopt, &
                     idd,tu,nu,tv,nv,nuest,nvest,p,step,c,nc,fp,fpintu,fpintv,nru, &
                     nrv,wrk,lwrk)

         step = abs(step)
         if (ier==fitpack_leastsquares_ok) fp0 = fp

         ! test whether the least-squares spline is an acceptable solution.
         if (iopt(1)<0) return
         fpms = fp-s; if (abs(fpms) < acc) return ! success

         ! if f(p=inf) < s, we accept the choice of knots.
         if (fpms<zero) exit iterations

         ! if nu=numax and nv=nvmax, sinf(u,v) is an interpolating spline
         if(nu==numax .and. nv==nvmax) then
            ier = fitpack_interpolating_ok
            fp  = zero
            return
         endif

         ! increase the number of knots.
         ! if nu=nue and nv=nve we cannot further increase the number of knots
         ! because of the storage capacity limitation.
         if (nu==nue .and. nv==nve) then
            ier = fitpack_insufficient_storage
            return
         endif

         if(ider(1)==0) fpintu(1)      = fpintu(1)     +(r0-dr(1))**2
         if(ider(3)==0) fpintu(nrintu) = fpintu(nrintu)+(r1-dr(4))**2
         ier = fitpack_ok

         ! adjust the parameter reducu or reducv according to the direction
         ! in which the last added knots were located.
         first_knot: if (lastdi==knot_dim_none) then

            nplv = 3
            idd(2) = ider(2)
            idd(4) = ider(4)
            fpold = fp

         else first_knot

            select case (lastdi)
               case (knot_dim_1); reducu = fpold-fp
               case (knot_dim_2); reducv = fpold-fp
            end select

            ! store the sum of squared residuals for the current set of knots.
            fpold = fp

            ! find nplu, the number of knots we should add in the u-direction.
            nplu = 1
            if (nu/=8) then
               rn   = nplusu
               npl1 = nplusu*2
! guard the division: skip it when reducu<=acc
               if (reducu>acc) npl1 = fp_knots_to_add(rn,fpms,reducu)
               nplu = min(nplusu*2,max(npl1,nplusu/2,1))
            endif

            ! find nplv, the number of knots we should add in the v-direction.
            nplv = 3
            if (nv/=8) then
               rn   = nplusv
               npl1 = nplusv*2
! guard the division: skip it when reducv<=acc
               if (reducv>acc) npl1 = fp_knots_to_add(rn,fpms,reducv)
               nplv = min(nplusv*2,max(npl1,nplusv/2,1))
            endif

         endif first_knot

         ! test whether we are going to add knots in the u- or v-direction.
         lastdi = new_knot_dimension(nu,nplu,nue,nv,nplv,nve,lastdi)

         choose_dim: if (lastdi==knot_dim_2) then

            ! addition in the v-direction.
            nplusv = nplv
            ifsv   = 0

            add_v_knots: do l=1,nplusv

               ! add a new knot in the v-direction.
               istart = 1
               call fpknot(v,mv,tv,nv,fpintv,nrdatv,nrintv,nvest,istart)

               ! test whether we cannot further increase the number of knots in the v-direction.
               if (nv==nve) exit add_v_knots

            end do add_v_knots

        else choose_dim

            ! addition in the u-direction.
            nplusu = nplu
            ifsu   = 0
            istart = merge(1,0,iopt(2)==0)
            add_u_knots: do l=1,nplusu

               ! add a new knot in the u-direction
               call fpknot(u,mu,tu,nu,fpintu,nrdatu,nrintu,nuest,istart)

               ! test whether we cannot further increase the number of knots in the u-direction.
               if (nu==nue) exit add_u_knots

            end do add_u_knots

        endif choose_dim

      ! restart the computations with the new set of knots.
      end do iterations

      ! test whether the least-squares polynomial is a solution of our approximation problem.
      if (ier==fitpack_leastsquares_ok) return

! **********************************************************************************************
      ! part 2: determination of the smoothing spline sp(u,v)
! **********************************************************************************************
      !  we have determined the number of knots and their position. we now compute the b-spline
! coefficients of the smoothing spline sp(u,v). this smoothing spline depends on the parameter
      !  p in such a way that
      !    f(p) = sumi=1,mu(sumj=1,mv((z(i,j)-sp(u(i),v(j)))**2)
      !  is a continuous, strictly decreasing function of p. moreover the least-squares polynomial
      !  corresponds to p=0 and the least-squares spline to p=infinity. then iteratively we have to
! determine the positive value of p such that f(p)=s. the process which is proposed here makes
! use of rational interpolation. f(p) is approximated by a rational function r(p)=(u*p+v)/(p+w);
      !  three values of p (p1,p2,p3) with corresponding values of f(p) (f1=f(p1)-s,f2=f(p2)-s,
      !  f3=f(p3)-s) are used to calculate the new value of p such that r(p)=s.
      !  convergence is guaranteed by taking f1 > 0 and f3 < 0.
! **********************************************************************************************

      !  initial value for p.
      p1   = zero
      f1   = fp0-s
      p3   = -one
      f3   = fpms
      p    = one
      drr  = dr
      check1 = fp_false
      check3 = fp_false

      !  iteration process to find the root of f(p)=s.
      find_root: do iter = 1,maxit

         ! find the smoothing spline sp(u,v) and the corresponding sum f(p).
         call fpopsp(ifsu,ifsv,ifbu,ifbv,u,mu,v,mv,r,mr,r0,r1,drr,iopt, &
                     idd,tu,nu,tv,nv,nuest,nvest,p,step,c,nc,fp,fpintu, &
                     fpintv,nru,nrv,wrk,lwrk)

         ! test whether the approximation sp(u,v) is an acceptable solution.
         fpms = fp-s; if (abs(fpms)<=acc) return ! success

         ! find the new value of p and carry out one more step.
         call root_finding_iterate(p1,f1,p2,f2,p3,f3,p,fpms,acc,check1,check3,success)
         if (.not.success) then
            ier = fitpack_s_too_small
            return
         end if

      end do find_root

      ! max number of iterations exceeded
      ier = fitpack_maxit
      return

      end subroutine fpspgr

      pure subroutine fpknot(x,m,t,n,fpint,nrdata,nrint,nest,istart)
        implicit none

      !  ..scalar arguments..
        integer, intent(in) :: m !! Number of coordinate sites.
        integer, intent(in) :: nest !! Allocated knot and interval capacity.
        integer, intent(in) :: istart !! Offset of the first eligible coordinate.
        integer, intent(inout) :: n !! Active knot count.
        integer, intent(inout) :: nrint !! Number of active knot intervals.
      !  ..array arguments..
        real(dp), intent(in) :: x(m) !! Ordered coordinates for eligible knot insertion.
        real(dp), intent(inout) :: t(nest) !! Knot vector updated by insertion.
        real(dp), intent(inout) :: fpint(nest) !! Interval residual estimates updated by insertion.
        integer,  intent(inout) :: &
            nrdata(nest) !! Eligible site counts per interval updated by insertion.

      !  ..local scalars..
      real(dp) :: an,am,fpmax
      integer :: ihalf,j,jbegin,jj,jk,jpoint,k,maxbeg,maxpt,next,nrx,number
      !  ..
      number = 0
      maxpt  = 0
      maxbeg = 0
      k      = (n-nrint-1)/2
! search for knot interval t(number+k) <= x <= t(number+k+1) where fpint(number) is maximal on
! the
! condition that nrdata(number)/=0. an interval is eligible only if it holds at least one
! interior
! data point (jpoint/=0); the first eligible interval seeds the search (number==0), after which
! any
! interval with a strictly larger residual takes over. seeding on the first eligible interval,
! rather
! than only on fpmax<fpint(j), guarantees a valid selection even when every residual is exactly
! zero
! (an all-zero fpint would otherwise leave number==0 and fall through to use of uninitialized
! am).
      fpmax  = zero
      jbegin = istart
      do j=1,nrint
        jpoint = nrdata(j)

        if (jpoint/=0 .and. (number==0 .or. fpint(j)>fpmax)) then
           fpmax  = fpint(j)
           number = j
           maxpt  = jpoint
           maxbeg = jbegin
        endif

        jbegin = jbegin+jpoint+1
      end do

! no interval carries an interior data point: there is nothing to refine, so leave the knot set
! as is
      if (number==0) return

      !  let coincide the new knot t(number+k+1) with a data point x(nrx)
      !  inside the old knot interval t(number+k) <= x <= t(number+k+1).
      ihalf = maxpt/2+1
      nrx   = maxbeg+ihalf
      next  = number+1

      !  adjust the different parameters.
      if (next<=nrint) then
         do j=next,nrint
            jj = next+nrint-j
            fpint(jj+1) = fpint(jj)
            nrdata(jj+1) = nrdata(jj)
            jk = jj+k
            t(jk+1) = t(jk)
         end do
      endif

      !  split the selected interval in two, apportioning its residual between the halves.
      !  number>0 is guaranteed by the early return above, so am is always initialized here.
      nrdata(number) = ihalf-1
      nrdata(next)   = maxpt-ihalf
      am = maxpt
      an = nrdata(number)
      fpint(number) = fpmax*an/am

      an = nrdata(next)
      fpint(next) = fpmax*an/am
      jk = next+k
      t(jk) = x(nrx)

      n     = n+1
      nrint = nrint+1

      end subroutine fpknot

      elemental integer function fp_knots_to_add(rn,fpms,reduc) result(npl1)
        implicit none
        real(dp), intent(in) :: rn !! Previous knot addition count, represented as real.
        real(dp), intent(in) :: fpms !! Current squared residual minus target.
        real(dp), intent(in) :: reduc !! Positive last residual reduction.

         real(dp), parameter :: npl_max = real(huge(0),dp)
         real(dp) :: ratio

         ratio = rn*fpms/reduc
         if (ratio<npl_max) then
            npl1 = int(max(ratio,zero))
         else
            npl1 = huge(0) ! saturate; a nan ratio lands here too
         end if

      end function fp_knots_to_add

      elemental integer function new_knot_dimension(n1,n1add,n1max,n2,n2add,n2max,last) result(dir)
        implicit none
        integer, intent(in) :: n1 !! Current first-axis knot count.
        integer, intent(in) :: n2 !! Current second-axis knot count.
        integer, intent(in) :: n1add !! Proposed first-axis knot additions.
        integer, intent(in) :: n2add !! Proposed second-axis knot additions.
        integer, intent(in) :: n1max !! First-axis knot limit.
        integer, intent(in) :: n2max !! Second-axis knot limit.
        integer, intent(in) :: last !! Previously refined axis, zero if none.

         ! Alternate equal-cost axis choices, following the spherical FITPACK path.
         if (n2<n2max .and. (n2add<n1add .or. n1>=n1max .or. last==knot_dim_none .or. &
             (n2add==n1add .and. last==knot_dim_1))) then
            dir = knot_dim_2
         else
            dir = knot_dim_1
         end if

      end function new_knot_dimension

      elemental subroutine root_finding_iterate(p1,f1,p2,f2,p3,f3,p,fpms,acc,check1,check3,success)
        implicit none
        real(dp), intent(inout) :: p1 !! Lower smoothing-parameter bracket.
        real(dp), intent(inout) :: f1 !! Residual minus target at the lower bracket.
        real(dp), intent(inout) :: p2 !! Current smoothing-parameter sample.
        real(dp), intent(inout) :: f2 !! Residual minus target at the current sample.
        real(dp), intent(inout) :: p3 !! Upper bracket; negative denotes infinity.
        real(dp), intent(inout) :: f3 !! Residual minus target at the upper bracket.
        real(dp), intent(inout) :: p !! Positive smoothing-parameter iterate.
        real(dp), intent(in) :: fpms !! Current squared residual minus target.
        real(dp), intent(in) :: acc !! Absolute residual-target tolerance.
        logical, intent(inout) :: check1 !! Whether the lower bracket has been established.
        logical, intent(inout) :: check3 !! Whether the upper bracket has been established.
        logical, intent(out) :: &
            success !! False when residual samples violate the required ordering.

          !  set constants
          real(dp), parameter :: con1 = 0.1e0_dp
          real(dp), parameter :: con9 = 0.9e0_dp
          real(dp), parameter :: con4 = 0.4e-01_dp

          success = .true.

          p2 = p
          f2 = fpms
          if (.not.check3) then
             if ((f2-f3)>acc) then
                check3=f2<zero
             else
                ! our initial choice of p is too large.
                p3 = p2
                f3 = f2
                p  = p*con4
                if (p<=p1) p=p1*con9 + p2*con1
                return
             endif
          endif

          if (.not.check1) then
             if ((f1-f2)>acc) then
                 check1 = f2>zero
             else
                 ! our initial choice of p is too small
                 p1 = p2
                 f1 = f2
                 p  = p/con4
                 if (p3>=zero .and. p>=p3) p = p2*con1 + p3*con9
                 return
             endif
          endif

          ! test whether the iteration process proceeds as theoretically expected.
          if (f2>=f1 .or. f2<=f3) then
             success = .false.
             return
          else
             ! find the new value of p.
             call fprati(p1,f1,p2,f2,p3,f3,p)
          endif

      end subroutine root_finding_iterate

      elemental subroutine fprati(p1,f1,p2,f2,p3,f3,p)
        implicit none
      !  ..scalar arguments..
        real(dp), intent(inout) :: p1 !! Lower smoothing-parameter bracket.
        real(dp), intent(inout) :: f1 !! Residual minus target at the lower bracket.
        real(dp), intent(inout) :: p3 !! Upper bracket; negative denotes infinity.
        real(dp), intent(inout) :: f3 !! Residual minus target at the upper bracket.
        real(dp), intent(in) :: p2 !! Current smoothing-parameter sample.
        real(dp), intent(in) :: f2 !! Residual minus target at the current sample.
        real(dp), intent(out) :: p !! Positive smoothing-parameter iterate.

      !  ..local scalars..
      real(dp) :: h1,h2,h3
      !  ..
      if (p3>zero) then
         h1 = f1*(f2-f3)
         h2 = f2*(f3-f1)
         h3 = f3*(f1-f2)
         p = -(p1*p2*h3+p2*p3*h1+p3*p1*h2)/(p1*h1+p2*h2+p3*h3)
      else
         !  value of p in case p3 = infinity.
         p = (p1*(f1-f3)*f2-p2*(f2-f3)*f1)/((f1-f2)*f3)
      end if

      !  adjust the value of p1,f1,p3 and f3 such that f1 > 0 and f3 < 0.
      if (f2>=zero) then
        p1 = p2
        f1 = f2
      else
        p3 = p2
        f3 = f2
      endif

      return
      end subroutine fprati
end module scifort_interpolate_sphere_grid_core
