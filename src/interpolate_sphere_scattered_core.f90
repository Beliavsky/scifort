! SPDX-License-Identifier: BSD-3-Clause
! Copyright (c) 2022, Federico Perini. All rights reserved.
! Copyright (c) 2026 SciFort contributors
!
! Adapted from perazz/fitpack src/fitpack_core.F90, commit
! 1907c7bad9f5c4557970e906bdd4909eb52b2d4c, based on Paul Dierckx's FITPACK.
! Full terms: LICENSES/perazz-fitpack-BSD-3-Clause.txt.
module scifort_interpolate_sphere_scattered_core
    use scifort_kinds, only : dp
    use scifort_constants, only : pi => scifort_pi
    use scifort_interpolate_sphere_grid_core, only : equal, not_equal, fp_knot_interval, &
        fpback, fpbspl, fpdisc, fpgivs, fprota, root_finding_iterate
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    implicit none
    private
    public :: fitpack_sphere_scattered
    integer, parameter :: max_order = 19, degree_3 = 3, degree_5 = 5
    integer, parameter :: fitpack_ok = 0, fitpack_interpolating_ok = -1
    integer, parameter :: fitpack_leastsquares_ok = -2, fitpack_insufficient_storage = 1
    integer, parameter :: fitpack_s_too_small = 2, fitpack_maxit = 3
    integer, parameter :: fitpack_too_many_knots = 4, fitpack_overlapping_knots = 5
    real(dp), parameter :: zero = 0.0_dp, one = 1.0_dp, two = 2.0_dp, three = 3.0_dp
    real(dp), parameter :: half = 0.5_dp, onep5 = 1.5_dp, pi2 = 2.0_dp*pi
    real(dp), parameter :: ten = 10.0_dp
    logical, parameter :: fp_false = .false.
contains
    subroutine fitpack_sphere_scattered(theta,phi,r,w,s,eps,ntest,npest,lwork1,lwork2, &
        tt,tp,c,fp,rank,ier)
        implicit none
        real(dp), intent(in) :: theta(:) !! Scattered colatitudes in [0,pi], radians.
        real(dp), intent(in) :: phi(:) !! Paired longitudes in [0,2*pi], radians.
        real(dp), intent(in) :: r(:) !! Finite scaled observations, at least two.
        real(dp), intent(in) :: w(:) !! Positive scaled residual weights.
        real(dp), intent(in) :: s !! Nonnegative scaled squared-residual target.
        real(dp), intent(in) :: eps !! Squared relative diagonal threshold in (0,1).
        integer, intent(in) :: ntest !! Maximum full colatitude knot count, at least eight.
        integer, intent(in) :: npest !! Maximum full longitude knot count, at least eight.
        integer, intent(in) :: lwork1 !! Primary real workspace size, checked by caller.
        integer, intent(in) :: lwork2 !! Rank-reduction real workspace size, checked by caller.
        real(dp), allocatable, intent(out) :: tt(:) !! Full active colatitude knot vector.
        real(dp), allocatable, intent(out) :: tp(:) !! Full active longitude knot vector.
        real(dp), allocatable, intent(out) :: c(:) !! Longitude-fastest physical tensor layout.
        real(dp), intent(out) :: fp !! Weighted squared residual in scaled units.
        integer, intent(out) :: rank !! Estimated final factor rank; -1 for polynomial fit.
        integer, intent(out) :: ier !! Raw FITPACK code, nonpositive means a valid solution.
        real(dp), allocatable :: tx(:),ty(:),coef(:),work1(:),work2(:)
        integer, allocatable :: iwork(:)
        integer :: m,ntt,npp,ncc,nreg,nrint,nc,ib1,ib3,nt,np,la,lf,lff,lfp,lco,lh
        integer :: lq,lbt,lbp,lro,lcc,lcs,lst,lsp
        real(dp) :: sup
        allocate(tt(0),tp(0),c(0))
        m = size(theta)
        ntt = ntest-7
        npp = npest-7
        ncc = 6+npp*(ntt-1)
        nreg = ntt*npp
        nrint = ntt+npp
        nc = (ntest-4)*(npest-4)
        ib1 = max(4*npp,6+3*npp)
        ib3 = max(ib1+3,6+3*npp)
        allocate(tx(ntest),ty(npest),coef(nc),work1(lwork1),work2(lwork2),iwork(m+nreg))
        tx = zero
        ty = zero
        coef = zero
        work1 = zero
        work2 = zero
        iwork = 0
        nt = 0
        np = 0
        fp = zero
        sup = zero
        lq = 2
        la = lq+ncc*ib3
        lf = la+ncc*ib1
        lff = lf+ncc
        lfp = lff+nc
        lco = lfp+nrint
        lh = lco+nrint
        lbt = lh+ib3
        lbp = lbt+5*ntest
        lro = lbp+5*npest
        lcc = lro+npest
        lcs = lcc+npest
        lst = lcs+npest
        lsp = lst+4*m
        call fpsphe(0,m,theta,phi,r,w,s,ntest,npest,eps,0.001_dp,20,ib1,ib3,nc,ncc, &
            nrint,nreg,nt,tx,np,ty,coef,fp,sup,work1(lfp),work1(lco),work1(lf), &
            work1(lff),work1(lro),work1(lcc),work1(lcs),work1(la),work1(lq), &
            work1(lbt),work1(lbp),work1(lst),work1(lsp),work1(lh),iwork(m+1), &
            iwork(1),work2,lwork2,ier,rank)
        if (ier > 0) return
        nc = (nt-4)*(np-4)
        if (.not. all(ieee_is_finite(coef(:nc))) .or. .not. ieee_is_finite(fp)) then
            ier = fitpack_maxit
            return
        end if
        tt = tx(:nt)
        tp = ty(:np)
        c = coef(:nc)
    end subroutine fitpack_sphere_scattered

    pure function weighted_sphere_residual(t,p,r,w,tt,tp,c,spt,spp) result(fp)
        implicit none
        real(dp), intent(in) :: t(:) !! Observation colatitudes in [0,pi], radians.
        real(dp), intent(in) :: p(:) !! Paired longitudes in [0,2*pi], radians.
        real(dp), intent(in) :: r(:) !! Finite scaled observations.
        real(dp), intent(in) :: w(:) !! Positive scaled residual weights.
        real(dp), intent(in) :: tt(:) !! Full active colatitude knot vector.
        real(dp), intent(in) :: tp(:) !! Full active longitude knot vector.
        real(dp), intent(in) :: c(:) !! Physical tensor coefficients, longitude varying fastest.
        real(dp), intent(in) :: spt(:,:) !! Four active latitude basis values per observation.
        real(dp), intent(in) :: spp(:,:) !! Four active longitude basis values per observation.
        real(dp) :: fp,value
        integer :: i,j,l,k,base,np4
        fp = zero
        np4 = size(tp)-4
        do i = 1, size(t)
            l = fp_knot_interval(tt,t(i),4,size(tt)-4)
            k = fp_knot_interval(tp,p(i),4,np4)
            base = (l-4)*np4+k-4
            value = zero
            do j = 1, 4
                value = value+spt(i,j)*dot_product(spp(i,:),c(base+1:base+4))
                base = base+np4
            end do
            fp = fp+(w(i)*(r(i)-value))**2
        end do
    end function weighted_sphere_residual

      pure subroutine fpsphe(iopt,m,teta,phi,r,w,s,ntest,npest,eta,tol,maxit, &
                             ib1,ib3,nc,ncc,intest,nrest,nt,tt,np,tp,c,fp,sup,fpint,coord,f, &
                             ff,row,coco,cosi,a,q,bt,bp,spt,spp,h,index,nummer,wrk,lwrk, &
            ier,rank_out)
        implicit none
        integer, intent(out) :: rank_out !! Factor rank; -1 for the initial polynomial.

      !  ..scalar arguments..
        integer, intent(in) :: iopt !! Fit mode; this interface uses zero for a fresh smoothing fit.
        integer, intent(in) :: m !! Number of coordinate sites.
        integer, intent(in) :: ntest !! Allocated full colatitude knot capacity.
        integer, intent(in) :: npest !! Allocated full longitude knot capacity.
        integer, intent(in) :: maxit !! Maximum smoothing-parameter iterations.
        integer, intent(in) :: ib1 !! Allocated observation factor bandwidth.
        integer, intent(in) :: ib3 !! Allocated penalized factor bandwidth.
        integer, intent(in) :: nc !! Allocated tensor coefficient capacity.
        integer, intent(in) :: ncc !! Allocated reduced coefficient capacity.
        integer, intent(in) :: intest !! Capacity of combined latitude/longitude interval storage.
        integer, intent(in) :: nrest !! Capacity of the coordinate panel index.
        integer, intent(in) :: lwrk !! Available numerical work entries.
        integer, intent(out) :: ier !! FITPACK termination code; nonpositive denotes success.
        integer, intent(inout) :: nt !! Active colatitude knot count.
        integer, intent(inout) :: np !! Active longitude knot count.
        real(dp), intent(in) :: s !! Nonnegative scaled squared-residual target.
        real(dp), intent(in) :: eta !! Relative squared rank threshold in (0,1).
        real(dp), intent(in) :: tol !! Relative residual target tolerance.
        real(dp), intent(inout) :: fp !! Squared grid residual of the current fit.
        real(dp), intent(inout) :: sup !! Squared residual of the initial two-parameter polynomial.

      !  ..array arguments..
        real(dp), intent(in) :: teta(m) !! Finite scattered colatitudes in [0,pi], radians.
        real(dp), intent(in) :: phi(m) !! Finite paired longitudes in [0,2*pi], radians.
        real(dp), intent(in) :: r(m) !! Longitude-fastest scaled grid observations.
        real(dp), intent(in) :: w(m) !! Finite strictly positive residual weights.
        real(dp), intent(inout) :: tt(ntest) !! Full colatitude knot vector and unused capacity.
        real(dp), intent(inout) :: &
            tp(npest) !! Full periodic longitude knot vector and unused capacity.
        real(dp), intent(inout) :: c(nc) !! Longitude-fastest tensor coefficients.
        real(dp), intent(inout) :: f(ncc) !! Transformed observation right hand side.
        real(dp), intent(inout) :: row(npest) !! Periodic longitude collocation row workspace.
        real(dp), intent(inout) :: coco(npest) !! Periodic spline coefficients representing cosine.
        real(dp), intent(inout) :: cosi(npest) !! Periodic spline coefficients representing sine.
        real(dp), intent(inout) :: &
            fpint(intest) !! Interval residual estimates updated by insertion.
        real(dp), intent(inout) :: &
            coord(intest) !! Residual-weighted coordinate sums by knot interval.
        real(dp), intent(inout) :: ff(nc) !! Transformed right hand side or coefficient workspace.
        real(dp), intent(inout) :: bt(ntest,5) !! Colatitude derivative-jump constraint matrix.
        real(dp), intent(inout) :: bp(npest,5) !! Longitude derivative-jump constraint matrix.
        real(dp), intent(inout) :: a(ncc,ib1) !! Banded triangular observation factor.
        real(dp), intent(inout) :: q(ncc,ib3) !! Banded triangular penalized factor.
        real(dp), intent(inout) :: &
            spt(m,4) !! Four active colatitude basis values for each observation.
        real(dp), intent(inout) :: &
            spp(m,4) !! Four active longitude basis values for each observation.
        real(dp), intent(inout) :: h(ib3) !! Active factor row workspace.
        real(dp), intent(inout) :: wrk(lwrk) !! Caller-owned numerical workspace.
        integer, intent(inout) :: index(nrest) !! Head observation index for each coordinate panel.
        integer, intent(inout) :: nummer(m) !! Next observation in each coordinate panel.

      !  ..local scalars..
      real(dp) :: aa,acc,arg,cn,co,c1,dmax,d1,d2,eps,fac1,fac2,facc,facs,fn,fpmax,fpms,f1, &
            f2,f3,htj, &
                     p,pinv,piv,p1,p2,p3,ri,si,sigma,sq,store,wi,rn
      integer :: i,iband,iband1,iband3,iband4,ii,ij,il,in,irot,iter,i1,i2,j,jlt,jrot,j1, &
            j2,l,la,lf,lh,ll,lp,&
                 lt,lwest,l1,l2,l4,ncof,ncoff,npp,np4,nreg,nrint,nrr,nr1,ntt,nt4,nt6,num,num1,rank
      logical :: check1,check3,success
      !  ..local arrays..
      real(dp), dimension(max_order+1) :: hp,ht

      rank_out = -1
      eps  = sqrt(eta)
      fpms = huge(zero)

      ! calculation of acc, the absolute tolerance for the root of f(p)=s.
      acc = tol*s

      ! initializations
      ier    = fitpack_ok
      lwest  = 0
      ntt    = 0
      iband1 = 0

      bootstrap: if (iopt>=0) then

          !  if iopt=0 we begin by computing the weighted least-squares polynomial of the form
          !     s(teta,phi) = c1*f1(teta) + cn*fn(teta)
          !  where f1(teta) and fn(teta) are the cubic polynomials satisfying
          !     f1(0) = 1, f1(pi) = f1'(0) = f1'(pi) = 0 ; fn(teta) = 1-f1(teta).
          !  the corresponding weighted sum of squared residuals gives the upper
          !  bound sup for the smoothing factor s.
          if (iopt==0 .or. (iopt>0 .and. s>=sup)) then

              sup  = zero
              d1   = zero
              d2   = zero
              c1   = zero
              cn   = zero
              aa   = zero
              fac1 = pi*(one + half)
              fac2 = two/pi**3
              initial_poly: do i=1,m
                  wi  = w(i)
                  ri  = r(i)*wi
                  arg = teta(i)
                  fn  = fac2*(fac1-arg)*arg**2
                  f1  = (one-fn)*wi
                  fn  = fn*wi
                  if (not_equal(fn,zero)) then
                      call fpgivs(fn,d1,co,si)
                      call fprota(co,si,f1,aa)
                      call fprota(co,si,ri,cn)
                  endif
                  if (not_equal(f1,zero)) then
                      call fpgivs(f1,d2,co,si)
                      call fprota(co,si,ri,c1)
                  endif
                  sup = sup+ri*ri
              end do initial_poly

              if (not_equal(d2,zero)) c1 = c1/d2
              if (not_equal(d1,zero)) cn = (cn-aa*c1)/d1

              !  find the b-spline representation of this least-squares polynomial
              nt       = 8
              np       = 8
              c(1:8)   = c1
              c(9:16)  = cn
              tt(1:4)  = zero
              tt(5:8)  = pi
              tp(1:4)  = zero
              tp(5:8)  = pi2
              fp       = sup

              ! test whether the least-squares polynomial is an acceptable solution
              fpms = sup-s
              if (fpms<acc) then
                  ier = fitpack_leastsquares_ok
                  return
              end if

          endif

          !  test whether we cannot further increase the number of knots.
          sixty: if (iopt==0 .or. (iopt>0 .and. (s>=sup .or. np<11))) then

              if (npest<11 .or. ntest<9) then
                  ier = fitpack_insufficient_storage
                  return
              end if

              !  find the initial set of interior knots of the spherical spline in case iopt = 0.
              np      = 11
              tp(5:7) = pi*[half,one,onep5]
              nt      = 9
              tt(5)   = tp(5)
          endif sixty

      endif bootstrap


! **********************************************************************************************
! **************
      !  part 1 : computation of least-squares spherical splines.
! **********************************************************************************************
! **************
! if iopt1<0 we compute the least-squares spherical spline according to the given set of knots.
! if iopt1>=0 we compute least-squares spherical splines with increasing numbers of knots until
! the
      !  corresponding sum f(p=inf)<=s.
      !  the initial set of knots then depends on the value of iopt1:
! if iopt1=0 we start with one interior knot in the teta-direction (pi/2) and and three in the
      !               phi-direction (pi/2,pi,3*pi/2).
      !    if iopt1>0 we start with the set of knots found at the last call of the routine.
! **********************************************************************************************
! **************

! main loop for the different sets of knots. m is a safe upper bound for the number of trials.
      compute_knots: do iter=1,m

          ! find the position of the additional knots which are needed for the
          ! b-spline representation of s(teta,phi).
          tp(1:4)     = [tp(np-6:np-4)-pi2, zero]
          tp(np-3:np) = [pi2, tp(5:7)+pi2]

          tt(1:4)     = zero
          tt(nt-3:nt) = pi

          ! find nrint, the total number of knot intervals and nreg, the number
          ! of panels in which the approximation domain is subdivided by the
          ! intersection of knots.
          ntt   = nt-7
          npp   = np-7
          nrr   = npp/2
          nr1   = nrr+1
          nrint = ntt+npp
          nreg  = ntt*npp

          ! arrange the data points according to the panel they belong to.
          call fporde(teta,phi,m,degree_3,degree_3,tt,nt,tp,np,nummer,index,nreg)

          ! find the b-spline coefficients coco and cosi of the cubic spline
          ! approximations sc(phi) and ss(phi) for cos(phi) and sin(phi).
          coco(1:npp) = zero
          cosi(1:npp) = zero
          a(1:npp,1:npp) = zero

          ! the coefficients coco and cosi are obtained from the conditions
          ! sc(tp(i))=cos(tp(i)),resp. ss(tp(i))=sin(tp(i)),i=4,5,...np-4.
          get_coefs: do i=1,npp
              l2 = i+3
              arg = tp(l2)
              hp = fpbspl(tp,np,degree_3,arg,l2)

              row(1:npp) = zero
              ll = i
              do j=1,3
                 if (ll>npp) ll = 1
                 row(ll) = row(ll)+hp(j)
                 ll = ll+1
              end do

              facc = cos(arg)
              facs = sin(arg)

              do j=1,npp
                  piv = row(j)
                  if (equal(piv,zero)) cycle

                  call fpgivs(piv,a(j,1),co,si)
                  call fprota(co,si,facc,coco(j))
                  call fprota(co,si,facs,cosi(j))
                  if (j<npp) then
                      j1 = j+1
                      i2 = 1
                      do l=j1,npp
                         i2 = i2+1
                         call fprota(co,si,row(l),a(j,i2))
                      end do
                  endif
              end do
          end do get_coefs

          coco(:npp) = fpback(a,coco,npp,npp,ncc)
          cosi(:npp) = fpback(a,cosi,npp,npp,ncc)

          ! find ncof, the dimension of the spherical spline and ncoff, the
          ! number of coefficients in the standard b-spline representation.
          nt4 = nt-4
          np4 = np-4
          ncoff = nt4*np4
          ncof = 6+npp*(ntt-1)

          ! find the bandwidth of the observation matrix a.
          if (ntt<4) then
             iband = ncof
          elseif (ntt==4) then
             iband = 3*(npp+1)
          else
             iband = 4*npp
          end if
          iband1 = iband-1

          ! initialize the observation matrix a.
          f(1:ncof) = zero
          a(1:ncof,1:iband) = zero

          ! initialize the sum of squared residuals.
          fp = zero

          ! fetch the data points in the new order. main loop for the different panels.
          panels: do num=1,nreg
              ! fix certain constants for the current panel; jrot records the column
              ! number of the first non-zero element in a row of the observation
              ! matrix according to a data point of the panel.
              num1 = num-1
              lt   = num1/npp
              l1   = lt+4
              lp   = num1-lt*npp+1
              l2   = lp+3
              lt   = lt+1
              jrot = merge(3+(lt-3)*npp,0,lt>2)

              ! test whether there are still data points in the current panel.
              in = index(num)
              points_left: do while (in/=0)
                  ! fetch a new data point.
                  wi = w(in)
                  ri = r(in)*wi

                  ! evaluate for the teta-direction, the 4 non-zero b-splines at teta(in)
                  ht = fpbspl(tt,nt,degree_3,teta(in),l1)

                  ! evaluate for the phi-direction, the 4 non-zero b-splines at phi(in)
                  hp = fpbspl(tp,np,degree_3,phi(in),l2)

                  ! store the value of these b-splines in spt and spp resp.
                  spp(in,1:4) = hp(1:4)
                  spt(in,1:4) = ht(1:4)

                  ! initialize the new row of observation matrix.
                  h(1:iband) = zero

                  ! calculate the non-zero elements of the new row by making the cross
                  ! products of the non-zero b-splines in teta- and phi-direction and
                  ! by taking into account the conditions of the spherical splines.
                  row(1:npp) = zero

                  ! take into account the condition (3) of the spherical splines.
                  ll = lp
                  do i=1,4
                     if (ll>npp) ll=1
                     row(ll) = row(ll)+hp(i)
                     ll = ll+1
                  end do

                  ! take into account the other conditions of the spherical splines.
                  if (lt<=2 .or. lt>=(ntt-1)) then
                     facc = dot_product(row(:npp),coco(:npp))
                     facs = dot_product(row(:npp),cosi(:npp))
                  else
                     facc = zero
                     facs = zero
                  endif

                  ! fill in the non-zero elements of the new row.
                  j1 = 0
                  new_row: do j =1,4
                     jlt = j+lt
                     htj = ht(j)
                     if (jlt==3) then
                         h(1) = h(1)+htj
                         h(2) = facc*htj
                         h(3) = facs*htj
                         j1 = 3
                     elseif (jlt==nt4) then
                         h(j1+1:j1+3) = htj*[facc,facs,one]
                         j1 = j1+2
                     elseif (jlt>2 .and. jlt<=nt4) then
                         h(j1+1:j1+npp) = row(1:npp)*htj
                         j1 = j1+npp
                     else
                         j1 = j1+1
                         h(j1) = h(j1)+htj
                     endif
                  end do new_row
                  h(:iband) = h(:iband)*wi

                  ! rotate the row into triangle by givens transformations.
                  call fp_rotate_row(h, iband, a, ri, f, jrot)

                  ! add the contribution of the row to the sum of squares of residual
                  ! right hand sides.
                  fp = fp+ri**2

                  ! find the number of the next data point in the panel.
                  in = nummer(in)
              end do points_left
          end do panels

          !  find dmax, the maximum value for the diagonal elements in the reduced triangle.
          dmax = max(zero,maxval(a(:ncof,1)))

          ! check whether the observation matrix is rank deficient.
          sigma = eps*dmax

          if (any(a(1:ncof,1)<=sigma)) then

                ! in case of rank deficiency, find the minimum norm solution.
                lwest = ncof*iband+ncof+iband
                if (lwest>lwrk) then
                   ier = lwest
                   return
                end if

                lf = 1
                lh = lf+ncof
                la = lh+iband
                ff(1:ncof) = f(1:ncof)
                q(1:ncof,1:iband) = a(1:ncof,1:iband)
                call fprank(q,ff,ncof,iband,ncc,sigma,c,sq,rank,wrk(la),wrk(lf),wrk(lh))
                q(1:ncof,1) = q(1:ncof,1)/dmax

                ! add to the sum of squared residuals, the contribution of reducing the rank.
                fp = fp+sq
            else

                ! backward substitution in case of full rank.
                c(:ncof) = fpback(a,f,ncof,iband,ncc)
                rank = ncof
                q(1:ncof,1) = a(1:ncof,1)/dmax

            endif

            ! find the coefficients in the standard b-spline representation of the spherical spline.
            rank_out = rank
          call fprpsp(nt,np,coco,cosi,c,ff,ncoff)

            ! Rank reduction can lose accuracy in the accumulated residual correction.
            ! Recompute the actual objective before selecting knots or accepting a fit.
            if (rank < ncof) fp = weighted_sphere_residual(teta,phi,r,w,tt(:nt), &
                tp(:np),c(:ncoff),spt,spp)

            ! test whether the least-squares spline is an acceptable solution.
            fpms = fp-s
            if (iopt<0 .or. abs(fpms)<=acc) then
              if (fp<=zero) then
                  ier = fitpack_interpolating_ok
                  fp  = zero
              else
                  if(ncof/=rank) ier = -rank
              endif
              return
            endif

            ! if f(p=inf) < s, accept the choice of knots.
            if (fpms<zero) exit compute_knots

            ! test whether we cannot further increase the number of knots.
            if (m<ncof) then
                ier = fitpack_too_many_knots
                return
            end if

            ! search where to add a new knot.
! find for each interval the sum of squared residuals fpint for the data points having the
! coordinate belonging to that knot interval. calculate also coord which is the same sum,
            ! weighted by the position of the data points considered.
            fpint(:nrint) = zero
            coord(:nrint) = zero

            do num=1,nreg
                num1 = num-1
                lt   = num1/npp
                l1   = lt+1
                lp   = num1-lt*npp
                l2   = lp+1+ntt
                jrot = lt*np4+lp
                in   = index(num)

                do while (in/=0)
                  store = zero
                  i1 = jrot
                  do i=1,4
                     store = store + spt(in,i)*dot_product(spp(in,1:4),c(i1+1:i1+4))
                     i1 = i1+np4
                  end do

                  store = (w(in)*(r(in)-store))**2
                  fpint(l1) = fpint(l1)+store
                  coord(l1) = coord(l1)+store*teta(in)
                  fpint(l2) = fpint(l2)+store
                  coord(l2) = coord(l2)+store*phi(in)
                  in = nummer(in)
                end do
            end do

            ! find the interval for which fpint is maximal on the condition that
            ! there still can be added a knot.
            l1 = merge(ntt+1,    1,ntest<nt+1)
            l2 = merge(ntt  ,nrint,npest<np+2)

            ! test whether we cannot further increase the number of knots.
            if (l1>l2) then
                ier = fitpack_insufficient_storage
                return
            end if

            add_knot: do

                i     = maxloc(fpint(l1:l2),1,fpint(l1:l2)>zero)
                l     = l1+i-1

                if (i==0) then
                    ier = fitpack_overlapping_knots
                    return
                end if

                fpmax = fpint(l)

                ! calculate the position of the new knot.
                arg = coord(l)/fpint(l)

                ! test in what direction the new knot is going to be added.
                choose_direction: if (l<=ntt) then

                   ! addition in the teta-direction
                   l4 = l+4
                   fpint(l) = zero
                   fac1 = tt(l4)-arg
                   fac2 = arg-tt(l4-1)

                else choose_direction

                   ! addition in the phi-direction
                   l4 = l+4-ntt
                   if (arg>=pi) then
                      arg = arg-pi
                      l4 = l4-nrr
                   endif

                   fpint(l) = zero
                   fac1 = tp(l4)-arg
                   fac2 = arg-tp(l4-1)

                endif choose_direction

                ! suitable location
                if (fac1<=(ten*fac2) .and. fac2<=(ten*fac1)) then

                    ! place knot
                    place_knot: if (l<=ntt) then

                        ! place in teta coordinate
                        j = nt
                        do i=l4,nt
                           tt(j+1) = tt(j)
                           j = j-1
                        end do
                        tt(l4) = arg
                        nt = nt+1

                    else place_knot

                        ! place in phi coordinate
                        ll = nrr+4
                        j = ll
                        do i=l4,ll
                           tp(j+1) = tp(j)
                           j = j-1
                        end do
                        tp(l4) = arg
                        np = np+2
                        nrr = nrr+1
                        do i=5,ll
                           j = i+nrr
                           tp(j) = tp(i)+pi
                        end do

                    endif place_knot

                    exit add_knot

                end if

            end do add_knot

      !  restart the computations with the new set of knots.
      end do compute_knots

! **********************************************************************************************
! **************
      ! part 2: determination of the smoothing spherical spline.
! **********************************************************************************************
! **************
! we have determined the number of knots and their position. we now compute the coefficients of
! the
! smoothing spline sp(teta,phi). the observation matrix a is extended by the rows of a matrix,
! expres-
! sing that sp(teta,phi) must be a constant function in the variable phi and a cubic polynomial
! in the
! variable teta. the corresponding weights of these additional rows are set to 1/(p).
! iteratively we then have
! to determine the value of p such that f(p) = sum((w(i)*(z(i)-sp(teta(i),phi(i))))**2)  be = s.
! we already know that the least-squares polynomial corresponds to p=0, and that the least-
! squares
! spherical spline corresponds to p=infinity. the iteration process makes use of rational
! interpolation.
! since f(p) is a convex and strictly decreasing function of p, it can be approximated by a
! rational
! function of the form r(p) = (u*p+v)/(p+w). three values of p (p1,p2,p3) with corresponding
! values of
! f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s) are used to calculate the new value of p such that
! r(p)=s.
      ! convergence is guaranteed by taking f1>0,f3<zero
! **********************************************************************************************
! **************

      ! evaluate the discontinuity jumps of the 3-th order derivative of the b-splines at the knots
      call fpdisc(tt,nt,degree_5,bt,ntest) ! tt(l),l=5,...,nt-4.
      call fpdisc(tp,np,degree_5,bp,npest) ! tp(l),l=5,...,np-4.

      !  initial value for p.
      p1 = zero
      f1 = sup-s
      p3 = -one
      f3 = fpms
      p  = sum(a(1:ncof,1))
      rn = ncof
      p  = rn/p

      ! find the bandwidth of the extended observation matrix.
      iband4 = merge(ncof,iband+3,ntt<=4)
      iband3 = iband4 -1
      check1 = fp_false
      check3 = fp_false

      ! iteration process to find the root of f(p)=s.
      iterations: do iter=1,maxit

          pinv = one/p

          ! store the triangularized observation matrix into q.
          ff(1:ncof) = f(1:ncof)
          q(1:ncof,1:iband4) = zero
          q(1:ncof,1:iband)  = a(1:ncof,1:iband)

          ! extend the observation matrix with the rows of a matrix, expressing
          ! that for teta=cst. sp(teta,phi) must be a constant function.
          nt6 = nt-6
          const_teta: do i=5,np4
              ii = i-4
              row(1:npp) = zero
              ll = ii
              do l=1,5
                 if (ll>npp) ll=1
                 row(ll) = row(ll)+bp(ii,l)
                 ll = ll+1
              end do
              facc = dot_product(row(1:npp),coco(1:npp))
              facs = dot_product(row(1:npp),cosi(1:npp))

              do j=1,nt6
                  ! initialize the new row.
                  h(1:iband) = zero

                  ! fill in the non-zero elements of the row. jrot records the column
                  ! number of the first non-zero element in the row.
                  jrot = 4+(j-2)*npp
                  if (j>1 .and. j<nt6) then
                      h(1:npp) = row(1:npp)
                  else
                      h(1) = facc
                      h(2) = facs
                      if(j==1) jrot = 2
                  endif

                  h(1:iband) = h(1:iband)*pinv
                  ri = zero

                  ! rotate the new row into triangle by givens transformations.
                  call fp_rotate_shifted(h, iband, q, ncc, ri, ff, jrot, ncof)
              end do
          end do const_teta

          ! extend the observation matrix with the rows of a matrix expressing
          ! that for phi=cst. sp(teta,phi) must be a cubic polynomial.
          const_phi: do i=5,nt4
              ii = i-4
              do j=1,npp
                  ! initialize the new row
                  h(1:iband4) = zero

                  ! fill in the non-zero elements of the row. jrot records the column
                  ! number of the first non-zero element in the row.
                  j1 = 1
                  fill_in: do l=1,5
                      il = ii+l
                      ij = npp

                      if (il==3 .or. il==nt4) then
                          if (il==3) then
                              j1 = 1
                              j2 = 2
                              ij = j+2
                          else
                              j1 = j1+3-j
                              j2 = j1-2
                              ij = 0
                          endif
                          h(j2:j2+1) = bt(ii,l)*[coco(j),cosi(j)]
                      endif
                      h(j1) = h(j1)+bt(ii,l)
                      j1 = j1+ij
                  end do fill_in

                  h(:iband4) = pinv*h(:iband4)
                  ri = zero
                  jrot = merge(3+j+(ii-3)*npp,1,ii>2)

                  ! rotate the new row into triangle by givens transformations.
                  call fp_rotate_shifted(h, iband4, q, ncc, ri, ff, jrot, ncof)

              end do
          end do const_phi

          ! find dmax, the maximum value for the diagonal elements in the reduced triangle.
          dmax  = max(zero,maxval(q(1:ncof,1)))

          ! check whether the matrix is rank deficient.
          sigma = eps*dmax

          if (all(q(1:ncof,1)>sigma)) then

              ! backward substitution in case of full rank.
              c(:ncof) = fpback(q,ff,ncof,iband4,ncc)
              rank = ncof

          else

              ! in case of rank deficiency, find the minimum norm solution
              lwest = ncof*iband4+ncof+iband4
              if (lwest>lwrk) then
                 ier = lwest
                 return
              end if
              lf = 1
              lh = lf+ncof
              la = lh+iband4
              call fprank(q,ff,ncof,iband4,ncc,sigma,c,sq,rank,wrk(la),wrk(lf),wrk(lh))

          end if

          q(1:ncof,1) = q(1:ncof,1)/dmax

          ! find the coefficients in the standard b-spline representation of the spherical spline.
          rank_out = rank
          call fprpsp(nt,np,coco,cosi,c,ff,ncoff)

          ! compute f(p).
          fp = zero
          get_fp: do num = 1,nreg
              num1 = num-1
              lt   = num1/npp
              lp   = num1-lt*npp
              jrot = lt*np4+lp
              in   = index(num)

              do while (in/=0)
                  store = zero
                  i1 = jrot
                  do i=1,4
                     store = store+spt(in,i)*dot_product(spp(in,1:4),c(i1+1:i1+4))
                     i1 = i1+np4
                  end do
                  fp = fp+(w(in)*(r(in)-store))**2
                  in = nummer(in)
              end do
          end do get_fp

          ! test whether the approximation sp(teta,phi) is an acceptable solution
          fpms = fp-s
          if (abs(fpms)<=acc) then
              if (ncof/=rank) ier = -rank
              return
          end if

          ! find the new value of p and carry out one more step.
          call root_finding_iterate(p1,f1,p2,f2,p3,f3,p,fpms,acc,check1,check3,success)
          if (.not.success) then
             ier = fitpack_s_too_small
             return
          end if

      end do iterations

      ! the maximum allowable number of iterations has been reached.
      ier = fitpack_maxit
      return

      end subroutine fpsphe

      pure subroutine fp_rotate_row(h, band, a, yi, z, j_start)
        implicit none
        real(dp),    intent(inout) :: h(:) !! Active row eliminated into the triangular factor.
        integer, intent(in) :: band !! Number of active factor diagonals.
        real(dp),    intent(inout) :: a(:,:) !! Banded triangular factor updated by rotations.
        real(dp),    intent(inout) :: yi !! Right-hand-side value eliminated by plane rotations.
        real(dp),    intent(inout) :: z(:) !! Transformed right-hand-side vector.
        integer, intent(in) :: j_start !! First factor row offset for the plane rotations.

          real(dp) :: cos, sin, piv
          integer :: i, j

          j = j_start
          do i = 1, band
              j = j + 1
              piv = h(i)
              if (equal(piv, zero)) cycle

              ! calculate parameters of givens transformation
              call fpgivs(piv, a(j,1), cos, sin)

              ! apply to right hand side
              call fprota(cos, sin, yi, z(j))

              ! apply to remaining row elements and matrix
              if (i < band) call fprota(cos, sin, h(i+1:band), a(j, 2:band-i+1))
          end do
      end subroutine fp_rotate_row

      pure subroutine fp_rotate_shifted(h, band, a, nest, yi, z, j_start, j_end)
        implicit none
        integer, intent(in) :: band !! Number of active factor diagonals.
        integer, intent(in) :: nest !! Leading factor storage dimension.
        integer, intent(in) :: j_start !! First factor row offset for the plane rotations.
        integer, intent(in) :: j_end !! Last factor row for the plane rotations.
        real(dp),    intent(inout) :: &
            h(band) !! Active row shifted while being eliminated into the factor.
        real(dp),    intent(inout) :: a(nest,*) !! Banded triangular factor updated by rotations.
        real(dp),    intent(inout) :: yi !! Right-hand-side value eliminated by plane rotations.
        real(dp),    intent(inout) :: z(*) !! Transformed right-hand-side vector.

          real(dp) :: cos, sin, piv
          integer :: noff, j

          do j = j_start, j_end
              piv = h(1)
              noff = min(j_end - j, band - 1)
              if (equal(piv, zero)) then
                  if (noff <= 0) return
              else
                  call fpgivs(piv, a(j,1), cos, sin)
                  call fprota(cos, sin, yi, z(j))
                  if (noff == 0) return
                  call fprota(cos, sin, h(2:noff+1), a(j, 2:noff+1))
              endif
              h(1:noff+1) = [h(2:noff+1), zero]
          end do
      end subroutine fp_rotate_shifted

      pure subroutine fporde(x,y,m,kx,ky,tx,nx,ty,ny,nummer,index,nreg)
        implicit none
      !  ..scalar arguments..
        integer, intent(in) :: m !! Number of scattered observations.
        integer, intent(in) :: kx !! Colatitude B-spline degree.
        integer, intent(in) :: ky !! Longitude B-spline degree.
        integer, intent(in) :: nx !! Full colatitude knot count.
        integer, intent(in) :: ny !! Full longitude knot count.
        integer, intent(in) :: nreg !! Allocated panel index capacity.
      !  ..array arguments..
        real(dp), intent(in) :: x(m) !! Scattered colatitudes in radians.
        real(dp), intent(in) :: y(m) !! Scattered longitudes in radians.
        real(dp), intent(in) :: tx(nx) !! Full colatitude knot vector.
        real(dp), intent(in) :: ty(ny) !! Full longitude knot vector.
        integer, intent(out) :: nummer(m) !! Next observation in each coordinate panel.
        integer, intent(out) :: index(nreg) !! Head observation index for each coordinate panel.
      !  ..local scalars..
      real(dp) :: xi,yi
      integer :: im,k,kx1,ky1,k1,l,l1,nk1x,nk1y,num,nyy

      !  ..
      kx1  = kx+1
      ky1  = ky+1
      nk1x = nx-kx1
      nk1y = ny-ky1
      nyy  = nk1y-ky
      index = 0

      points: do im=1,m
         xi = x(im)
         yi = y(im)
         l = fp_knot_interval(tx, xi, kx1, nk1x)
         l1 = l + 1

         k = fp_knot_interval(ty, yi, ky1, nk1y)
         k1 = k + 1
         num = (l-kx1)*nyy+k-ky
         nummer(im) = index(num)
         index(num) = im
      end do points
      return
      end subroutine fporde

      pure subroutine fprank(a,f,n,m,na,tol,c,sq,rank,aa,ff,h)
        implicit none
      !  ..scalar arguments..
        integer, intent(in) :: n !! Number of columns in the reduced linear system.
        integer, intent(in) :: m !! Bandwidth of the reduced triangular factor.
        integer, intent(in) :: na !! Leading dimension of the banded factor.
        integer, intent(out) :: rank !! Estimated factor rank after diagonal thresholding.
        real(dp), intent(in) :: tol !! Nonnegative diagonal threshold for rank reduction.
        real(dp), intent(out) :: sq !! Additional squared residual from rank reduction.

      !  ..array arguments..
        real(dp), intent(out) :: c(n) !! Minimum-norm reduced coefficient solution.
        real(dp), intent(inout) :: &
            a(na,m) !! Banded triangular factor, overwritten during rank reduction.
                                            ! triangularization by givens transformations.
        real(dp), intent(inout) :: &
            f(n) !! Transformed observations, overwritten during rank reduction.

        real(dp), intent(inout) :: aa(n,m) !! Rank-reduction factor workspace.
        real(dp), intent(inout) :: ff(n) !! Rank-reduction right-hand-side workspace.
        real(dp), intent(inout) :: h(m) !! Active factor row workspace.
      !  ..local scalars..
      integer     :: i,ii,ij,i1,i2,j,jj,j1,j2,j3,k,kk,m1,nl
      real(dp) :: cos,fac,piv,sin,yi
      real(dp) :: store,stor1,stor2,stor3

      !  ..
      m1 = m-1
! the rank deficiency nl is considered to be the number of sufficient small diagonal elements of
! a.
      nl = 0
      sq = zero
      rows: do i=1,n
         if (a(i,1)>tol) cycle rows

! if a sufficient small diagonal element is found, we put it to zero. the remainder of the row
! corresponding to that zero diagonal element is then rotated into triangle by givens rotations
! .
         ! the rank deficiency is increased by one.
         nl = nl+1
         if (i==n) cycle rows
         yi     = f(i)
         h(1:m) = [a(i,2:m),zero]
         i1 = i+1
         in_row: do ii=i1,n
            i2 = min(n-ii,m1)
            piv = h(1)
            if (equal(piv,zero)) then
               if (i2==0) exit in_row
               h(1:i2) = h(2:i2+1)
            else
               call fpgivs(piv,a(ii,1),cos,sin)
               call fprota(cos,sin,yi,f(ii))
               if (i2==0) exit in_row
               do j=1,i2
                 j1 = j+1
                 call fprota(cos,sin,h(j1),a(ii,j1))
                 h(j) = h(j1)
               end do
            endif
            h(i2+1) = zero
         end do in_row
         ! add to the sum of squared residuals the contribution of deleting
         ! the row with small diagonal element.
         sq = sq+yi**2
      end do rows

      !  rank denotes the rank of a.
      rank = n-nl

! let b denote the (rank*n) upper trapezoidal matrix which can be obtained from the (n*n) upper
! triangular matrix a by deleting the rows and interchanging the columns corresponding to a zero
! diagonal element. if this matrix is factorized using givens transformations as  b = (r) (u)
! where
      !    r is a (rank*rank) upper triangular matrix,
      !    u is a (rank*n) orthonormal matrix
! then the minimal least-squares solution c is given by c = b' v, where v is the solution of the
! system  (r) (r)' v = g  and g denotes the vector obtained from the old right hand side f, by
      !  removing the elements corresponding to a zero diagonal element of a.

      ! initialization.
      aa(1:rank,1:m) = zero

! form in aa the upper triangular matrix obtained from a by removing rows and columns with zero
! diagonal elements. form in ff the new right hand side by removing the elements of the old
! right
      !  hand side corresponding to a deleted row.
      ii = 0
      make_u: do i=1,n
         if (a(i,1)<=tol) cycle make_u
         ii = ii+1
         ff(ii) = f(i)
         aa(ii,1) = a(i,1)
         jj = ii
         kk = 1
         j  = i
         j1 = min(j-1,m1)
         if (j1/=0) then  ! probably unnecessary
            do k=1,j1
              j = j-1
              if (a(j,1)>tol) then
                 kk = kk+1
                 jj = jj-1
                 aa(jj,kk) = a(j,k+1)
              endif
            end do
         endif
      end do make_u

      !  form successively in h the columns of a with a zero diagonal element.
      ii = 0
      make_h: do i=1,n
        ii = ii+1
        if (a(i,1)>tol) cycle make_h
        ii = ii-1
        if (ii==0) cycle make_h
        jj = 1
        j  = i
        j1 = min(j-1,m1)
        do k=1,j1
           j = j-1
           if (a(j,1)>tol) then
             h(jj) = a(j,k+1)
             jj = jj+1
           endif
        end do
        h(jj:m) = zero

        ! rotate this column into aa by givens transformations.
        jj = ii
        rotate_col: do i1=1,ii
           j1 = min(jj-1,m1)
           piv = h(1)
           if (not_equal(piv,zero)) then
              call fpgivs(piv,aa(jj,1),cos,sin)
              if (j1==0) cycle make_h
              kk = jj
              do j2=1,j1
                 j3 = j2+1
                 kk = kk-1
                 call fprota(cos,sin,h(j3),aa(kk,j3))
                 h(j2) = h(j3)
              end do
           else
              if (j1==0) cycle make_h
              do j2=1,j1
                j3 = j2+1
                h(j2) = h(j3)
              end do
           endif
           jj = jj-1
           h(j3) = zero
        end do rotate_col
      end do make_h

      !  solve the system (aa) (f1) = ff
      ff(rank) = ff(rank)/aa(rank,1)
      i = rank-1
      if (i/=0) then
         do j=2,rank
            store = ff(i)
            i1 = min0(j-1,m1)
            k = i
            do ii=1,i1
               k = k+1
               stor1 = ff(k)
               stor2 = aa(i,ii+1)
               store = store-stor1*stor2
            end do
            stor1 = aa(i,1)
            ff(i) = store/stor1
            i = i-1
         end do
      endif

      !  solve the system  (aa)' (f2) = f1
      ff(1) = ff(1)/aa(1,1)
      if (rank>1) then
          do j=2,rank
            store = ff(j)
            i1 = min0(j-1,m1)
            k = j
            do ii=1,i1
               k = k-1
               stor1 = ff(k)
               stor2 = aa(k,ii+1)
               store = store-stor1*stor2
            end do
            stor1 = aa(j,1)
            ff(j) = store/stor1
         end do
      endif

      !  premultiply f2 by the transpoze of a.
      k = 0
      do i=1,n
        store = zero
        if (a(i,1)>tol) k = k+1
        j1 = min(i,m)
        kk = k
        ij = i+1
        do j=1,j1
           ij = ij-1
           if (a(ij,1)>tol) then
              stor1 = a(ij,j)
              stor2 = ff(kk)
              store = store+stor1*stor2
              kk = kk-1
           endif
        end do
        c(i) = store
      end do

      !  add to the sum of squared residuals the contribution of putting
      !  to zero the small diagonal elements of matrix (a).
      stor3 = zero
      add_to_res: do i=1,n
         if (a(i,1)<=tol) then
            store = f(i)
            i1 = min(n-i,m1)
            if (i1>0) then
               do j=1,i1
                  ij = i+j
                  stor1 = c(ij)
                  stor2 = a(i,j+1)
                  store = store-stor1*stor2
               end do
            endif
            fac = a(i,1)*c(i)
            stor1 = a(i,1)
            stor2 = c(i)
            stor1 = stor1*stor2
            stor3 = stor3+stor1*(stor1-store-store)
         endif
      end do add_to_res

      fac = stor3
      sq  = sq+fac
      return
      end subroutine fprank

      pure subroutine fprpsp(nt,np,co,si,c,f,ncoff)
        implicit none
          !  ..
          !  ..scalar arguments
        integer, intent(in) :: nt !! Active colatitude knot count.
        integer, intent(in) :: np !! Active longitude knot count.
        integer, intent(in) :: ncoff !! Number of physical tensor coefficients.

          !  ..array arguments
        real(dp), intent(in) :: co(np) !! Periodic cosine spline coefficients.
        real(dp), intent(in) :: si(np) !! Periodic sine spline coefficients.
        real(dp), intent(inout) :: &
            c(ncoff) !! Reduced coefficients, replaced by physical tensor coefficients.
        real(dp), intent(inout) :: f(ncoff) !! Physical tensor coefficient workspace.

          !  ..local scalars
          real(dp) :: cn,c1,c2,c3
          integer :: i,ii,j,k,l,ncof,npp,np4,nt4
          !  ..
          nt4 = nt-4
          np4 = np-4
          npp = np4-3
          ncof = 6+npp*(nt4-4)
          c1 = c(1)
          cn = c(ncof)
          f(1:np4)             = c1
          f(ncoff-np4+1:ncoff) = cn

          i = np4
          j = 1
          do l=3,nt4
             ii = i
             if (l==3 .or. l==nt4) then
                if (l==nt4) c1 = cn
                c2 = c(j+1)
                c3 = c(j+2)
                j = j+2
                do k=1,npp
                   i = i+1
                   f(i) = c1+c2*co(k)+c3*si(k)
                end do
             else
                do k=1,npp
                   i = i+1
                   j = j+1
                   f(i) = c(j)
                end do
             endif

             do k=1,3
                ii = ii+1
                i = i+1
                f(i) = f(ii)
             end do
          end do

          c(1:ncoff) = f(1:ncoff)

      return
      end subroutine fprpsp
end module scifort_interpolate_sphere_scattered_core
