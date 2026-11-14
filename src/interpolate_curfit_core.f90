! SPDX-License-Identifier: MIT AND BSD-3-Clause
! Copyright (c) 2022, Federico Perini. All rights reserved.
! Copyright (c) 2026 SciFort contributors
!
! Scalar, periodic and parametric fitting adapted from perazz/fitpack src/fitpack_core.F90,
! commit 1907c7bad9f5c4557970e906bdd4909eb52b2d4c, based on Paul Dierckx's FITPACK.
! License: LICENSES/perazz-fitpack-BSD-3-Clause.txt. See CODE_PROVENANCE.md.
module scifort_interpolate_curfit_core
    use scifort_kinds, only: dp
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: fitpack_smooth_curve, fitpack_smooth_parametric, fitpack_periodic_curve
    integer, parameter :: max_order = 19, degree_5 = 5, ione = 1
    integer, parameter :: fitpack_ok = 0, fitpack_interpolating_ok = -1
    integer, parameter :: fitpack_leastsquares_ok = -2, fitpack_insufficient_storage = 1
    integer, parameter :: fitpack_s_too_small = 2, fitpack_maxit = 3
    real(dp), parameter :: zero = 0.0_dp, one = 1.0_dp, half = 0.5_dp
    logical, parameter :: fp_true = .true., fp_false = .false.
contains

    subroutine fitpack_periodic_curve(x,y,w,k,s,knots,coeff,status)
        implicit none
        real(dp), intent(in) :: x(:) !! Finite increasing coordinates, including period endpoint.
        real(dp), intent(in) :: y(:) !! Observations; the final entry is ignored.
        real(dp), intent(in), optional :: w(:) !! Positive weights, default one; last entry ignored.
        integer, intent(in) :: k !! Degree, one through five.
        real(dp), intent(in) :: s !! Nonnegative weighted residual target; zero interpolates.
        real(dp), allocatable, intent(out) :: knots(:) !! Full periodic knots on success.
        real(dp), allocatable, intent(out) :: coeff(:) !! Coefficients padded to knot count.
        integer, intent(out) :: status !! Zero success, three numerical or allocation failure.
        real(dp), allocatable :: xx(:),yy(:),ww(:),t(:),c(:),fpint(:),z(:),q(:,:)
        real(dp), allocatable :: a1(:,:),a2(:,:),b(:,:),g1(:,:),g2(:,:)
        integer, allocatable :: nrdata(:)
        real(dp) :: span,yscale,wscale,target,logtarget,fp,residual,basis(max_order+1)
        integer :: m,nest,n,nc,info,ier,i,l
        status=3
        m=size(x)
        if (m>huge(m)/64) return
        ! Extended periodic knots require exponent headroom beyond the data interval.
        if (maxval(abs(x))>huge(one)/real(4*(k+1),dp)) return
        span=x(m)-x(1)
        if (span<tiny(one)) return
        nest=m+2*k
        allocate(xx(m),yy(m),ww(m),t(nest),c(nest),fpint(nest),z(nest),q(m,k+1), &
            a1(nest,k+1),a2(nest,k),b(nest,k+2),g1(nest,k+2),g2(nest,k+1), &
            nrdata(nest),stat=info)
        if (info/=0) return
        xx=(x-x(1))/span
        if (any(xx(2:)<=xx(:m-1))) return
        if (minval(xx(2:)-xx(:m-1))<tiny(one)**(one/real(4*(k+1),dp))) return
        yscale=maxval(abs(y(:m-1)))
        if (yscale==zero) yscale=one
        yy(:m-1)=y(:m-1)/yscale
        yy(m)=yy(1)
        if (any(y(:m-1)/=zero .and. yy(:m-1)==zero)) return
        ww=one
        if (present(w)) ww(:m-1)=w(:m-1)
        wscale=maxval(ww(:m-1))
        ww(:m-1)=ww(:m-1)/wscale
        ww(m)=ww(1)
        if (any(ww<tiny(one)**0.125_dp)) return
        target=zero
        if (s>zero) then
            logtarget=log(s)-2*log(yscale)-2*log(wscale)
            if (logtarget>=log(real(m,dp)+one)) then
                target=real(m,dp)+one
            else
                if (logtarget<log(tiny(one))) return
                target=exp(logtarget)
            end if
        end if
        t=zero
        c=zero
        fpint=zero
        z=zero
        q=zero
        a1=zero
        a2=zero
        b=zero
        g1=zero
        g2=zero
        nrdata=0
        n=2*(k+1)
        fp=zero
        ier=0
        call fpperi(0,xx,yy,ww,m,k,target,nest,0.001_dp,20,k+1,k+2,n,t,c,fp, &
            fpint,z,a1,a2,b,g1,g2,q,nrdata,ier)
        if (ier>0) return
        nc=n-k-1
        if (.not.all(ieee_is_finite(c(:nc)))) return
        if (yscale>one) then
            if (maxval(abs(c(:nc)))>huge(one)/yscale) return
        end if
        ! Validate the actual physical knots and rescaled coefficients, not just work arrays.
        t(:n)=x(1)+t(:n)*span
        t(k+1)=x(1)
        t(n-k)=x(m)
        if (.not.all(ieee_is_finite(t(:n)))) return
        if (any(t(2:n)<=t(:n-1))) return
        if (minval(t(2:n)-t(:n-1))<tiny(one)) return
        c(:nc)=c(:nc)*yscale
        fp=zero
        l=k+1
        do i=1,m-1
            l=fp_knot_interval(t,x(i),l,nc)
            basis=fpbspl(t,n,k,x(i),l)
            residual=ww(i)*(dot_product(c(l-k:l)/yscale,basis(:k+1))-yy(i))
            if (abs(residual)>sqrt(huge(one)/real(m,dp))) return
            fp=fp+residual**2
        end do
        if (target>zero) then
            if (fp>1.0011_dp*target) return
        else
            if (fp>real(m,dp)*(128*epsilon(one))**2) return
        end if
        allocate(knots(n),coeff(n),stat=info)
        if (info/=0) then
            if (allocated(knots)) deallocate(knots)
            if (allocated(coeff)) deallocate(coeff)
            return
        end if
        knots=t(:n)
        coeff=zero
        coeff(:nc)=c(:nc)
        status=0
    end subroutine fitpack_periodic_curve

    subroutine fitpack_smooth_curve(x, y, w, k, s, bounds, knots, coeff, status)
        implicit none
        real(dp), intent(in) :: x(:) !! Finite strictly increasing sample coordinates.
        real(dp), intent(in) :: y(:) !! Finite sample values, length size(x).
        real(dp), intent(in), optional :: w(:) !! Positive residual multipliers, default one.
        integer, intent(in) :: k !! Polynomial degree, from one through five.
        real(dp), intent(in) :: s !! Positive weighted squared-residual target.
        real(dp), intent(in) :: bounds(2) !! Finite enclosing left and right bounds.
        real(dp), allocatable, intent(out) :: knots(:) !! Full active knots on success.
        real(dp), allocatable, intent(out) :: coeff(:) !! Coefficients padded to knot count.
        integer, intent(out) :: status !! Zero success, three numerical/allocation failure.
        real(dp), allocatable :: xx(:), yy(:), ww(:), t(:), c(:), a(:, :), b(:, :), g(:, :)
        real(dp), allocatable :: q(:, :), z(:), fpint(:)
        integer, allocatable :: nrdata(:)
        real(dp) :: ax, lo, span, yscale, wscale, target, logtarget, fp, residual
        real(dp) :: basis(max_order+1)
        integer :: m, nest, n, ier, info, nc, i, l
        status = 3
        m = size(x)
        if (m > huge(m)/32) return
        nest = m + k + 1
        allocate (xx(m), yy(m), ww(m), t(nest), c(nest), a(nest, k + 1), b(nest, k + 2), &
                  g(nest, k + 2), q(m, k + 1), z(nest), fpint(nest), nrdata(nest), stat=info)
        if (info /= 0) return
        ax = maxval(abs(bounds))
        if (ax == zero) return
        lo = bounds(1)/ax
        span = bounds(2)/ax - lo
        if (span <= zero) return
        xx = (x/ax - lo)/span
        if (any(xx(2:) <= xx(:m - 1))) return
        ! Reserve exponent range for degree-dependent derivative-jump products.
        if (minval(xx(2:)-xx(:m-1)) < tiny(one)**(one/real(4*(k+1),dp))) return
        yscale = maxval(abs(y))
        if (yscale == zero) yscale = one
        yy = y/yscale
        if (any(y /= zero .and. yy == zero)) return
        ww = one
        if (present(w)) ww = w
        wscale = maxval(ww)
        ww = ww/wscale
        if (any(ww < tiny(one)**0.125_dp)) return
        logtarget = log(s) - 2*log(yscale) - 2*log(wscale)
        ! A normalized polynomial residual is at most m; avoid huge targets.
        if (logtarget >= log(real(m, dp) + one)) then
            target = real(m, dp) + one
        else
            ! An unrepresentable positive target is an explicit numerical failure.
            if (logtarget < log(tiny(one))) return
            target = exp(logtarget)
        end if
        t = zero
        c = zero
        a = zero
        b = zero
        g = zero
        q = zero
        z = zero
        fpint = zero
        nrdata = 0
        n = 2*(k + 1)
        call fpcurf(0, xx, yy, ww, m, zero, one, k, target, nest, 0.001_dp, 20, k + 1, k + 2, &
                    n, t, c, fp, fpint, z, a, b, g, q, nrdata, ier)
        if (ier > 0) return
        nc = n - k - 1
        if (.not. all(ieee_is_finite(c(:nc)))) return
        ! Check the returned representation, including interpolation-roundoff limits.
        fp=zero
        l=k+1
        do i=1,m
            l=fp_knot_interval(t,xx(i),l,nc)
            basis=fpbspl(t,n,k,xx(i),l)
            residual=ww(i)*(dot_product(c(l-k:l),basis(:k+1))-yy(i))
            if (abs(residual)>sqrt(huge(one)/real(m,dp))) return
            fp=fp+residual**2
        end do
        if (fp>1.001_dp*target) return
        if (yscale > one) then
            if (maxval(abs(c(:nc))) > huge(one)/yscale) return
        end if
        c(:nc) = c(:nc)*yscale
        allocate (knots(n), coeff(n), stat=info)
        if (info /= 0) then
            if (allocated(knots)) deallocate (knots)
            if (allocated(coeff)) deallocate (coeff)
            return
        end if
        knots = (one - t(:n))*bounds(1) + t(:n)*bounds(2)
        knots(:k + 1) = bounds(1)
        knots(n - k:) = bounds(2)
        coeff = zero
        coeff(:nc) = c(:nc)
        if (.not. all(ieee_is_finite(knots))) then
            deallocate (knots, coeff)
            return
        end if
        status = 0
    end subroutine fitpack_smooth_curve

    pure subroutine fpcurf(iopt, x, y, w, m, xb, xe, k, s, nest, tol, &
                           maxit, k1, k2, n, t, c, fp, fpint, z, a, b, g, q, nrdata, ier)
        implicit none

        !  ..
        !  ..scalar arguments..
        real(dp), intent(in) :: xb !! Left bound of the normalized approximation interval.
        real(dp), intent(in) :: xe !! Right bound of the normalized approximation interval.
        real(dp), intent(in) :: s !! Nonnegative normalized weighted squared-residual target.
        real(dp), intent(in) :: tol !! Relative residual-target tolerance, 0.001.
        real(dp), intent(out) :: fp !! Achieved normalized weighted squared residual.
        integer, intent(in) :: iopt !! FITPACK task selector; this wrapper always uses zero.
        integer, intent(in) :: m !! Number of strictly increasing samples.
        integer, intent(in) :: k !! Spline degree, from one through five.
        integer, intent(in) :: nest !! Allocated knot and coefficient capacity.
        integer, intent(in) :: maxit !! Maximum rational smoothing-parameter iterations, 20.
        integer, intent(in) :: k1 !! Degree plus one.
        integer, intent(in) :: k2 !! Degree plus two.
        integer, intent(inout) :: n !! Active knot count, updated during automatic refinement.
        integer, intent(out) :: ier !! FITPACK success codes <=0, positive numerical failure.

        !  ..array arguments..
        real(dp), intent(in) :: x(m) !! Strictly increasing normalized sample coordinates.
        real(dp), intent(in) :: y(m) !! Normalized scalar observations.
        real(dp), intent(in) :: w(m) !! Positive normalized residual multipliers.
        real(dp), intent(inout) :: t(nest) !! Full knot vector with repeated boundary knots.
        real(dp), intent(inout) :: c(nest) !! B-spline coefficients; active length n-k-1.
        real(dp), intent(inout) :: &
            fpint(nest) !! Residual sums per interval and saved knot-pass diagnostics.
        real(dp), intent(inout) :: z(nest) !! Triangularized observation right-hand side.
        real(dp), intent(inout) :: &
            a(nest, k1) !! Upper triangular observation factor in band storage.
        real(dp), intent(inout) :: &
            b(nest, k2) !! Derivative-jump rows used by the smoothing penalty.
        real(dp), intent(inout) :: g(nest, k2) !! Penalized triangular band factor workspace.
        real(dp), intent(inout) :: q(m, k1) !! Nonzero basis values per observation.
        integer, intent(inout) :: nrdata(nest) !! Interior sample counts per knot interval.
        !  ..local scalars..
        real(dp) :: acc, fpart, fpms, fpold, fp0, f1, f2, f3, p, pinv, p1, p2, p3, rn, store, &
                    term, wi, xi, yi
        integer :: i, it, iter, j, k3, l, l0, mk1, nk1, nmax, nmin, nplus, npl1, nrint, n8
        !  ..local arrays..
        real(dp) :: h(max_order + 1)
        logical :: new, check1, check3, success, accepted

        ier = fitpack_ok
        accepted = fp_false
        fpold = zero
        fp0 = zero
        nplus = 0
        fpms = huge(zero)
        nk1 = n - k1

        ! calculation of acc, the absolute tolerance for the root of f(p)=s.
        acc = tol*s

        ! determine nmax, the number of knots for spline interpolation.
        nmax = m + k1

        ! *****
        !  part 1: determination of the number of knots and their position
        ! *****
        ! given a set of knots we compute the least-squares spline sinf(x), and the corresponding
        ! sum
        !  of squared residuals fp=f(p=inf).
        !  if iopt=-1 sinf(x) is the requested approximation.
        !  if iopt=0 or iopt=1 we check whether we can accept the knots:
        !    if fp <=s we will continue with the current set of knots.
        !    if fp > s we will increase the number of knots and compute the corresponding least-
        !    squares spline until finally fp<=s.
        !  the initial choice of knots depends on the value of s and iopt.
        ! if s=0 we have spline interpolation; in that case the number of knots equals nmax = m+k+1.
        !    if (s>0 and iopt=0) we first compute the least-squares polynomial curve of degree k;
        !      n = nmin = 2*k+2
        !    iopt=1 we start with the set of knots found at the last call of the routine, except for
        ! the case that s > fp0; then we compute directly the least-squares polynomial of degree k.
        ! *****

        !  determine nmin, the number of knots for polynomial approximation.
        nmin = 2*k1
        bootstrap: if (iopt >= 0) then

            interpolating: if (s <= zero) then

                !  if s=0, s(x) is an interpolating spline.
                !  test whether the required storage space exceeds the available one.
                n = nmax
                if (nmax > nest) then
                    ier = fitpack_insufficient_storage
                    return
                end if

                !  find the position of the interior knots in case of interpolation.
                mk1 = m - k1
                if (mk1 /= 0) then
                    k3 = k/2
                    i = k2
                    j = k3 + 2
                    do l = 1, mk1
                        t(i) = merge(x(j), (x(j) + x(j - 1))*half, k3*2 /= k)
                        i = i + 1
                        j = j + 1
                    end do
                end if

                else interpolating

                !  if s>0 our initial choice of knots depends on the value of iopt.
                !  if iopt=0 or iopt=1 and s>=fp0, we start computing the least-squares
                !  polynomial of degree k which is a spline without interior knots.
                !  if iopt=1 and fp0>s we start computing the least squares spline
                !  according to the set of knots found at the last call of the routine.
                use_last_call: if (iopt /= 0 .and. n /= nmin) then
                    fp0 = fpint(n)
                    fpold = fpint(n - 1)
                    nplus = nrdata(n)
                    if (fp0 <= s) then
                        n = nmin
                        fpold = zero
                        nplus = 0
                        nrdata(1) = m - 2
                    end if
                    else use_last_call
                    n = nmin
                    fpold = zero
                    nplus = 0
                    nrdata(1) = m - 2
                end if use_last_call
            end if interpolating

        end if bootstrap

        !  main loop for the different sets of knots. m is a save upper bound
        !  for the number of trials.
        iter = 0
        main_loop: do while (iter <= m)

            iter = iter + 1

            if (n == nmin) ier = fitpack_leastsquares_ok

            ! find nrint, tne number of knot intervals.
            nrint = n - nmin + 1

            ! find the position of the additional knots which are needed for
            ! the b-spline representation of s(x).
            nk1 = n - k1
            t(1:k1) = xb
            t(nk1 + 1:n) = xe

            ! compute the b-spline coefficients of the least-squares spline
            ! sinf(x). the observation matrix a is built up row by row and
            ! reduced to upper triangular form by givens transformations.
            ! at the same time fp=f(p=inf) is computed.
            fp = zero

            ! initialize the observation matrix a.
            z(1:nk1) = zero
            a(1:nk1, 1:k1) = zero
            l = k1
            coefs: do it = 1, m

                ! fetch the current data point x(it),y(it).
                xi = x(it)
                wi = w(it)
                yi = y(it)*wi

                ! search for knot interval t(l) <= xi < t(l+1).
                l = fp_knot_interval(t, xi, l, nk1)

                ! evaluate the (k+1) non-zero b-splines at xi and store them in q.
                h = fpbspl(t, n, k, xi, l)

                q(it, 1:k1) = h(1:k1)
                h(:k1) = wi*h(:k1)

                ! rotate the new row of the observation matrix into triangle.
                call fp_rotate_row(h, k1, a, yi, z, l - k1)

                !  add contribution of this row to the sum of squares of residual
                !  right hand sides.
                fp = fp + yi*yi

            end do coefs

            if (ier == fitpack_leastsquares_ok) fp0 = fp
            fpint(n - 1:n) = [fpold, fp0]
            nrdata(n) = nplus

            ! backward substitution to obtain the b-spline coefficients.
            c(:nk1) = fpback(a, z, nk1, k1, nest)

            ! test whether the approximation sinf(x) is an acceptable solution. the three acceptance
            ! exits leave through the epilogue below, which screens the accepted factor once.
            if (iopt < 0) then
                accepted = fp_true
                exit main_loop
            end if

            fpms = fp - s
            if (abs(fpms) < acc) then
                accepted = fp_true
                exit main_loop
            end if

            ! if f(p=inf) < s accept the choice of knots and go on to the smoothing spline.
            if (fpms < zero) exit main_loop

            ! if n = nmax, sinf(x) is an interpolating spline.
            if (n == nmax) then
                ier = fitpack_interpolating_ok
                accepted = fp_true
                exit main_loop
            end if

            ! increase the number of knots.
            ! if n=nest we cannot increase the number of knots because of the storage capacity
            ! limitation.
            if (n == nest) then
                ier = fitpack_insufficient_storage
                return
            end if

            ! determine the number of knots nplus we are going to add.
            if (ier == fitpack_ok) then
                npl1 = nplus*2
                rn = nplus
                if (fpold - fp > acc) npl1 = fp_knots_to_add(rn, fpms, fpold - fp)
                nplus = min(nplus*2, max(npl1, nplus/2, 1))
            else
                nplus = 1
                ier = fitpack_ok
            end if

            ! initialize update
            fpold = fp

            ! compute the sum((w(i)*(y(i)-s(x(i))))**2) for each knot interval
            ! t(j+k) <= x(i) <= t(j+k+1) and store it in fpint(j),j=1,2,...nrint.
            fpart = zero
            i = 1
            l = k2
            new = .false.
            square_residuals: do it = 1, m
                if (x(it) >= t(l) .and. l <= nk1) then
                    new = .true.
                    l = l + 1
                end if
                l0 = l - k2
                term = dot_product(c(l0 + 1:l0 + k1), q(it, 1:k1))

                term = (w(it)*(term - y(it)))**2
                fpart = fpart + term

                if (new) then
                    store = term*half
                    fpint(i) = fpart - store
                    i = i + 1
                    fpart = store
                    new = .false.
                end if
            end do square_residuals

            fpint(nrint) = fpart

            add_new_knots: do l = 1, nplus

                ! add a new knot.
                call fpknot(x, m, t, n, fpint, nrdata, nrint, nest, ione)

                ! if n=nmax we locate the knots as for interpolation.
                if (n == nmax) then
                    mk1 = m - k1
                    if (mk1 /= 0) then
                        k3 = k/2
                        i = k2
                        j = k3 + 2
                        do l0 = 1, mk1
                            t(i) = merge(x(j), (x(j) + x(j - 1))*half, k3*2 /= k)
                            i = i + 1
                            j = j + 1
                        end do
                    end if

                    ! restart main loop
                    iter = 0
                    cycle main_loop

                end if

                ! test whether we cannot further increase the number of knots.
                if (n == nest) exit add_new_knots

            end do add_new_knots
            !  restart the computations with the new set of knots.
        end do main_loop

        !  sinf(x) is the answer whenever the knot pass accepted it, and also when the least-squares
        !  kth degree polynomial already solves our approximation problem. screen the factor that
        !  produced those coefficients: a zero diagonal means the back substitution had to leave a
        !  coefficient unconstrained, so the fit is not determined.
        if (accepted .or. ier == fitpack_leastsquares_ok) then
            if (fp_rank_deficient(a, nk1, nest)) ier = fitpack_s_too_small
            return
        end if

        ! *****
        !  part 2: determination of the smoothing spline sp(x).
        ! *****
        !  we have determined the number of knots and their position.
        ! we now compute the b-spline coefficients of the smoothing spline sp(x). the observation
        ! matrix a
        ! is extended by the rows of matrix b expressing that the kth derivative discontinuities of
        ! sp(x)
        !  at the interior knots t(k+2),...t(n-k-1) must be zero. the corresponding weights of these
        !  additional rows are set to 1/p.
        !  iteratively we then have to determine the value of p such that f(p), the sum of squared
        ! residuals be = s. we already know that the least squares kth degree polynomial corresponds
        ! to p=0, and that the least-squares spline corresponds to p=infinity. the iteration process
        ! which is proposed here, makes use of rational interpolation. since f(p) is a convex and
        ! strictly
        ! decreasing function of p, it can be approximated by a rational function r(p) =
        ! (u*p+v)/(p+w).
        ! three values of p(p1,p2,p3) with corresponding values of f(p)
        ! (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s)
        ! are used to calculate the new value of p such that r(p)=s. convergence is guaranteed by
        ! taking
        !  f1>0 and f3<zero
        ! *****

        !  evaluate the discontinuity jump of the kth derivative of the
        !  b-splines at the knots t(l),l=k+2,...n-k-1 and store in b.
        call fpdisc(t, n, k2, b, nest)

        !  initial value for p.
        p1 = zero
        f1 = fp0 - s
        p3 = -one
        f3 = fpms
        p = sum(a(1:nk1, 1))
        rn = nk1
        p = rn/p
        check1 = fp_false
        check3 = fp_false
        n8 = n - nmin

        !  iteration process to find the root of f(p) = s.
        iter = 0
        find_root: do while (iter < maxit)

            iter = iter + 1

            !  the rows of matrix b with weight 1/p are rotated into the
            !  triangularised observation matrix a which is stored in g.
            pinv = one/p
            c(1:nk1) = z(1:nk1)
            g(1:nk1, 1:k1) = a(1:nk1, 1:k1)
            g(1:nk1, k2) = zero
            b_rows: do it = 1, n8
                ! the row of matrix b is rotated into triangle by givens transformation
                h(1:k2) = b(it, 1:k2)*pinv
                yi = zero
                call fp_rotate_shifted(h, k2, g, nest, yi, c, it, nk1)
            end do b_rows

            ! backward substitution to obtain the b-spline coefficients.
            c(:nk1) = fpback(g, c, nk1, k2, nest)

            ! computation of f(p).
            fp = zero
            l = k2

            get_fp: do it = 1, m
                if (x(it) >= t(l) .and. l <= nk1) l = l + 1
                l0 = l - k2
                term = dot_product(c(l0 + 1:l0 + k1), q(it, :))
                fp = fp + (w(it)*(term - y(it)))**2
            end do get_fp

            ! success! the approximation sp(x) is an acceptable solution: screen its factor once.
            fpms = fp - s
            if (abs(fpms) < acc) then
                if (fp_rank_deficient(g, nk1, nest)) ier = fitpack_s_too_small
                return
            end if

            ! find the new value of p and carry out one more step.
            call root_finding_iterate(p1, f1, p2, f2, p3, f3, p, fpms, acc, check1, check3, success)
            if (.not. success) then
                ier = fitpack_s_too_small
                return
            end if

        end do find_root

        ! maximum number of iterations reached
        ier = fitpack_maxit
        return

    end subroutine fpcurf

    pure subroutine fpknot(x, m, t, n, fpint, nrdata, nrint, nest, istart)
        implicit none

        !  ..scalar arguments..
        integer, intent(in) :: m !! Number of strictly increasing samples.
        integer, intent(in) :: nest !! Allocated knot and coefficient capacity.
        integer, intent(in) :: istart !! First sample index for interval bookkeeping.
        integer, intent(inout) :: n !! Active knot count, updated during automatic refinement.
        integer, intent(inout) :: nrint !! Number of knot intervals, increased after insertion.
        !  ..array arguments..
        real(dp), intent(in) :: x(m) !! Strictly increasing normalized sample coordinates.
        real(dp), intent(inout) :: t(nest) !! Full knot vector with repeated boundary knots.
        real(dp), intent(inout) :: &
            fpint(nest) !! Residual sums per interval and saved knot-pass diagnostics.
        integer, intent(inout) :: nrdata(nest) !! Interior sample counts per knot interval.

        !  ..local scalars..
        real(dp) :: an, am, fpmax
        integer :: ihalf, j, jbegin, jj, jk, jpoint, k, maxbeg, maxpt, next, nrx, number
        !  ..
        number = 0
        maxpt = 0
        maxbeg = 0
        k = (n - nrint - 1)/2
        ! search for knot interval t(number+k) <= x <= t(number+k+1) where fpint(number) is maximal
        ! on
        ! the
        ! condition that nrdata(number)/=0. an interval is eligible only if it holds at least one
        ! interior
        ! data point (jpoint/=0); the first eligible interval seeds the search (number==0), after
        ! which any
        ! interval with a strictly larger residual takes over. seeding on the first eligible
        ! interval,
        ! rather
        ! than only on fpmax<fpint(j), guarantees a valid selection even when every residual is
        ! exactly zero
        ! (an all-zero fpint would otherwise leave number==0 and fall through to use of
        ! uninitialized
        ! am).
        fpmax = zero
        jbegin = istart
        do j = 1, nrint
            jpoint = nrdata(j)

            if (jpoint /= 0 .and. (number == 0 .or. fpint(j) > fpmax)) then
                fpmax = fpint(j)
                number = j
                maxpt = jpoint
                maxbeg = jbegin
            end if

            jbegin = jbegin + jpoint + 1
        end do

        ! no interval carries an interior data point: there is nothing to refine, so leave the knot
        ! set as is
        if (number == 0) return

        !  let coincide the new knot t(number+k+1) with a data point x(nrx)
        !  inside the old knot interval t(number+k) <= x <= t(number+k+1).
        ihalf = maxpt/2 + 1
        nrx = maxbeg + ihalf
        next = number + 1

        !  adjust the different parameters.
        if (next <= nrint) then
            do j = next, nrint
                jj = next + nrint - j
                fpint(jj + 1) = fpint(jj)
                nrdata(jj + 1) = nrdata(jj)
                jk = jj + k
                t(jk + 1) = t(jk)
            end do
        end if

        !  split the selected interval in two, apportioning its residual between the halves.
        !  number>0 is guaranteed by the early return above, so am is always initialized here.
        nrdata(number) = ihalf - 1
        nrdata(next) = maxpt - ihalf
        am = maxpt
        an = nrdata(number)
        fpint(number) = fpmax*an/am

        an = nrdata(next)
        fpint(next) = fpmax*an/am
        jk = next + k
        t(jk) = x(nrx)

        n = n + 1
        nrint = nrint + 1

    end subroutine fpknot

    elemental integer function fp_knots_to_add(rn, fpms, reduc) result(npl1)
        implicit none
        real(dp), intent(in) :: rn !! Previous number of knots added, represented as real.
        real(dp), intent(in) :: fpms !! Residual minus requested smoothing target.
        real(dp), intent(in) :: reduc !! Positive residual reduction from the preceding knot pass.

        real(dp), parameter :: npl_max = real(huge(0), dp)
        real(dp) :: ratio

        ratio = rn*fpms/reduc
        if (ratio < npl_max) then
            npl1 = int(max(ratio, zero))
        else
            npl1 = huge(0) ! saturate; a nan ratio lands here too
        end if

    end function fp_knots_to_add

    pure logical function fp_rank_deficient(a, n, nest) result(singular)
        implicit none
        integer, intent(in) :: n !! Number of active coefficient rows.
        integer, intent(in) :: nest !! Leading dimension of the band factor.
        real(dp), intent(in) :: &
            a(nest, *) !! Triangular band factor; its diagonal is stored in column one.
        integer :: i
        singular = fp_true
        do i = 1, n
            if (equal(a(i, 1), zero)) return
        end do
        singular = fp_false
    end function fp_rank_deficient

    pure subroutine fp_rotate_row(h, band, a, yi, z, j_start)
        implicit none
        real(dp), intent(inout) :: h(:) !! Working row to eliminate by Givens rotations.
        integer, intent(in) :: band !! Number of entries in the observation or penalty row.
        real(dp), intent(inout) :: a(:, :) !! Upper triangular band factor updated by rotations.
        real(dp), intent(inout) :: &
            yi !! Observation or penalty right-hand side updated by rotations.
        real(dp), intent(inout) :: z(:) !! Right-hand side updated along with the factor.
        integer, intent(in) :: j_start !! First factor row (minus one for unshifted rotations).

        real(dp) :: cos, sin, piv
        integer :: i, j

        j = j_start
        do i = 1, band
            j = j + 1
            piv = h(i)
            if (equal(piv, zero)) cycle

            ! calculate parameters of givens transformation
            call fpgivs(piv, a(j, 1), cos, sin)

            ! apply to right hand side
            call fprota(cos, sin, yi, z(j))

            ! apply to remaining row elements and matrix
            if (i < band) call fprota(cos, sin, h(i + 1:band), a(j, 2:band - i + 1))
        end do
    end subroutine fp_rotate_row

    pure subroutine fp_rotate_shifted(h, band, a, nest, yi, z, j_start, j_end)
        implicit none
        integer, intent(in) :: band !! Number of entries in the observation or penalty row.
        integer, intent(in) :: nest !! Allocated knot and coefficient capacity.
        integer, intent(in) :: j_start !! First factor row (minus one for unshifted rotations).
        integer, intent(in) :: j_end !! Last active coefficient row in the band factor.
        real(dp), intent(inout) :: h(band) !! Working row to eliminate by Givens rotations.
        real(dp), intent(inout) :: a(nest, *) !! Upper triangular band factor updated by rotations.
        real(dp), intent(inout) :: &
            yi !! Observation or penalty right-hand side updated by rotations.
        real(dp), intent(inout) :: z(*) !! Right-hand side updated along with the factor.

        real(dp) :: cos, sin, piv
        integer :: noff, j

        do j = j_start, j_end
            piv = h(1)
            noff = min(j_end - j, band - 1)
            if (equal(piv, zero)) then
                if (noff <= 0) return
            else
                call fpgivs(piv, a(j, 1), cos, sin)
                call fprota(cos, sin, yi, z(j))
                if (noff == 0) return
                call fprota(cos, sin, h(2:noff + 1), a(j, 2:noff + 1))
            end if
            h(1:noff + 1) = [h(2:noff + 1), zero]
        end do
    end subroutine fp_rotate_shifted

    elemental logical function equal(a, b) result(ok)
        implicit none
        real(dp), intent(in) :: a !! First value in the spacing-based comparison.
        real(dp), intent(in) :: b !! Second value in the spacing-based comparison.
        ok = abs(a - b) < spacing(merge(a, b, abs(a) < abs(b)))
    end function equal

    elemental logical function not_equal(a, b) result(ok)
        implicit none
        real(dp), intent(in) :: a !! First value in the spacing-based comparison.
        real(dp), intent(in) :: b !! Second value in the spacing-based comparison.
        ok = .not. equal(a, b)
    end function not_equal

    pure function fp_knot_interval(t, x, l_start, l_max) result(l)
        implicit none
        real(dp), intent(in) :: t(:) !! Nondecreasing B-spline knot vector.
        real(dp), intent(in) :: x !! Query coordinate in the active knot span.
        integer, intent(in) :: l_start !! Initial interval index for the monotone search.
        integer, intent(in) :: l_max !! Last permitted knot interval index.
        integer :: l

        integer, parameter :: linear_threshold = 8
        integer :: lo, hi, mid

        ! for small ranges, use linear search (also optimal for warm starts)
        if (l_max - l_start <= linear_threshold) then
            l = l_start
            do while (x >= t(l + 1) .and. l < l_max)
                l = l + 1
            end do
            return
        end if

        ! handle boundary cases
        if (x < t(l_start + 1)) then
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
            mid = (lo + hi)/2
            if (t(mid + 1) <= x) then
                lo = mid + 1
            else
                hi = mid
            end if
        end do
        l = lo
    end function fp_knot_interval

    pure function fpback(a, z, n, k, nest) result(c)
        implicit none

        !  ..scalar arguments..
        integer, intent(in) :: n !! Order of the upper-triangular system.
        integer, intent(in) :: k !! Number of diagonals stored in the upper band.
        integer, intent(in) :: nest !! Leading storage dimension of the triangular factor.
        !  ..array arguments..
        real(dp), intent(in) :: a(nest, k) !! Upper-triangular band factor, diagonal in column one.
        real(dp), intent(in) :: z(n) !! Right-hand side of the triangular solve.
        real(dp)              :: c(n)
        !  ..local scalars..
        real(dp) :: store
        integer :: i, i1, j, k1, l, m, jm1
        !  ..
        k1 = k - 1

        if (equal(a(n, 1), zero)) then
            c(n) = zero
        else
            c(n) = z(n)/a(n, 1)
        end if

        i = n - 1
        if (i > 0) then

            jm1 = 1
            rows: do j = 2, n
                store = z(i)
                i1 = merge(jm1, k1, j <= k1)
                m = i
                do l = 1, i1
                    m = m + 1
                    store = store - c(m)*a(i, l + 1)
                end do
                if (equal(a(i, 1), zero)) then
                    c(i) = zero
                else
                    c(i) = store/a(i, 1)
                end if
                i = i - 1
                jm1 = j
            end do rows

        end if

    end function fpback

    pure function fpbspl(t, n, k, x, l) result(h)
        implicit none
        integer, intent(in) :: n !! Length of the nondecreasing knot vector.
        integer, intent(in) :: k !! B-spline polynomial degree, at most max_order.
        integer, intent(in) :: l !! Knot interval index containing the query coordinate.
        real(dp), intent(in) :: x !! Query coordinate in the active knot span.
        real(dp), intent(in) :: t(n) !! Nondecreasing B-spline knot vector.
        real(dp) :: h(max_order + 1)

        ! local variables
        real(dp) :: f, hh(max_order + 1)
        integer :: i, j, li, lj

        h(1) = one
        do j = 1, k
            hh(1:j) = h(1:j)
            h(1) = zero
            do i = 1, j
                li = l + i
                lj = li - j
                if (not_equal(t(li), t(lj))) then
                    f = hh(i)/(t(li) - t(lj))
                    h(i) = h(i) + f*(t(li) - x)
                    h(i + 1) = f*(x - t(lj))
                else
                    h(i + 1) = zero
                end if
            end do
        end do

    end function fpbspl

    pure subroutine fpdisc(t, n, k2, b, nest)
        implicit none

        !  ..scalar arguments..
        integer, intent(in) :: n !! Length of the nondecreasing knot vector.
        integer, intent(in) :: k2 !! Polynomial degree plus two.
        integer, intent(in) :: nest !! Leading storage dimension of the triangular factor.
        !  ..array arguments..
        real(dp), intent(in) :: t(n) !! Nondecreasing B-spline knot vector.
        real(dp), intent(inout) :: b(nest, k2) !! Derivative-jump values at the interior knots.
        !  ..local scalars..
        real(dp) :: an, fac, prod
        integer :: i, ik, j, jk, k, k1, l, lj, lk, lmk, lp, nk1, nrint
        !  ..local array..
        real(dp) :: h(2*(degree_5 + 1))
        !  ..
        k1 = k2 - 1
        k = k1 - 1
        nk1 = n - k1
        nrint = nk1 - k
        an = nrint
        fac = an/(t(nk1 + 1) - t(k1))
        do l = k2, nk1
            lmk = l - k1
            do j = 1, k1
                ik = j + k1
                lj = l + j
                lk = lj - k2
                h(j) = t(l) - t(lk)
                h(ik) = t(l) - t(lj)
            end do
            lp = lmk
            do j = 1, k2
                jk = j
                prod = h(j)
                do i = 1, k
                    jk = jk + 1
                    prod = prod*h(jk)*fac
                end do
                lk = lp + k1
                b(lmk, j) = (t(lk) - t(lp))/prod
                lp = lp + 1
            end do
        end do
        return
    end subroutine fpdisc

    elemental subroutine fpgivs(piv, ww, cos, sin)
        implicit none
        real(dp), intent(in) :: piv !! Pivot to eliminate with a Givens rotation.
        real(dp), intent(inout) :: ww !! Diagonal entry updated by the Givens rotation.
        real(dp), intent(out) :: cos !! Cosine of the plane rotation.
        real(dp), intent(out) :: sin !! Sine of the plane rotation.
        !  ..local scalars..
        real(dp) :: dd, store

        store = abs(piv)
        !  evaluate only the taken branch: an if (not merge) avoids a spurious
        !  divide-by-zero in the unselected expression (ww/piv or piv/ww) when
        !  the other operand is zero, which would raise ieee_invalid under fpe traps.
        if (store >= ww) then
            dd = store*sqrt(one + (ww/piv)**2)
        else
            dd = ww*sqrt(one + (piv/ww)**2)
        end if
        cos = ww/dd
        sin = piv/dd
        ww = dd
        return
    end subroutine fpgivs

    elemental subroutine fprota(cos, sin, a, b)
        implicit none

        !  ..scalar arguments..
        real(dp), intent(in) :: cos !! Cosine of the plane rotation.
        real(dp), intent(in) :: sin !! Sine of the plane rotation.
        real(dp), intent(inout) :: a !! First component, overwritten by its rotation.
        real(dp), intent(inout) :: b !! Second component, overwritten by its rotation.

        ! ..local scalars..
        real(dp) :: stor1, stor2

        !  ..
        stor1 = a
        stor2 = b
        b = cos*stor2 + sin*stor1
        a = cos*stor1 - sin*stor2
        return

    end subroutine fprota

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
        if (.not. check3) then
            if ((f2 - f3) > acc) then
                check3 = f2 < zero
            else
                ! our initial choice of p is too large.
                p3 = p2
                f3 = f2
                p = p*con4
                if (p <= p1) p = p1*con9 + p2*con1
                return
            end if
        end if

        if (.not. check1) then
            if ((f1 - f2) > acc) then
                check1 = f2 > zero
            else
                ! our initial choice of p is too small
                p1 = p2
                f1 = f2
                p = p/con4
                if (p3 >= zero .and. p >= p3) p = p2*con1 + p3*con9
                return
            end if
        end if

        ! test whether the iteration process proceeds as theoretically expected.
        if (f2 >= f1 .or. f2 <= f3) then
            success = .false.
            return
        else
            ! find the new value of p.
            call fprati(p1, f1, p2, f2, p3, f3, p)
        end if

    end subroutine root_finding_iterate

    elemental subroutine fprati(p1, f1, p2, f2, p3, f3, p)
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
        real(dp) :: h1, h2, h3
        !  ..
        if (p3 > zero) then
            h1 = f1*(f2 - f3)
            h2 = f2*(f3 - f1)
            h3 = f3*(f1 - f2)
            p = -(p1*p2*h3 + p2*p3*h1 + p3*p1*h2)/(p1*h1 + p2*h2 + p3*h3)
        else
            !  value of p in case p3 = infinity.
            p = (p1*(f1 - f3)*f2 - p2*(f2 - f3)*f1)/((f1 - f2)*f3)
        end if

        !  adjust the value of p1,f1,p3 and f3 such that f1 > 0 and f3 < 0.
        if (f2 >= zero) then
            p1 = p2
            f1 = f2
        else
            p3 = p2
            f3 = f2
        end if

        return
    end subroutine fprati

    subroutine fitpack_smooth_parametric(u, x, w, k, s, bounds, knots, coeff, status)
        implicit none
        real(dp), intent(in) :: u(:) !! Finite strictly increasing parameter sites.
        real(dp), intent(in) :: x(:, :) !! Finite coordinates, shape (size(u),channels).
        real(dp), intent(in), optional :: w(:) !! Positive residual multipliers, default one.
        integer, intent(in) :: k !! Polynomial degree, from one through five.
        real(dp), intent(in) :: s !! Positive total weighted squared-residual target.
        real(dp), intent(in) :: bounds(2) !! Finite enclosing parameter bounds.
        real(dp), allocatable, intent(out) :: knots(:) !! Full active knots on success.
        real(dp), allocatable, intent(out) :: coeff(:, :) !! Unpadded coefficients by channel.
        integer, intent(out) :: status !! Zero success, three numerical/allocation failure.
        real(dp), allocatable :: uu(:), xx(:), ww(:), t(:), c(:), a(:, :), b(:, :), g(:, :)
        real(dp), allocatable :: q(:, :), z(:), fpint(:)
        integer, allocatable :: nrdata(:)
        real(dp) :: span, xscale, wscale, target, logtarget, fp, residual, limit
        real(dp) :: basis(max_order + 1)
        integer :: m, ndim, nest, ncc, n, ncoef, ier, info, i, j, l, offset
        status = 3
        m = size(u)
        ndim = size(x, 2)
        if (m > huge(m)/32) return
        nest = m + k + 1
        if (ndim > huge(m)/nest) return
        ncc = nest*ndim
        allocate (uu(m), xx(m*ndim), ww(m), t(nest), c(ncc), a(nest, k + 1), b(nest, k + 2), &
                  g(nest, k + 2), q(m, k + 1), z(ncc), fpint(nest), nrdata(nest), stat=info)
        if (info /= 0) return
        ! Subtract the origin before scaling to preserve large-offset spacing.
        ! A nonrepresentable physical span cannot be evaluated by the public spline API.
        if (bounds(1)<zero .and. bounds(2)>zero) then
            if (bounds(2)>huge(one)+bounds(1)) return
        end if
        span = bounds(2)-bounds(1)
        if (span <= zero) return
        uu = (u-bounds(1))/span
        if (minval(uu(2:) - uu(:m - 1)) < tiny(one)**(one/real(4*(k + 1), dp))) return
        xscale = maxval(abs(x))
        if (xscale == zero) xscale = one
        if (any(x /= zero .and. x/xscale == zero)) return
        xx = reshape(transpose(x/xscale), [m*ndim])
        ww = one
        if (present(w)) ww = w
        wscale = maxval(ww)
        ww = ww/wscale
        if (any(ww < tiny(one)**0.125_dp)) return
        logtarget = log(s) - 2*log(xscale) - 2*log(wscale)
        limit = real(m, dp)*real(ndim, dp) + one
        if (logtarget >= log(limit)) then
            target = limit
        else
            if (logtarget < log(tiny(one))) return
            target = exp(logtarget)
        end if
        t = zero
        c = zero
        a = zero
        b = zero
        g = zero
        q = zero
        z = zero
        fpint = zero
        nrdata = 0
        n = 2*(k + 1)
        ier = fitpack_ok
        fp = zero
        call fppara(0, ndim, m, uu, m*ndim, xx, ww, zero, one, k, target, nest, 0.001_dp, 20, &
                    k + 1, k + 2, n, t, ncc, c, fp, fpint, z, a, b, g, q, nrdata, ier)
        if (ier > 0) return
        ncoef = n - k - 1
        if (any(g(:ncoef, 1) /= zero)) then
            if (fp_rank_deficient(g, ncoef, nest)) return
        else
            if (fp_rank_deficient(a, ncoef, nest)) return
        end if
        if (.not. all(ieee_is_finite(c))) return
        fp = zero
        l = k + 1
        limit = sqrt(huge(one)/(real(m, dp)*real(ndim, dp)))
        do i = 1, m
            l = fp_knot_interval(t, uu(i), l, ncoef)
            basis = fpbspl(t, n, k, uu(i), l)
            do j = 1, ndim
                offset = (j - 1)*n
                residual = ww(i)*(dot_product(c(offset + l - k:offset + l), basis(:k + 1)) - &
                                  xx((i - 1)*ndim + j))
                if (abs(residual) > limit) return
                fp = fp + residual**2
            end do
        end do
        if (fp > 1.001_dp*target) return
        if (xscale > one) then
            if (maxval(abs(c)) > huge(one)/xscale) return
        end if
        allocate (knots(n), coeff(ncoef, ndim), stat=info)
        if (info /= 0) then
            if (allocated(knots)) deallocate (knots)
            if (allocated(coeff)) deallocate (coeff)
            return
        end if
        knots = bounds(1)+t(:n)*span
        knots(:k + 1) = bounds(1)
        knots(n - k:) = bounds(2)
        do j = 1, ndim
            offset = (j - 1)*n
            coeff(:, j) = c(offset + 1:offset + ncoef)*xscale
        end do
        if (.not. all(ieee_is_finite(knots))) then
            deallocate (knots, coeff)
            return
        end if
        if (minval(knots(k+2:n-k)-knots(k+1:n-k-1))<tiny(one)) then
            deallocate(knots,coeff)
            return
        end if
        ! Verify physical knots and the actual rounded/rescaled coefficients.
        fp=zero
        l=k+1
        do i=1,m
            l=fp_knot_interval(knots,u(i),l,ncoef)
            basis=fpbspl(knots,n,k,u(i),l)
            do j=1,ndim
                residual=ww(i)*(dot_product(coeff(l-k:l,j)/xscale,basis(:k+1))- &
                    xx((i-1)*ndim+j))
                if (abs(residual)>limit) then
                    deallocate(knots,coeff)
                    return
                end if
                fp=fp+residual**2
            end do
        end do
        if (fp>1.001_dp*target) then
            deallocate(knots,coeff)
            return
        end if
        status = 0
    end subroutine fitpack_smooth_parametric

    pure subroutine fppara(iopt, idim, m, u, mx, x, w, ub, ue, k, s, nest, tol, maxit, &
                           k1, k2, n, t, nc, c, fp, fpint, z, a, b, g, q, nrdata, ier)
        implicit none
        !  ..
        !  ..scalar arguments..
        real(dp), intent(in) :: ub !! Left normalized parameter bound, zero.
        real(dp), intent(in) :: ue !! Right normalized parameter bound, one.
        real(dp), intent(in) :: s !! Positive normalized total weighted squared-residual target.
        real(dp), intent(in) :: tol !! Relative residual-target tolerance, 0.001.
        real(dp), intent(inout) :: fp !! Total normalized weighted squared residual.
        integer, intent(in) :: idim !! Number of coordinate channels.
        integer, intent(in) :: maxit !! Maximum smoothing iterations, 20.
        integer, intent(in) :: iopt !! FITPACK task selector, zero for this wrapper.
        integer, intent(in) :: m !! Number of strictly increasing parameter sites.
        integer, intent(in) :: mx !! Length of coordinate-fastest flattened observations.
        integer, intent(in) :: k !! Spline polynomial degree, from one through five.
        integer, intent(in) :: k1 !! Degree plus one.
        integer, intent(in) :: k2 !! Degree plus two.
        integer, intent(in) :: nest !! Allocated knot capacity.
        integer, intent(inout) :: n !! Active full knot count, updated by refinement.
        integer, intent(inout) :: nc !! Allocated flattened coefficient capacity.
        integer, intent(inout) :: ier !! FITPACK success codes <=0 or positive numerical failure.
        !  ..array arguments..
        real(dp), intent(in) :: u(m) !! Normalized strictly increasing parameter sites.
        real(dp), intent(in) :: x(mx) !! Normalized coordinate-fastest flattened observations.
        real(dp), intent(in) :: w(m) !! Positive normalized residual multipliers.
        real(dp), intent(inout) :: t(nest) !! Full knot vector with repeated boundary knots.
        real(dp), intent(inout) :: &
            c(nc) !! Channel-major coefficients with stride n, updated by the solver.
        real(dp), intent(inout) :: &
            fpint(nest) !! Interval residuals and saved knot-pass diagnostics.
        real(dp), intent(inout) :: &
            z(nc) !! Channel-major triangularized observation right-hand sides, stride n.
        real(dp), intent(inout) :: &
            a(nest, k1) !! Upper triangular observation factor in band storage.
        real(dp), intent(inout) :: b(nest, k2) !! Derivative-jump penalty rows.
        real(dp), intent(inout) :: &
          g(nest, k2) !! Penalized upper triangular band factor workspace.
        real(dp), intent(inout) :: q(m, k1) !! Nonzero basis values per observation.
        integer, intent(inout) :: nrdata(nest) !! Interior observation counts per knot interval.

        !  ..local scalars..
        real(dp) :: acc,fac,fpart,fpms,fpold,fp0,f1,f2,f3,p,pinv,p1,p2,p3,rn,store,term,ui,wi
  integer :: i, it, iter, j, jj, j1, j2, k3, l, l0, mk1, nk1, nmax, nmin, nplus, npl1, nrint, n8
        logical :: new, check1, check3, success
        !  ..local arrays..
        real(dp) :: h(max_order + 1), xi(idim)

        fpold = zero
        fp0 = zero
        nplus = 0
        fpms = huge(zero)

        ! calculation of acc, the absolute tolerance for the root of f(p)=s.
        acc = tol*s

        ! determine nmax, the number of knots for spline interpolation.
        nmax = m + k1

        ! *****
        !  part 1: determination of the number of knots and their position
        ! *****
        ! given a set of knots we compute the least-squares curve sinf(u), and the corresponding
        ! sum
        !  of squared residuals fp=f(p=inf).
        !  if iopt=-1 sinf(u) is the requested curve.
        !  if iopt=0 or iopt=1 we check whether we can accept the knots:
        !    if fp <=s we will continue with the current set of knots.
        !    if fp > s we will increase the number of knots and compute the corresponding least-
        !    squares curve until finally fp<=s.
        !  the initial choice of knots depends on the value of s and iopt.
        ! if s=0 we have spline interpolation; in that case the number of knots equals nmax =
        ! m+k+1.
        !    if (s>0 and iopt=0) we first compute the least-squares polynomial curve of degree k;
        !      n = nmin = 2*k+2
        ! iopt=1 we start with the set of knots found at the last call of the routine, except for
        !    the case that s > fp0; then we compute directly the polynomial curve of degree k.
        ! *****

        !  determine nmin, the number of knots for polynomial approximation.
        nmin = 2*k1
        bootstrap: if (iopt >= 0) then

            interpolating: if (s <= zero) then

                !  if s=0, s(u) is an interpolating curve.
                !  check that the required storage space exceeds the available one.
                n = nmax
                if (nmax > nest) then
                    ier = fitpack_insufficient_storage
                    return
                end if

                !  find the position of the interior knots in case of interpolation.
                mk1 = m - k1
                if (mk1 /= 0) then
                    k3 = k/2
                    i = k2
                    j = k3 + 1
                    odd_k: if (k3*2 /= k) then
                        t(k2:k2 + mk1 - 1) = u(j + 1:j + mk1)
                        else odd_k
                        t(k2:k2 + mk1 - 1) = half*(u(j + 1:j + mk1) + u(j:j + mk1 - 1))
                    end if odd_k

                end if

                else interpolating

                !  if s>0 our initial choice of knots depends on the value of iopt.
                !  if iopt=0 or iopt=1 and s>=fp0, we start computing the least-squares
                !  polynomial curve which is a spline curve without interior knots.
                !  if iopt=1 and fp0>s we start computing the least squares spline curve
                !  according to the set of knots found at the last call of the routine.
                use_last_call: if (iopt /= 0 .and. n /= nmin) then
                    fp0 = fpint(n)
                    fpold = fpint(n - 1)
                    nplus = nrdata(n)
                end if use_last_call

                if (fp0 <= s .or. iopt == 0 .or. n == nmin) then
                    n = nmin
                    fpold = zero
                    nplus = 0
                    nrdata(1) = m - 2
                end if
            end if interpolating

        end if bootstrap

        !  main loop for the different sets of knots. m is a safe upper bound
        !  for the number of trials.
        iter = 0
        main_loop: do while (iter <= m)

            iter = iter + 1

            if (n == nmin) ier = fitpack_leastsquares_ok

            ! find nrint, tne number of knot intervals.
            nrint = n - nmin + 1
            ! find the position of the additional knots which are needed for
            ! the b-spline representation of s(u).
            nk1 = n - k1
            t(1:k1) = ub
            t(n - k:n) = ue

            ! compute the b-spline coefficients of the least-squares spline curve
            ! sinf(u). the observation matrix a is built up row by row and
            ! reduced to upper triangular form by givens transformations.
            ! at the same time fp=f(p=inf) is computed.
            fp = zero

            ! initialize the b-spline coefficients and the observation matrix a.
            z(1:nc) = zero
            a(1:nk1, 1:k1) = zero
            l = k1
            coefs: do it = 1, m

                ! fetch the current data point u(it),x(it).
                ui = u(it)
                wi = w(it)
                xi = x((it - 1)*idim + 1:it*idim)*wi

                ! search for knot interval t(l) <= ui < t(l+1).
                l = fp_knot_interval(t, ui, l, nk1)

                ! evaluate the (k+1) non-zero b-splines at ui and store them in q.
                h = fpbspl(t, n, k, ui, l)
                q(it, 1:k1) = h(1:k1)
                h(:k1) = wi*h(:k1)

                ! rotate the new row of the observation matrix into triangle.
                call fp_rotate_row_vec(h, k1, a, nest, xi, z, l - k1, n, idim)

                ! add contribution of this row to the sum of squares of residual right hand sides.
                fp = fp + sum(xi**2)

            end do coefs

            if (ier == fitpack_leastsquares_ok) fp0 = fp
            fpint(n - 1:n) = [fpold, fp0]
            nrdata(n) = nplus

            ! backward substitution to obtain the b-spline coefficients.
            j1 = 1
            do j2 = 1, idim
                c(j1:j1 + nk1 - 1) = fpback(a, z(j1), nk1, k1, nest)
                j1 = j1 + n
            end do

            ! test whether the approximation sinf(u) is an acceptable solution.
            if (iopt < 0) return ! was done already

            fpms = fp - s
            if (abs(fpms) < acc) return

            ! if f(p=inf) < s accept the choice of knots.
            if (fpms < zero) exit main_loop

            ! if n = nmax, sinf(u) is an interpolating spline curve.
            if (n == nmax) then
                ier = fitpack_interpolating_ok
                return
            end if

            ! increase the number of knots.
            ! if n=nest we cannot increase the number of knots because of the storage capacity
            ! limitation.
            if (n == nest) then
                ier = fitpack_insufficient_storage
                return
            end if

            ! determine the number of knots nplus we are going to add.
            if (ier == fitpack_ok) then
                npl1 = nplus*2
                rn = nplus
                if (fpold - fp > acc) npl1 = fp_knots_to_add(rn, fpms, fpold - fp)
                nplus = min(nplus*2, max(npl1, nplus/2, 1))
            else
                nplus = 1
                ier = fitpack_ok
            end if

            ! initialize iterate
            fpold = fp

            ! compute the sum of squared residuals for each knot interval
            ! t(j+k) <= u(i) <= t(j+k+1) and store it in fpint(j),j=1,2,...nrint.
            fpart = zero
            i = 1
            l = k2
            new = .false.
            jj = 0
            square_residuals: do it = 1, m
                if (u(it) >= t(l) .and. l <= nk1) then
                    new = .true.
                    l = l + 1
                end if
                term = zero
                l0 = l - k2
                do j2 = 1, idim
                    fac = dot_product(c(l0 + 1:l0 + k1), q(it, 1:k1))
                    jj = jj + 1
                    term = term + (w(it)*(fac - x(jj)))**2
                    l0 = l0 + n
                end do
                fpart = fpart + term
                if (new) then
                    store = term*half
                    fpint(i) = fpart - store
                    i = i + 1
                    fpart = store
                    new = .false.
                end if
            end do square_residuals

            fpint(nrint) = fpart
            add_new_knots: do l = 1, nplus

                ! add a new knot.
                call fpknot(u, m, t, n, fpint, nrdata, nrint, nest, ione)

                ! if n=nmax we locate the knots as for interpolation
                if (n == nmax) then
                    !  find the position of the interior knots in case of interpolation.
                    mk1 = m - k1
                    if (mk1 /= 0) then
                        k3 = k/2
                        i = k2
                        j = k3 + 2
                        do jj = 1, mk1
                            t(i) = merge(u(j), (u(j) + u(j - 1))*half, k3*2 /= k)
                            i = i + 1
                            j = j + 1
                        end do
                    end if

                    ! restart main loop
                    iter = 0
                    cycle main_loop

                end if

                ! test whether we cannot further increase the number of knots.
                if (n == nest) exit add_new_knots
            end do add_new_knots
            !  restart the computations with the new set of knots.
        end do main_loop

        ! test whether the least-squares kth degree polynomial curve is a solution
        ! of our approximation problem.
        if (ier == fitpack_leastsquares_ok) return

        ! *****
        !  part 2: determination of the smoothing spline curve sp(u).
        ! *****
        !  we have determined the number of knots and their position.
        ! we now compute the b-spline coefficients of the smoothing curve sp(u). the observation
        ! matrix a
        ! is extended by the rows of matrix b expressing that the kth derivative discontinuities
        ! of
        ! sp(u)
        ! at the interior knots t(k+2),...t(n-k-1) must be zero. the corresponding weights of
        ! these
        !  additional rows are set to 1/p.
        !  iteratively we then have to determine the value of p such that f(p), the sum of squared
        ! residuals be = s. we already know that the least squares kth degree polynomial curve
        ! corresponds
        ! to p=0, and that the least-squares spline curve corresponds to p=infinity. the iteration
        ! process
        ! which is proposed here, makes use of rational interpolation. since f(p) is a convex and
        ! strictly
        ! decreasing function of p, it can be approximated by a rational function r(p) =
        ! (u*p+v)/(p+w).
        ! three values of p(p1,p2,p3) with corresponding values of f(p)
        ! (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s)
        ! are used to calculate the new value of p such that r(p)=s. convergence is guaranteed by
        ! taking
        !  f1>0 and f3<zero
        ! *****

        ! evaluate the discontinuity jump of the kth derivative of the b-splines at the knots
        ! t(l),l=k+2,...n-k-1 and store in b.
        call fpdisc(t, n, k2, b, nest)

        !  initial value for p.
        p1 = zero
        f1 = fp0 - s
        p3 = -one
        f3 = fpms
        p = sum(a(1:nk1, 1))/nk1
        check1 = fp_false
        check3 = fp_false
        n8 = n - nmin

        !  iteration process to find the root of f(p) = s.
        iter = 0
        find_root: do while (iter < maxit)

            iter = iter + 1

            ! the rows of matrix b with weight 1/p are rotated into the
            ! triangularised observation matrix a which is stored in g.
            pinv = one/p
            c(1:nc) = z(1:nc)
            g(1:nk1, 1:k1) = a(1:nk1, 1:k1)
            g(1:nk1, k2) = zero
            b_rows: do it = 1, n8
                ! the row of matrix b is rotated into triangle by givens transformation
                h(1:k2) = b(it, 1:k2)*pinv
                xi = zero
                call fp_rotate_shifted_vec(h, k2, g, nest, xi, c, it, nk1, n, idim)
            end do b_rows

            ! backward substitution to obtain the b-spline coefficients.
            j1 = 1
            do j2 = 1, idim
                c(j1:j1 + nk1 - 1) = fpback(g, c(j1), nk1, k2, nest)
                j1 = j1 + n
            end do

            ! computation of f(p).
            fp = zero
            l = k2
            jj = 0
            get_fp: do it = 1, m
                if (u(it) >= t(l) .and. l <= nk1) l = l + 1
                l0 = l - k2
                term = zero
                do j2 = 1, idim
                    fac = dot_product(c(l0 + 1:l0 + k1), q(it, 1:k1))
                    jj = jj + 1
                    term = term + (fac - x(jj))**2
                    l0 = l0 + n
                end do
                fp = fp + term*w(it)**2
            end do get_fp

            ! success! the approximation sp(u) is an acceptable solution.
            fpms = fp - s
            if (abs(fpms) < acc) return

            ! find the new value of p and carry out one more step.
        call root_finding_iterate(p1, f1, p2, f2, p3, f3, p, fpms, acc, check1, check3, success)
            if (.not. success) then
                ier = fitpack_s_too_small
                return
            end if

        end do find_root

        ! maximum number of iterations reached
        ier = fitpack_maxit
        return

    end subroutine fppara

    pure subroutine fp_rotate_row_vec(h, band, a, nest, xi, z, j_start, n, idim)
        implicit none
        integer, intent(in) :: band !! Number of entries in the row to rotate.
        integer, intent(in) :: nest !! Allocated knot capacity.
        integer, intent(in) :: n !! Leading dimension of the coordinate right-hand side array.
        integer, intent(in) :: idim !! Number of coordinate channels.
        integer, intent(in) :: j_start !! First factor row, minus one for unshifted rotation.
        real(dp), intent(inout) :: h(band) !! Working row eliminated by Givens rotations.
        real(dp), intent(inout) :: &
            a(nest, *) !! Upper triangular band factor updated by rotations.
        real(dp), intent(inout) :: &
            xi(idim) !! Coordinate right-hand side contributions, updated by rotations.
        real(dp), intent(inout) :: &
            z(n, idim) !! Right-hand sides, one column per coordinate channel.

        real(dp) :: cos, sin, piv
        integer :: i, j

        j = j_start
        do i = 1, band
            j = j + 1
            piv = h(i)
            if (equal(piv, zero)) cycle

            ! calculate parameters of givens transformation
            call fpgivs(piv, a(j, 1), cos, sin)

            ! apply to all idim right-hand sides (elemental fprota)
            call fprota(cos, sin, xi, z(j, 1:idim))

            ! apply to remaining row elements and matrix
            if (i < band) call fprota(cos, sin, h(i + 1:band), a(j, 2:band - i + 1))
        end do
    end subroutine fp_rotate_row_vec

    pure subroutine fp_rotate_shifted_vec(h, band, a, nest, xi, z, j_start, j_end, n, idim)
        implicit none
        integer, intent(in) :: band !! Number of entries in the row to rotate.
        integer, intent(in) :: nest !! Allocated knot capacity.
        integer, intent(in) :: n !! Leading dimension of the coordinate right-hand side array.
        integer, intent(in) :: idim !! Number of coordinate channels.
        integer, intent(in) :: j_start !! First factor row, minus one for unshifted rotation.
        integer, intent(in) :: j_end !! Last active factor row to update.
        real(dp), intent(inout) :: h(band) !! Working row eliminated by Givens rotations.
        real(dp), intent(inout) :: &
            a(nest, *) !! Upper triangular band factor updated by rotations.
        real(dp), intent(inout) :: &
            xi(idim) !! Coordinate right-hand side contributions, updated by rotations.
        real(dp), intent(inout) :: &
            z(n, idim) !! Right-hand sides, one column per coordinate channel.

        real(dp) :: cos, sin, piv
        integer :: i2, j

        do j = j_start, j_end
            piv = h(1)
            call fpgivs(piv, a(j, 1), cos, sin)
            call fprota(cos, sin, xi, z(j, 1:idim))
            if (j < j_end) then
                i2 = min(j_end - j, band - 1) + 1
                call fprota(cos, sin, h(2:i2), a(j, 2:i2))
                h(1:i2) = [h(2:i2), zero]
            end if
        end do
    end subroutine fp_rotate_shifted_vec



      pure subroutine fpperi(iopt,x,y,w,m,k,s,nest,tol,maxit, &
                            k1,k2,n,t,c,fp,fpint,z,a1,a2,b,g1,g2,q,nrdata,ier)
          implicit none

      !  ..scalar arguments..
      real(dp), intent(in) :: s !! Nonnegative normalized residual target.
      real(dp), intent(in) :: tol !! Relative residual tolerance, 0.001.
      real(dp), intent(inout) :: fp !! Achieved weighted squared residual.
      integer, intent(in) :: iopt !! FITPACK task selector; zero for a fresh periodic fit.
      integer, intent(in) :: m !! Number of coordinates, including the period endpoint.
      integer, intent(in) :: k !! Spline degree, one through five.
      integer, intent(in) :: k1 !! Degree plus one.
      integer, intent(in) :: k2 !! Degree plus two.
      integer, intent(in) :: nest !! Allocated knot and coefficient capacity.
      integer, intent(in) :: maxit !! Maximum smoothing iterations, 20.
      integer, intent(inout) :: n !! Active knot or triangular-system size.
      integer, intent(inout) :: ier !! FITPACK success <=0 or failure >0.

      !  ..array arguments..
      real(dp), intent(in) :: &
          x(m) !! Strictly increasing normalized sample coordinates, including the period endpoint.
      real(dp), intent(in) :: y(m) !! Normalized observations; the final entry is ignored.
      real(dp), intent(in) :: &
          w(m) !! Positive normalized residual multipliers; the final entry is ignored.
      real(dp), intent(inout) :: t(nest) !! Full periodic knot vector.
      real(dp), intent(inout) :: c(nest) !! Periodic spline coefficients.
      real(dp), intent(inout) :: &
          fpint(nest) !! Per-interval residual sums and saved refinement diagnostics.
      real(dp), intent(inout) :: z(nest) !! Triangularized right-hand side.
      real(dp), intent(inout) :: a1(nest,k1) !! Banded observation factor.
      real(dp), intent(inout) :: a2(nest,k) !! Wrapped columns of the observation factor.
      real(dp), intent(inout) :: b(nest,k2) !! Derivative-jump rows or wrapped factor.
      real(dp), intent(inout) :: g1(nest,k2) !! Banded penalized factor.
      real(dp), intent(inout) :: g2(nest,k1) !! Wrapped penalized factor.
      real(dp), intent(inout) :: q(m,k1) !! Nonzero basis values at each observation.
      integer, intent(inout) :: nrdata(nest) !! Interior data counts per knot interval.

      !  ..local scalars..
      real(dp) :: acc,cos,c1,d1,fpart,fpms,fpold,fp0,f1,f2,f3,p,per,pinv,p1,p2,p3, &
                  sin,store,term,wi,xi,yi,rn
      integer :: i,ik,it,iter,i1,i2,i3,j,jk,jper,j1,j2,kk,kk1,k3,l,l0,l1,l5,mm,m1,new, &
                 nk1,nk2,nmax,nmin,nplus,npl1,nrint,n10,n11,n7,n8
      !  ..local arrays..
      real(dp) :: h(max_order+1),h1(degree_5+2),h2(degree_5+1)
      logical :: done,check1,check3,success

      fpold = zero
      fp0   = zero
      nplus = 0
      fpms  = huge(zero)

      !  calculation of acc, the absolute tolerance for the root of f(p)=s.
      acc = tol*s

      !  determine nmax, the number of knots for periodic spline interpolation
      nmax = m+2*k

      ! *****
      !  part 1: determination of the number of knots and their position
      ! *****
      !  given a set of knots we compute the least-squares periodic spline sinf(x). if the
      !  sum f(p=inf) <= s we accept the choice of knots. the initial choice of knots depends on the
      !  value of s and iopt.
      ! For s=0, interpolation uses nmax=m+2*k knots.
      !    if s > 0 and iopt=0 iopt=0 we first compute the least-squares polynomial of degree k;
      ! For n=nmin=2*k+2, the only periodic polynomial is a constant.
      !      iopt=1 we start with the set of knots found at the last call of the routine, except for
      !      the case that s > fp0; then we compute directly the least-squares periodic polynomial.
      ! *****

      m1   = m-1
      kk   = k
      kk1  = k1
      k3   = 3*k+1
      nmin = 2*k1
      ier  = fitpack_ok

      ! determine the length of the period of s(x).
      per = x(m)-x(1)
      if (iopt>=0) then

         interp_or_fit: if (s<=zero .and. nmax/=nmin) then

            !  if s=0, s(x) is an interpolating spline.
            n = nmax

            !  test whether the required storage space exceeds the available one.
            if (n>nest) then
               ier = fitpack_insufficient_storage
               return
            end if

            !  find the position of the interior knots in case of interpolation.
            call fpperi_reset_interp(k,m,n,nest,kk,kk1,x,y,t,c,fp,per,fp0,s,fpint,nrdata,done)
            if (done) then
               ier = fitpack_interpolating_ok
               return
            end if

          else interp_or_fit

             !  if s > 0 our initial choice depends on the value of iopt.
             !  if iopt=0 or iopt=1 and s>=fp0, we start computing the least-squares
             !  periodic polynomial. (i.e. a constant function).
             !  if iopt=1 and fp0>s we start computing the least-squares periodic
             !  spline according the set of knots found at the last call of the routine.
             if (iopt/=0 .and. n/=nmin) then
                fp0   = fpint(n)
                fpold = fpint(n-1)
                nplus = nrdata(n)
             endif

             ! the case that s(x) is a constant function is treated separetely.
             ! find the least-squares constant c1 and compute fp0 at the same time.
             if (iopt==0 .or. n==nmin .or. (iopt/=0 .and. n/=nmin .and. s>=fp0)) then
                fp0 = zero
                d1 = zero
                c1 = zero
                do it=1,m1
                   wi = w(it)
                   yi = y(it)*wi
                   call fpgivs(wi,d1,cos,sin)
                   call fprota(cos,sin,yi,c1)
                   fp0 = fp0+yi**2
                end do
                c1 = c1/d1

                ! test whether that constant function is a solution of our problem.
                fpms = fp0-s

                is_constant: if (fpms<acc .or. nmax==nmin) then

                   ier = fitpack_leastsquares_ok
                   ! the least-squares constant function c1 is a solution of our problem. a constant
                   ! function is a spline of degree k with all b-spline coefficients equal to c1.
                   do i=1,k1
                      rn = k1-i
                      t(i) = x(1)-rn*per
                      c(i) = c1
                      j = i+k1
                      rn = i-1
                      t(j) = x(m)+rn*per
                   end do
                   n = nmin
                   fp = fp0
                   fpint(n-1:n) = [zero,fp0]
                   nrdata(n) = 0
                   return

                endif is_constant

                fpold = fp0
                ! test whether the required storage space exceeds the available one.
                if (nmin>=nest) then
                   ier = fitpack_insufficient_storage
                   return
                endif

                ! start computing the least-squares periodic spline with one interior knot.
                nplus = 1
                n = nmin+1
                mm = (m+1)/2
                t(k2) = x(mm)
                nrdata(1) = mm-2
                nrdata(2) = m1-mm
             endif

          endif interp_or_fit

      endif

      ! Refine the knots, with at most m trials before an interpolation reset.
      iter = 0
      update_knots: do while (iter<m)

         iter = iter+1

         ! find nrint, the number of knot intervals.
         nrint = n-nmin+1

         ! Extend the knot vector periodically for the B-spline representation.
         ! of s(x). if we take
         !      t(k+1) = x(1), t(n-k) = x(m)
         !      t(k+1-j) = t(n-k-j) - per, j=1,2,...k
         !      t(n-k+j) = t(k+1+j) + per, j=1,2,...k
         ! then s(x) is a periodic spline with period per if the b-spline coefficients satisfy the
         ! following conditions
         !      c(n7+j) = c(j), j=1,...k   (**)   with n7=n-2*k-1.
         t(k1) = x(1)
         nk1 = n-k1
         nk2 = nk1+1
         t(nk2) = x(m)
         do j=1,k
            i1 = nk2+j
            i2 = nk2-j
            j1 = k1+j
            j2 = k1-j
            t(i1) = t(j1)+per
            t(j2) = t(i2)-per
         end do

         ! compute the b-spline coefficients c(j),j=1,...n7 of the least-squares periodic spline
         ! sinf(x). the observation matrix a is built up row by row while taking into account
         ! condition (**) and is reduced to triangular form by givens transformations .
         ! at the same time fp=f(p=inf) is computed.
         ! the n7 x n7 triangularised upper matrix a has the form
         !            ! a1 '    !
         !        a = !    ' a2 !
         !            ! 0  '    !
         ! a2 is n7 by k; a1 is n10 by n10, upper triangular with bandwidth k+1.

         ! initialization.
         z (1:nk1) = zero
         a1(1:nk1,1:kk1) = zero
         n7   = nk1-k
         n10  = n7-kk
         jper = 0
         fp   = zero
         l    = k1
         get_coefs: do it=1,m1
            ! fetch the current data point x(it),y(it)
            xi = x(it)
            wi = w(it)
            yi = y(it)*wi

            ! search for knot interval t(l) <= xi < t(l+1).
            l = fp_knot_interval(t, xi, l, nk1)

            ! evaluate the (k+1) non-zero b-splines at xi and store them in q.
            h = fpbspl(t,n,k,xi,l)

            q(it,:k1) = h(:k1)
            h(:k1)    = h(:k1)*wi

            l5 = l-k1

            ! test whether the b-splines nj,k+1(x),j=1+n7,...nk1 are all zero at xi
            if (l5>=n10) then

               if (jper==0) then
                  ! initialize the matrix a2.
                  a2(1:n7,1:kk) = zero
                  jk = n10+1
                  do i=1,kk
                     ik = jk
                     do j=1,kk1
                        if (ik<=0) exit
                        a2(ik,i) = a1(ik,j)
                        ik = ik-1
                     end do
                     jk = jk+1
                  end do
                  jper = 1
               endif

               ! if one of the b-splines nj,k+1(x),j=n7+1,...nk1 is not zero at xi
               ! we take account of condition (**) for setting up the new row
               ! of the observation matrix a. this row is stored in the arrays h1
               ! (the part with respect to a1) and h2 (the part with respect to a2).
               h1 = zero
               h2 = zero
               j  = l5-n10
               do i=1,kk1
                  j  = j+1
                  l0 = j
                  l1 = l0-kk
                  do while (l1>max(0,n10))
                    l0 = l1-n10
                    l1 = l0-kk
                  end do
                  if (l1>0) then
                     h1(l1) = h(i)
                  else
                     h2(l0) = h2(l0)+h(i)
                  end if
               end do

               ! rotate the new row of the observation matrix into triangle
               ! by givens transformations.
               call fp_rotate_row_2mat(h1, kk1, h2, kk, a1, a2, nest, yi, z, 1, n10)

          else ! l5<n10

             ! rotation of the new row of the observation matrix into triangle in case the b-splines
             ! nj,k+1(x),j=n7+1,...n-k-1 are all zero at xi.
             ! rotate the new row of the observation matrix into triangle.
             call fp_rotate_row(h, kk1, a1, yi, z, l5)

          endif

          ! add contribution of this row to the sum of squares of residual rhs
          fp = fp+yi**2

        end do get_coefs

        fpint(n-1:n) = [fpold,fp0]
        nrdata(n)    = nplus

        ! backward substitution to obtain the b-spline coefficients c(j),j=1,.n
        if (fp_periodic_rank_deficient(a1,a2,n7,kk)) then
            ier=fitpack_s_too_small
            return
        end if
        c(1:n7) = fpbacp(a1,a2,z,n7,kk,kk1,nest)

        ! calculate from condition (**) the coefficients c(j+n7),j=1,2,...k.
        c(n7+1:n7+k) = c(1:k)

        ! success! if the approximation sinf(x) is an acceptable solution.
        fpms = fp-s
        if (iopt<0 .or. abs(fpms)<acc) return

        ! if f(p=inf) < s accept the choice of knots.
        if (fpms<zero) exit update_knots

        ! if n=nmax, sinf(x) is an interpolating spline.
        if (n==nmax) then
           ier = fitpack_interpolating_ok
           return
        end if

        ! increase the number of knots.
        ! if n=nest we cannot increase the number of knots because of the
        ! storage capacity limitation.
        if (n==nest) then
           ier = fitpack_insufficient_storage
           return
        end if

        ! determine the number of knots nplus we are going to add.
        npl1  = nplus*2
        if (fpold-fp>acc) npl1 = fp_knots_to_add(real(nplus,dp),fpms,fpold-fp)
        nplus = min(nplus*2,max(npl1,nplus/2,1))
        fpold = fp

        ! compute the sum(wi*(yi-s(xi))**2) for each knot interval
        ! t(j+k) <= xi <= t(j+k+1) and store it in fpint(j),j=1,2,...nrint.
        fpart = zero
        i     = 1
        l     = k1
        do it=1,m1
           if (x(it)>=t(l)) then
              new = 1
              l = l+1
           endif
           term = zero
           l0 = l-k2
           do j=1,k1
              l0 = l0+1
              term = term+c(l0)*q(it,j)
           end do
           term = (w(it)*(term-y(it)))**2
           fpart = fpart+term
           if (new/=0) then
              if (l>k2) then
                 store    = term*half
                 fpint(i) = fpart-store
                 i        = i+1
                 fpart    = store
              else
                 fpint(nrint) = term
              endif
              new = 0
           endif
        end do
        fpint(nrint) = fpint(nrint)+fpart

        add_new_knots: do l=1,nplus

           ! add a new knot
           call fpknot(x,m,t,n,fpint,nrdata,nrint,nest,ione)

           ! if n=nmax we locate the knots as for interpolation.
           if (n==nmax) then

               ! find the position of the interior knots in case of interpolation.
               call fpperi_reset_interp(k,m,n,nest,kk,kk1,x,y,t,c,fp,per,fp0,s,fpint,nrdata,done)
               if (done) then
                  ier = fitpack_interpolating_ok
                  return
               else
                  ! restart iteration
                  iter = 0
                  cycle update_knots
               end if

           end if
           ! test whether we cannot further increase the number of knots.
           if (n==nest) exit add_new_knots

        end do add_new_knots
      !  restart the computations with the new set of knots.
      end do update_knots

      ! *****
      !  part 2: determination of the smoothing periodic spline sp(x).
      ! *****
      !  we have determined the number of knots and their position.
      ! Extend a by rows of b penalizing kth-derivative jumps at interior knots.
      ! The actual rotated rows below have multiplier 1/p. Find p such that
      ! f(p)=sum((w*(y-sp(x)))**2)=s. The constant fit corresponds to p=0;
      ! the unpenalized periodic spline corresponds to p=infinity.
      ! Rational interpolation uses three residuals f1,f2,f3 and the function
      ! r(p)=(u*p+v)/(p+w), with a bracket satisfying f1>0 and f3<0.
      ! *****

      !  evaluate the discontinuity jump of the kth derivative of the b-splines at the knots
      !  t(l),l=k+2,...n-k-1 and store in b.
      call fpdisc(t,n,k2,b,nest)

      !  initial value for p.
      p   = zero
      p1  = zero
      f1  = fp0-s
      p3  = -one
      f3  = fpms
      n11 = n10-1
      n8  = n7-1
      l   = n7
      do i=1,k
         j = k+1-i
         p = p+a2(l,j)
         l = l-1
         if (l==0) exit
      end do
      if (l/=0) p = p+sum(a1(1:n10,1))
      rn = n7
      p = rn/p
      check1 = fp_false
      check3 = fp_false

      ! iteration process to find the root of f(p) = s.
      find_root: do iter=1,maxit

         !  form the matrix g  as the matrix a extended by the rows of matrix b.
         !  the rows of matrix b with weight 1/p are rotated into
         !  the triangularised observation matrix a.
         !  after triangularisation our n7 x n7 matrix g takes the form
         !            ! g1 '    !
         !        g = !    ' g2 !
         !            ! 0  '    !
         !  with g2 a n7 x (k+1) matrix and g1 a n11 x n11 upper triangular
         !  matrix of bandwidth k+2. ( n11 = n7-k-1)

         pinv = one/p

         !  store matrix a into g
         c (1:n7)      = z(1:n7)
         g1(1:n7,1:k1) = a1(1:n7,1:k1)
         g1(1:n7,k2)   = zero
         g2(1:n7,1)    = zero
         g2(1:n7,2:k1) = a2(1:n7,1:k)

         l = n10
         do j=1,k1
            if (l<=0) exit
            g2(l,1) = a1(l,j)
            l = l-1
         end do

         n8_rows: do it=1,n8

            ! fetch a new row of matrix b and store it in the arrays h1 (the part
            ! with respect to g1) and h2 (the part with respect to g2).
            yi = zero
            h1(:k2) = zero
            h2(:k1) = zero
            if (it<=n11) then

                l  = it
                l0 = it
                storeh2: do j=1,k2
                   if (l0==n10) then
                      l0 = 1
                      do l1=j,k2
                         h2(l0) = b(it,l1)*pinv
                         l0 = l0+1
                      end do
                      exit storeh2
                   endif
                   h1(j) = b(it,j)*pinv
                   l0 = l0+1
                end do storeh2

            else
               l = 1
               i = it-n10
               do j=1,k2
                  i = i+1
                  l0 = i
                  l1 = l0-k1
                  do while (l1>max(0,n11))
                     l0 = l1-n11
                     l1 = l0-k1
                  end do
                  if (l1>0) then
                     h1(l1) = b(it,j)*pinv
                  else
                     h2(l0) = h2(l0)+b(it,j)*pinv
                  end if
               end do
            endif

            ! rotate this row into triangle by givens transformations without square roots.
            call fp_rotate_row_2mat(h1, k2, h2, k1, g1, g2, nest, yi, c, l, n11)
         end do n8_rows

         ! backward substitution to obtain the b-spline coefficients c(j),j=1,2,...n7 of sp(x).
         if (fp_periodic_rank_deficient(g1,g2,n7,k1)) then
             ier=fitpack_s_too_small
             return
         end if
         c(:n7) = fpbacp(g1,g2,c,n7,k1,k2,nest)

         ! calculate from condition (**) the b-spline coefficients c(n7+j),j=1,.
         c(n7+1:n7+k) = c(1:k)

         ! computation of f(p).
         fp = zero
         l = k1
         do it=1,m1
            if (x(it)>=t(l)) l = l+1
            l0 = l-k2
            term = dot_product(c(l0+1:l0+k1),q(it,1:k1))
            fp = fp+(w(it)*(term-y(it)))**2
         end do

         ! test whether the approximation sp(x) is an acceptable solution.
         fpms = fp-s
         if (abs(fpms)<acc) return

         ! test whether the maximal number of iterations is reached.
         if (iter==maxit) then
            ier = fitpack_maxit
            return
         end if

         ! find the new value of p and carry out one more step.
         call root_finding_iterate(p1,f1,p2,f2,p3,f3,p,fpms,acc,check1,check3,success)
         if (.not.success) then
            ier = fitpack_s_too_small
            return
         end if

      end do find_root

      return

      end subroutine fpperi

      pure subroutine fpperi_reset_interp(k,m,n,nest,kk,kk1,x,y,t,c,fp,per,fp0,s,fpint,nrdata,done)
          implicit none
         integer, intent(in) :: k !! Spline degree, one through five.
         integer, intent(in) :: m !! Number of coordinates, including the period endpoint.
         integer, intent(in) :: n !! Active knot or triangular-system size.
         integer, intent(in) :: nest !! Allocated knot and coefficient capacity.
         integer, intent(inout) :: kk !! Effective band degree.
         integer, intent(inout) :: kk1 !! Effective band order.
         real(dp), intent(in) :: &
          x(m) !! Strictly increasing normalized sample coordinates, including the period endpoint.
         real(dp), intent(in) :: y(m) !! Normalized observations; the final entry is ignored.
         real(dp), intent(in) :: per !! Normalized period length.
         real(dp), intent(in) :: fp0 !! Weighted residual of the constant fit.
         real(dp), intent(in) :: s !! Nonnegative normalized residual target.
         real(dp), intent(inout) :: t(nest) !! Full periodic knot vector.
         real(dp), intent(inout) :: c(nest) !! Periodic spline coefficients.
         real(dp), intent(inout) :: fp !! Achieved weighted squared residual.
         real(dp), intent(inout) :: &
          fpint(nest) !! Per-interval residual sums and saved refinement diagnostics.
         integer, intent(inout) :: nrdata(nest) !! Interior data counts per knot interval.
         logical, intent(out) :: done !! True when linear interpolation is complete.

         integer :: m1

         m1 = m-1
         done = .false.

         k_is_odd: if (mod(k,2)/=0) then

            t(k+2:k+m1) = x(2:m1)

            if (s<=zero) then
               kk  = k-1
               kk1 = k
               if (kk<=0) then
                   t(1:2)       = [t(m)-per,x(1)]
                   t(m+1:m+2)   = [x(m),t(3)+per]
                   c(:m)        = [y(1:m1),y(1)]
                   fp           = zero
                   fpint(n-1:n) = [zero,fp0]
                   nrdata(n)    = 0
                   done = .true.
               endif
            endif

         else k_is_odd

            t(k+2:k+m1) = half*(x(2:m1)+x(1:m1-1))

         endif k_is_odd

      end subroutine fpperi_reset_interp

      pure function fpbacp(a,b,z,n,k,k1,nest) result(c)
          implicit none

      !  ..scalar arguments..
      integer, intent(in) :: n !! Number of independent periodic coefficients.
      integer, intent(in) :: k !! Width of wrapped columns, adjusted for the current factor.
      integer, intent(in) :: k1 !! Banded factor width, one more than k.
      integer, intent(in) :: nest !! Allocated knot and coefficient capacity.
      !  ..array arguments..
      real(dp), intent(in) :: a(nest,k1) !! Banded part of the triangular factor.
      real(dp), intent(in) :: b(nest,k) !! Wrapped columns of the triangular factor.
      real(dp), intent(in) :: z(n) !! Triangularized right-hand side.
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

      pure subroutine fp_rotate_row_2mat(h1, band1, h2, band2, a1, a2, nest, &
                                        yi, z, j1_start, j1_end)
          implicit none
          integer, intent(in) :: band1 !! Width of the banded row.
          integer, intent(in) :: band2 !! Width of the wrapped row.
          integer, intent(in) :: nest !! Allocated knot and coefficient capacity.
          integer, intent(in) :: j1_start !! First banded row to eliminate.
          integer, intent(in) :: j1_end !! Last banded row to eliminate.
          real(dp),    intent(inout) :: h1(band1) !! Working banded row for elimination.
          real(dp),    intent(inout) :: h2(band2) !! Working wrapped row for elimination.
          real(dp),    intent(inout) :: a1(nest,*) !! Banded observation factor.
          real(dp),    intent(inout) :: a2(nest,*) !! Wrapped columns of the observation factor.
          real(dp),    intent(inout) :: yi !! Working right-hand-side value.
          real(dp),    intent(inout) :: z(*) !! Triangularized right-hand side.

          real(dp) :: cos, sin, piv
          integer :: i2, ij, j

          ! phase 1: rotate h1 through a1 with shifting, cross-rotating h2/a2
          do j = j1_start, j1_end
              piv = h1(1)
              if (equal(piv, zero)) then
                  ! zero pivot: shift h1 left without rotation
                  h1(1:band1-1) = h1(2:band1)
                  h1(band1) = zero
                  cycle
              end if
              call fpgivs(piv, a1(j,1), cos, sin)
              call fprota(cos, sin, yi, z(j))
              call fprota(cos, sin, h2(1:band2), a2(j, 1:band2))
              if (j < j1_end) then
                  i2 = min(j1_end - j, band1 - 1) + 1
                  call fprota(cos, sin, h1(2:i2), a1(j, 2:i2))
                  h1(1:i2) = [h1(2:i2), zero]
              endif
          end do

          ! phase 2: rotate h2 through a2 (diagonal at a2(ij,j))
          do j = 1, band2
              ij = j1_end + j
              if (ij <= 0) cycle
              piv = h2(j)
              if (equal(piv, zero)) cycle
              call fpgivs(piv, a2(ij,j), cos, sin)
              call fprota(cos, sin, yi, z(ij))
              if (j < band2) then
                  call fprota(cos, sin, h2(j+1:band2), a2(ij, j+1:band2))
              endif
          end do
      end subroutine fp_rotate_row_2mat
    pure logical function fp_periodic_rank_deficient(a,b,n,k) result(singular)
        implicit none
        real(dp), intent(in) :: a(:,:) !! Banded part of the periodic triangular factor.
        real(dp), intent(in) :: b(:,:) !! Wrapped columns of the periodic triangular factor.
        integer, intent(in) :: n !! Number of independent coefficients.
        integer, intent(in) :: k !! Width of the wrapped factor.
        integer :: i,j,l
        singular=.true.
        l=n
        do i=1,k
            j=k+1-i
            if (.not.ieee_is_finite(b(l,j))) return
            if (equal(b(l,j),zero)) return
            l=l-1
            if (l==0) exit
        end do
        do i=1,l
            if (.not.ieee_is_finite(a(i,1))) return
            if (equal(a(i,1),zero)) return
        end do
        singular=.false.
    end function fp_periodic_rank_deficient
end module scifort_interpolate_curfit_core
