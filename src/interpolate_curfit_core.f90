! SPDX-License-Identifier: MIT AND BSD-3-Clause
! Copyright (c) 2022, Federico Perini. All rights reserved.
! Copyright (c) 2026 SciFort contributors
!
! Adaptive scalar smoothing adapted from perazz/fitpack src/fitpack_core.F90,
! commit 1907c7bad9f5c4557970e906bdd4909eb52b2d4c, based on Paul Dierckx's FITPACK.
! License: LICENSES/perazz-fitpack-BSD-3-Clause.txt. See CODE_PROVENANCE.md.
module scifort_interpolate_curfit_core
    use scifort_kinds, only: dp
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: fitpack_smooth_curve
    integer, parameter :: max_order = 19, degree_5 = 5, ione = 1
    integer, parameter :: fitpack_ok = 0, fitpack_interpolating_ok = -1
    integer, parameter :: fitpack_leastsquares_ok = -2, fitpack_insufficient_storage = 1
    integer, parameter :: fitpack_s_too_small = 2, fitpack_maxit = 3
    real(dp), parameter :: zero = 0.0_dp, one = 1.0_dp, half = 0.5_dp
    logical, parameter :: fp_true = .true., fp_false = .false.
contains

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

end module scifort_interpolate_curfit_core
