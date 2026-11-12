! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Optional adapter to Federico Perini's BSD-3-Clause FITPACK modernization.
! The actual FITPACK source is supplied as a pinned FPM dependency; it is not
! copied into the MIT-licensed SciFort repository.
module scifort_optional_fitpack
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_interpolate_fitpack_compat, only : fitpack_tck_type, fitpack_parametric_tck_type
    use scifort_interpolate_fitpack_bivariate, only : fitpack_bivariate_tck_type
    use fitpack, only : fitpack_curve, fitpack_periodic_curve, &
        fitpack_parametric_curve, fitpack_closed_curve, fitpack_surface, &
        FP_REAL, FP_FLAG, FITPACK_SUCCESS
    implicit none
    private
    integer, parameter, public :: fitpack_smooth_success = 0
    integer, parameter, public :: fitpack_smooth_invalid_input = 1
    integer, parameter, public :: fitpack_smooth_upstream_error = 2
    public :: univariate_spline, periodic_univariate_spline, parametric_smoothing_spline
    public :: bivariate_smoothing_spline, lsq_bivariate_spline

contains
    ! SciPy's UnivariateSpline(x,y,w=...,k=...,s=...) equivalent interface.
    ! The returned TCK can be used with SciFort's splev, splint, splder, etc.
    function univariate_spline(x,y,w,k,s,status,residual) result(tck)
        real(dp), intent(in) :: x(:) !! strictly increasing finite sample coordinates
        real(dp), intent(in) :: y(:) !! sampled scalar ordinates
        real(dp), intent(in), optional :: w(:) !! positive FITPACK residual multipliers
        integer, intent(in), optional :: k !! degree 1..5 (default 3)
        real(dp), intent(in), optional :: s !! nonnegative residual target (default n)
        integer, intent(out), optional :: status !! 0 success, 1 invalid, 2 fit failure
        real(dp), intent(out), optional :: residual !! sum((w*(y-f(x)))**2)
        type(fitpack_tck_type) :: tck
        type(fitpack_curve) :: curve
        integer(FP_FLAG) :: ierr
        integer :: deg, code, n
        real(dp) :: target
        n=size(x)
        deg=3
        if (present(k)) deg=k
        target=real(n,dp)
        if (present(s)) target=s
        code=fitpack_smooth_invalid_input
        if (present(residual)) residual=0.0_dp
        if (deg < 1 .or. deg > 5 .or. n < deg+1 .or. size(y) /= n) goto 900
        if (n < 2) goto 900
        if (.not. all(ieee_is_finite(x)) .or. .not. all(ieee_is_finite(y))) goto 900
        if (any(x(2:) <= x(:n-1))) goto 900
        if (.not. ieee_is_finite(target) .or. target < 0.0_dp) goto 900
        if (present(w)) then
            if (size(w) /= n) goto 900
            if (.not. all(ieee_is_finite(w)) .or. any(w <= 0.0_dp)) goto 900
            ierr=curve%new_fit(real(x,FP_REAL),real(y,FP_REAL), &
                w=real(w,FP_REAL),smoothing=real(target,FP_REAL),order=deg)
        else
            ierr=curve%new_fit(real(x,FP_REAL),real(y,FP_REAL), &
                smoothing=real(target,FP_REAL),order=deg)
        end if
        code=fitpack_smooth_upstream_error
        if (.not. FITPACK_SUCCESS(ierr)) then
            call curve%destroy()
            goto 900
        end if
        if (curve%knots < 2*deg+2) then
            call curve%destroy()
            goto 900
        end if
        tck%t=real(curve%t(:curve%knots),dp)
        tck%c=real(curve%c(:curve%knots),dp)
        tck%k=deg
        if (present(residual)) residual=real(curve%mse(),dp)
        if (.not. tck%valid()) then
            if (allocated(tck%t)) deallocate(tck%t)
            if (allocated(tck%c)) deallocate(tck%c)
            tck%k=-1
        else
            code=fitpack_smooth_success
        end if
        call curve%destroy()
900     if (present(status)) status=code
    end function univariate_spline
    ! SciPy splrep(...,per=True,s=s) / periodic UnivariateSpline-style fit.
    ! A periodic FITPACK TCK retains the original support domain; SciPy's
    ! splev evaluates the polynomial continuation rather than wrapping inputs.
    function periodic_univariate_spline(x,y,w,k,s,status,residual) result(tck)
        real(dp), intent(in) :: x(:) !! strictly increasing observation abscissae
        real(dp), intent(in) :: y(:) !! ordinates, last value is ignored by FITPACK
        real(dp), intent(in), optional :: w(:) !! positive residual weights
        integer, intent(in), optional :: k !! degree, 1..5, default 3
        real(dp), intent(in), optional :: s !! nonnegative weighted residual target
        integer, intent(out), optional :: status !! success, invalid input or upstream error
        real(dp), intent(out), optional :: residual !! FITPACK weighted residual
        type(fitpack_tck_type) :: tck
        type(fitpack_periodic_curve) :: curve
        integer(FP_FLAG) :: ierr
        integer :: deg, code, n
        real(dp) :: target

        n = size(x)
        deg = 3
        if (present(k)) deg = k
        target = real(n,dp)
        if (present(s)) target = s
        code = fitpack_smooth_invalid_input
        if (present(residual)) residual = 0.0_dp
        if (deg < 1 .or. deg > 5 .or. n < deg+2 .or. size(y) /= n) goto 900
        if (.not. all(ieee_is_finite(x)) .or. .not. all(ieee_is_finite(y))) goto 900
        if (any(x(2:) <= x(:n-1))) goto 900
        if (.not. ieee_is_finite(target) .or. target < 0.0_dp) goto 900
        if (present(w)) then
            if (size(w) /= n) goto 900
            if (.not. all(ieee_is_finite(w)) .or. any(w <= 0.0_dp)) goto 900
        end if
        ! Absent optional dummy arguments are forwarded as absent (Fortran 2008).
        ierr = curve%new_fit(x,y,w=w,smoothing=target,order=deg)
        code = fitpack_smooth_upstream_error
        if (.not. FITPACK_SUCCESS(ierr)) goto 800
        if (curve%knots < 2*deg+2) goto 800
        tck%t = real(curve%t(:curve%knots),dp)
        tck%c = real(curve%c(:curve%knots),dp)
        tck%k = deg
        if (present(residual)) residual = real(curve%mse(),dp)
        if (tck%valid()) code = fitpack_smooth_success
800     call curve%destroy()
        if (code /= fitpack_smooth_success) then
            if (allocated(tck%t)) deallocate(tck%t)
            if (allocated(tck%c)) deallocate(tck%c)
            tck%k = -1
        end if
900     if (present(status)) status = code
    end function periodic_univariate_spline

    ! FITPACK parcur / clocur adaptive-knot parametric smoothing, exposed
    ! as the same (knots,coefficients,degree,parameters) TCK as SciFort splprep.
    function parametric_smoothing_spline(x,u,w,k,s,per,status,residual) result(tck)
        real(dp), intent(in) :: x(:,:) !! observations x(coordinates), (npoints, ndim)
        real(dp), intent(in), optional :: u(:) !! strictly increasing parameter values
        real(dp), intent(in), optional :: w(:) !! positive residual weights
        integer, intent(in), optional :: k !! degree 1..5, default 3
        real(dp), intent(in), optional :: s !! nonnegative smoothing residual target
        logical, intent(in), optional :: per !! use FITPACK's closed parametric fit
        integer, intent(out), optional :: status !! success, invalid or upstream failure
        real(dp), intent(out), optional :: residual !! sum of weighted squared residuals
        type(fitpack_parametric_tck_type) :: tck
        type(fitpack_parametric_curve) :: curve
        type(fitpack_closed_curve) :: closed
        integer(FP_FLAG) :: ierr
        integer :: code, deg, n, ndim
        real(dp) :: target
        logical :: is_closed
        real(dp), allocatable :: closed_points(:,:)

        n = size(x,1)
        ndim = size(x,2)
        deg = 3
        if (present(k)) deg = k
        target = real(n,dp)
        if (present(s)) target = s
        is_closed = .false.
        if (present(per)) is_closed = per
        code = fitpack_smooth_invalid_input
        if (present(residual)) residual = 0.0_dp
        if (deg < 1 .or. deg > 5 .or. n < deg+1 .or. ndim < 1 .or. ndim > 10) goto 900
        if (is_closed .and. n < deg+2) goto 900
        if (.not. all(ieee_is_finite(x))) goto 900
        if (.not. ieee_is_finite(target) .or. target < 0.0_dp) goto 900
        if (present(u)) then
            if (size(u) /= n) goto 900
            if (.not. all(ieee_is_finite(u))) goto 900
            if (any(u(2:) <= u(:n-1))) goto 900
        else
            ! Avoid zero chord lengths in the implicit chord-length parameterization.
            if (any(sum(abs(x(2:,:) - x(:n-1,:)),dim=2) <= 0.0_dp)) goto 900
        end if
        if (present(w)) then
            if (size(w) /= n) goto 900
            if (.not. all(ieee_is_finite(w)) .or. any(w <= 0.0_dp)) goto 900
        end if
        code = fitpack_smooth_upstream_error
        if (is_closed) then
            ! Match scipy.splprep(per=True): the final ordinate is ignored.
            ! FITPACK clocur requires bitwise equality at the endpoints.
            closed_points = transpose(x)
            closed_points(:,n) = closed_points(:,1)
            ierr = closed%new_fit(closed_points,u=u,w=w,smoothing=target,order=deg)
            if (FITPACK_SUCCESS(ierr)) then
                call copy_parametric_tck(closed,tck)
                if (present(residual)) residual = real(closed%mse(),dp)
            end if
            call closed%destroy()
        else
            ierr = curve%new_fit(transpose(x),u=u,w=w,smoothing=target,order=deg)
            if (FITPACK_SUCCESS(ierr)) then
                call copy_parametric_tck(curve,tck)
                if (present(residual)) residual = real(curve%mse(),dp)
            end if
            call curve%destroy()
        end if
        if (tck%valid()) code = fitpack_smooth_success
        if (code /= fitpack_smooth_success) then
            if (allocated(tck%t)) deallocate(tck%t)
            if (allocated(tck%c)) deallocate(tck%c)
            if (allocated(tck%u)) deallocate(tck%u)
            tck%k = -1
        end if
900     if (present(status)) status = code
    end function parametric_smoothing_spline

    ! Adaptive-knot scattered bivariate smoothing, using FITPACK surfit.
    ! This corresponds to SciPy's SmoothBivariateSpline (with equal degrees)
    ! or bisplrep(task=0); exact coefficient equality is not guaranteed.
    function bivariate_smoothing_spline(x,y,z,w,k,s,status,residual) result(tck)
        real(dp), intent(in) :: x(:) !! scattered x sites
        real(dp), intent(in) :: y(:) !! scattered y sites
        real(dp), intent(in) :: z(:) !! observations
        real(dp), intent(in), optional :: w(:) !! positive residual weights
        integer, intent(in), optional :: k !! common x/y spline degree, 1:5
        real(dp), intent(in), optional :: s !! nonnegative weighted residual target
        integer, intent(out), optional :: status !! 0 success, 1 invalid, 2 fitting failure
        real(dp), intent(out), optional :: residual !! squared weighted residual
        type(fitpack_bivariate_tck_type) :: tck
        type(fitpack_surface) :: surface
        integer(FP_FLAG) :: ierr
        integer :: n, deg, code
        real(dp) :: target

        n = size(x)
        deg = 3
        if (present(k)) deg = k
        target = real(n,dp)
        if (present(s)) target = s
        code = fitpack_smooth_invalid_input
        if (present(residual)) residual = 0.0_dp
        if (deg < 1 .or. deg > 5 .or. n < (deg+1)**2) goto 900
        if (size(y) /= n .or. size(z) /= n) goto 900
        if (.not.all(ieee_is_finite(x)) .or. .not.all(ieee_is_finite(y))) goto 900
        if (.not.all(ieee_is_finite(z))) goto 900
        if (minval(x) >= maxval(x) .or. minval(y) >= maxval(y)) goto 900
        if (.not.ieee_is_finite(target) .or. target < 0.0_dp) goto 900
        if (present(w)) then
            if (size(w) /= n) goto 900
            if (.not.all(ieee_is_finite(w)) .or. any(w <= 0.0_dp)) goto 900
        end if
        code = fitpack_smooth_upstream_error
        call surface%new_points(x,y,z,w=w)
        call prepare_bivariate_workspace(surface,[deg,deg])
        ierr = surface%fit(smoothing=target,order=deg)
        if (FITPACK_SUCCESS(ierr)) then
            call copy_bivariate_tck(surface,tck)
            if (present(residual)) residual = real(surface%mse(),dp)
            if (tck%valid()) code = fitpack_smooth_success
        end if
        call surface%destroy()
900     if (code /= fitpack_smooth_success) call clear_bivariate_tck(tck)
        if (present(status)) status = code
    end function bivariate_smoothing_spline

    ! Fixed-knot scattered bivariate least squares (SciPy LSQBivariateSpline).
    ! Interior knots are supplied; FITPACK's 'surfit(iopt=-1)' solves the fit.
    function lsq_bivariate_spline(x,y,z,tx,ty,w,kx,ky,status,residual) result(tck)
        real(dp), intent(in) :: x(:) !! scattered x coordinates
        real(dp), intent(in) :: y(:) !! scattered y coordinates
        real(dp), intent(in) :: z(:) !! observation ordinates
        real(dp), intent(in) :: tx(:) !! strictly increasing interior x knots
        real(dp), intent(in) :: ty(:) !! strictly increasing interior y knots
        real(dp), intent(in), optional :: w(:) !! positive residual weights
        integer, intent(in), optional :: kx !! x spline degree, default 3
        integer, intent(in), optional :: ky !! y spline degree, default 3
        integer, intent(out), optional :: status !! status
        real(dp), intent(out), optional :: residual !! weighted squared residual
        type(fitpack_bivariate_tck_type) :: tck
        type(fitpack_surface) :: surface
        integer(FP_FLAG) :: ierr
        integer :: nx,ny,n,degx,degy,code
        real(dp) :: xmin,xmax,ymin,ymax

        n = size(x)
        degx = 3
        degy = 3
        if (present(kx)) degx = kx
        if (present(ky)) degy = ky
        code = fitpack_smooth_invalid_input
        if (present(residual)) residual = 0.0_dp
        if (degx < 1 .or. degx > 5 .or. degy < 1 .or. degy > 5) goto 900
        if (size(y) /= n .or. size(z) /= n .or. n < (degx+1)*(degy+1)) goto 900
        if (.not.all(ieee_is_finite(x)) .or. .not.all(ieee_is_finite(y))) goto 900
        if (.not.all(ieee_is_finite(z))) goto 900
        if (.not.all(ieee_is_finite(tx)) .or. .not.all(ieee_is_finite(ty))) goto 900
        if (size(tx) > 1) then
            if (any(tx(2:) <= tx(:size(tx)-1))) goto 900
        end if
        if (size(ty) > 1) then
            if (any(ty(2:) <= ty(:size(ty)-1))) goto 900
        end if
        xmin = minval(x)
        xmax = maxval(x)
        ymin = minval(y)
        ymax = maxval(y)
        if (xmin >= xmax .or. ymin >= ymax) goto 900
        if (any(tx <= xmin) .or. any(tx >= xmax)) goto 900
        if (any(ty <= ymin) .or. any(ty >= ymax)) goto 900
        if (present(w)) then
            if (size(w) /= n) goto 900
            if (.not.all(ieee_is_finite(w)) .or. any(w <= 0.0_dp)) goto 900
        end if
        nx = size(tx)+2*(degx+1)
        ny = size(ty)+2*(degy+1)
        code = fitpack_smooth_upstream_error
        call surface%new_points(x,y,z,w=w)
        if (size(surface%t,1) < max(nx,ny)) goto 800
        call prepare_bivariate_workspace(surface,[degx,degy])
        surface%knots = [nx,ny]
        surface%t(:degx+1,1) = xmin
        surface%t(degx+2:degx+1+size(tx),1) = tx
        surface%t(nx-degx:nx,1) = xmax
        surface%t(:degy+1,2) = ymin
        surface%t(degy+2:degy+1+size(ty),2) = ty
        surface%t(ny-degy:ny,2) = ymax
        ! This calls surfit with its fixed-knots branch, iopt=-1.
        ! An explicit smoothing argument forces exactly one FITPACK call.
        ! Calling surface%least_squares() without it makes get_smoothing()
        ! run three passes and silently replaces fixed knots after pass one.
        surface%iopt = -1
        ierr = surface%fit(smoothing=0.0_dp,keep_knots=.true.)
        if (FITPACK_SUCCESS(ierr)) then
            call copy_bivariate_tck(surface,tck)
            if (present(residual)) residual = real(surface%mse(),dp)
            if (tck%valid()) code = fitpack_smooth_success
        end if
800     call surface%destroy()
900     if (code /= fitpack_smooth_success) call clear_bivariate_tck(tck)
        if (present(status)) status = code
    end function lsq_bivariate_spline

    ! The upstream surface constructor allocates workspace for its default
    ! bicubic order. When using other degrees, resize the workspaces to the
    ! documented surfit formulas before fitting (without vendoring FITPACK).
    subroutine prepare_bivariate_workspace(surface,degree)
        type(fitpack_surface), intent(inout) :: surface !! FITPACK surface object
        integer, intent(in) :: degree(2) !! per-axis spline degrees
        integer :: uv(2), nx(2), b1, b2, bx, by, km, nmax, n
        surface%order = degree
        nx = surface%nest
        n = surface%m
        nmax = surface%nmax
        uv = nx-degree-1
        km = maxval(degree)+1
        bx = degree(1)*uv(2)+degree(2)+1
        by = degree(2)*uv(1)+degree(1)+1
        b1 = min(bx,by)
        if (bx <= by) then
            b2 = b1+uv(2)-degree(2)
        else
            b2 = b1+uv(1)-degree(1)
        end if
        ! surfit's implementation uses ib3=(kx+1)*(nyest-ky-1)+1
        ! in its workspace bound, even for asymmetric degrees. Allocate at
        ! least that amount; the upstream convenience constructor otherwise
        ! underallocates for kx > ky and returns FITPACK_INPUT_ERROR.
        b2 = max(b2,(degree(1)+1)*uv(2)+1)
        surface%lwrk = product(uv)*(2+b1+b2) + &
            2*(sum(uv)+km*(n+nmax)+nmax-sum(degree))+b2+1
        surface%liwrk = n+product(nx-2*degree-1)
        surface%lwrk2 = product(uv)*(b2+1)+b2
        deallocate(surface%c,surface%wrk,surface%iwrk,surface%wrk2)
        allocate(surface%c(product(uv)),source=0.0_dp)
        allocate(surface%wrk(surface%lwrk),source=0.0_dp)
        allocate(surface%iwrk(surface%liwrk),source=0)
        allocate(surface%wrk2(surface%lwrk2),source=0.0_dp)
    end subroutine prepare_bivariate_workspace

    subroutine copy_bivariate_tck(surface,tck)
        type(fitpack_surface), intent(in) :: surface !! fitted upstream surface
        type(fitpack_bivariate_tck_type), intent(out) :: tck !! SciPy/FITPACK TCK
        integer :: nx,ny
        nx = surface%knots(1)-surface%order(1)-1
        ny = surface%knots(2)-surface%order(2)-1
        if (nx < 1 .or. ny < 1) return
        tck%tx = real(surface%t(:surface%knots(1),1),dp)
        tck%ty = real(surface%t(:surface%knots(2),2),dp)
        tck%c = real(surface%c(:nx*ny),dp)
        tck%kx = surface%order(1)
        tck%ky = surface%order(2)
    end subroutine copy_bivariate_tck

    subroutine clear_bivariate_tck(tck)
        type(fitpack_bivariate_tck_type), intent(inout) :: tck !! reset on error
        if (allocated(tck%tx)) deallocate(tck%tx)
        if (allocated(tck%ty)) deallocate(tck%ty)
        if (allocated(tck%c)) deallocate(tck%c)
        tck%kx = -1
        tck%ky = -1
    end subroutine clear_bivariate_tck

    subroutine copy_parametric_tck(curve,tck)
        class(fitpack_parametric_curve), intent(in) :: curve !! fitted upstream spline
        type(fitpack_parametric_tck_type), intent(out) :: tck !! SciFort tuple
        integer :: j,nc,nk
        nk = curve%knots
        nc = nk-curve%order-1
        if (nc < curve%order+1) return
        if (.not. allocated(curve%t) .or. .not. allocated(curve%c) &
            .or. .not. allocated(curve%u)) return
        if (size(curve%t) < nk .or. size(curve%c) < nk*curve%idim) return
        tck%t = real(curve%t(:nk),dp)
        tck%u = real(curve%u,dp)
        allocate(tck%c(nc,curve%idim))
        do j=1,curve%idim
            tck%c(:,j) = real(curve%c((j-1)*nk+1:(j-1)*nk+nc),dp)
        end do
        tck%k = curve%order
    end subroutine copy_parametric_tck
end module scifort_optional_fitpack
