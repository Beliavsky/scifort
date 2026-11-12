! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Dependency-free evaluation of tensor-product bivariate FITPACK B-splines.
! The coefficient ordering is the one used by scipy.interpolate.bisplrep:
! c((ix-1)*ny+iy), with the y coefficient varying fastest.
module scifort_interpolate_fitpack_bivariate
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_interpolate_bspline, only : bspline_interpolator_type, &
        bspline_array_interpolator_type, bspline, make_interp_spline, &
        bspline_status_success
    implicit none
    private

    integer, parameter, public :: bispline_success = 0
    integer, parameter, public :: bispline_invalid_input = 1
    integer, parameter, public :: bispline_singular = 2

    type, public :: fitpack_bivariate_tck_type
        real(dp), allocatable :: tx(:) !! x-direction knot vector
        real(dp), allocatable :: ty(:) !! y-direction knot vector
        real(dp), allocatable :: c(:)  !! y-fastest flattened tensor coefficients
        integer :: kx = -1 !! x-degree, in 1:5
        integer :: ky = -1 !! y-degree, in 1:5
    contains
        procedure :: valid => bivariate_valid
        procedure :: integral => bivariate_integral
        procedure :: partial_derivative => bivariate_partial_derivative
        procedure :: insert_knot => bivariate_insert_knot
    end type fitpack_bivariate_tck_type

    public :: bisplev, bisplev_points, bispline_integral, rect_bivariate_spline
    public :: bispline_partial_derivative, bispline_insert_knot
    interface bisplev
        module procedure bisplev_scalar
        module procedure bisplev_grid
    end interface bisplev

contains
    ! RectBivariateSpline(s=0) equivalent for a rectangular, complete real grid.
    ! Reuse the established one-dimensional collocation solver successively
    ! along the x and y axes. Fortran array ordering has x varying fastest in
    ! the observation grid and y varying fastest in the FITPACK coefficients.
    function rect_bivariate_spline(x, y, z, kx, ky, status) result(tck)
        real(dp), intent(in) :: x(:) !! strictly increasing x sites, at least kx+1
        real(dp), intent(in) :: y(:) !! strictly increasing y sites, at least ky+1
        real(dp), intent(in) :: z(:,:) !! real observations with shape (size(x),size(y))
        integer, intent(in), optional :: kx !! x spline degree in 1:5, default 3
        integer, intent(in), optional :: ky !! y spline degree in 1:5, default 3
        integer, intent(out), optional :: status !! 0 success, 1 invalid, 2 singular fit
        type(fitpack_bivariate_tck_type) :: tck
        type(bspline_array_interpolator_type) :: fit_x, fit_y
        real(dp), allocatable :: cx(:,:), cy(:,:)
        integer :: degx, degy, code

        code = bispline_invalid_input
        degx = 3
        degy = 3
        if (present(kx)) degx = kx
        if (present(ky)) degy = ky
        if (degx < 1 .or. degx > 5 .or. degy < 1 .or. degy > 5) goto 900
        if (size(x) < degx+1 .or. size(y) < degy+1) goto 900
        if (size(z,1) /= size(x) .or. size(z,2) /= size(y)) goto 900
        if (.not. all(ieee_is_finite(x)) .or. .not. all(ieee_is_finite(y))) goto 900
        if (.not. all(ieee_is_finite(z))) goto 900
        if (any(x(2:) <= x(:size(x)-1))) goto 900
        if (any(y(2:) <= y(:size(y)-1))) goto 900

        code = bispline_singular
        fit_x = make_interp_spline(x,z,k=degx,status=code)
        if (code /= bspline_status_success) then
            code = bispline_singular
            goto 900
        end if
        cx = fit_x%coefficients()
        fit_y = make_interp_spline(y,transpose(cx),k=degy,status=code)
        if (code /= bspline_status_success) then
            code = bispline_singular
            goto 900
        end if
        cy = fit_y%coefficients()
        tck%tx = fit_x%knots()
        tck%ty = fit_y%knots()
        tck%c = reshape(cy,[size(cy)])
        tck%kx = degx
        tck%ky = degy
        code = bispline_success
900     if (present(status)) status = code
    end function rect_bivariate_spline

    ! Construct the tensor-product spline for a mixed partial derivative.
    ! This is the transformation used by BivariateSpline.partial_derivative:
    ! differentiation lowers the polynomial degree, removes one knot from
    ! each end of that axis, and differences adjacent coefficient fibers.
    ! SciPy requires dx < kx and dy < ky; the remaining degrees are >= 1.
    function bivariate_partial_derivative(self, dx, dy, status) result(derivative)
        class(fitpack_bivariate_tck_type), intent(in) :: self !! original spline tuple
        integer, intent(in) :: dx !! number of x derivatives, 0 <= dx < kx
        integer, intent(in) :: dy !! number of y derivatives, 0 <= dy < ky
        integer, intent(out), optional :: status !! 0 success; 1 invalid input
        type(fitpack_bivariate_tck_type) :: derivative
        type(fitpack_bivariate_tck_type) :: work
        real(dp), allocatable :: cnew(:)
        real(dp) :: denominator
        integer :: nx, ny, ix, iy, r, code

        code = bispline_invalid_input
        if (.not. self%valid()) goto 900
        if (dx < 0 .or. dy < 0 .or. dx >= self%kx .or. dy >= self%ky) goto 900
        work = self
        nx = size(work%tx)-work%kx-1
        ny = size(work%ty)-work%ky-1
        ! FITPACK may return a padded coefficient vector. Only nx*ny entries
        ! participate in the actual spline; the derivative has no padding.
        if (size(work%c) > nx*ny) work%c = work%c(:nx*ny)
        do r = 1, dx
            allocate(cnew((nx-1)*ny))
            do ix = 1, nx-1
                denominator = work%tx(ix+work%kx+1)-work%tx(ix+1)
                if (denominator <= 0.0_dp) goto 900
                do iy = 1, ny
                    cnew((ix-1)*ny+iy) = real(work%kx,dp)* &
                        (work%c(ix*ny+iy)-work%c((ix-1)*ny+iy))/denominator
                end do
            end do
            call move_alloc(cnew,work%c)
            work%tx = work%tx(2:size(work%tx)-1)
            work%kx = work%kx-1
            nx = nx-1
        end do
        do r = 1, dy
            allocate(cnew(nx*(ny-1)))
            do iy = 1, ny-1
                denominator = work%ty(iy+work%ky+1)-work%ty(iy+1)
                if (denominator <= 0.0_dp) goto 900
                do ix = 1, nx
                    cnew((ix-1)*(ny-1)+iy) = real(work%ky,dp)* &
                        (work%c((ix-1)*ny+iy+1)-work%c((ix-1)*ny+iy))/denominator
                end do
            end do
            call move_alloc(cnew,work%c)
            work%ty = work%ty(2:size(work%ty)-1)
            work%ky = work%ky-1
            ny = ny-1
        end do
        if (.not. work%valid()) goto 900
        derivative = work
        code = bispline_success
900     if (present(status)) status = code
    end function bivariate_partial_derivative

    function bispline_partial_derivative(tck, dx, dy, status) result(derivative)
        type(fitpack_bivariate_tck_type), intent(in) :: tck !! input tensor-product spline
        integer, intent(in) :: dx !! x derivative count
        integer, intent(in) :: dy !! y derivative count
        integer, intent(out), optional :: status !! 0 success or 1 invalid input
        type(fitpack_bivariate_tck_type) :: derivative
        derivative = tck%partial_derivative(dx,dy,status)
    end function bispline_partial_derivative

    ! Shape-preserving tensor-product knot refinement via Boehm insertion.
    ! Each coefficient fiber is independently refined with the same 1D
    ! transformation, leaving spline values and all supported derivatives
    ! unchanged. There is no refitting or evaluation-grid approximation.
    function bivariate_insert_knot(self, axis, x, m, status) result(refined)
        class(fitpack_bivariate_tck_type), intent(in) :: self !! input surface
        integer, intent(in) :: axis !! 1 inserts x knots; 2 inserts y knots
        real(dp), intent(in) :: x !! new knot in the active domain of axis
        integer, intent(in), optional :: m !! repetitions, default one; zero gives a copy
        integer, intent(out), optional :: status !! 0 success, 1 invalid input
        type(fitpack_bivariate_tck_type) :: refined
        type(bspline_interpolator_type) :: fiber, inserted
        real(dp), allocatable :: new_coeffs(:), values(:)
        integer :: reps, nx, ny, ix, iy, code, nnew

        code = bispline_invalid_input
        if (.not. self%valid()) goto 900
        if (axis /= 1 .and. axis /= 2) goto 900
        if (.not. ieee_is_finite(x)) goto 900
        reps = 1
        if (present(m)) reps = m
        if (reps < 0) goto 900
        nx = size(self%tx)-self%kx-1
        ny = size(self%ty)-self%ky-1
        if (axis == 1) then
            if (x < self%tx(self%kx+1) .or. x > self%tx(nx+1)) goto 900
            if (reps > self%kx+1-count(self%tx == x)) goto 900
            ! Extract every x fiber from FITPACK's y-fastest coefficient
            ! vector. Reinsert as a 2D matrix with the same ordering.
            nnew = nx+reps
            allocate(new_coeffs(nnew*ny))
            do iy = 1,ny
                fiber = bspline(self%tx,self%c(iy:nx*ny:ny),self%kx)
                inserted = fiber%insert_knot(x,m=reps,status=code)
                if (code /= bspline_status_success) goto 900
                values = inserted%coefficients()
                new_coeffs(iy:nnew*ny:ny) = values
                if (iy == 1) refined%tx = inserted%knots()
            end do
            refined%ty = self%ty
        else
            if (x < self%ty(self%ky+1) .or. x > self%ty(ny+1)) goto 900
            if (reps > self%ky+1-count(self%ty == x)) goto 900
            nnew = ny+reps
            allocate(new_coeffs(nx*nnew))
            do ix = 1,nx
                fiber = bspline(self%ty,self%c((ix-1)*ny+1:ix*ny),self%ky)
                inserted = fiber%insert_knot(x,m=reps,status=code)
                if (code /= bspline_status_success) goto 900
                values = inserted%coefficients()
                new_coeffs((ix-1)*nnew+1:ix*nnew) = values
                if (ix == 1) refined%ty = inserted%knots()
            end do
            refined%tx = self%tx
        end if
        refined%c = new_coeffs
        refined%kx = self%kx
        refined%ky = self%ky
        if (.not. refined%valid()) goto 800
        code = bispline_success
        goto 900
800     code = bispline_invalid_input
900     if (code /= bispline_success) then
            if (allocated(refined%tx)) deallocate(refined%tx)
            if (allocated(refined%ty)) deallocate(refined%ty)
            if (allocated(refined%c)) deallocate(refined%c)
            refined%kx = -1
            refined%ky = -1
        end if
        if (present(status)) status = code
    end function bivariate_insert_knot

    function bispline_insert_knot(tck, axis, x, m, status) result(refined)
        type(fitpack_bivariate_tck_type), intent(in) :: tck !! original bivariate TCK
        integer, intent(in) :: axis !! axis 1=x or 2=y
        real(dp), intent(in) :: x !! knot location in selected axis
        integer, intent(in), optional :: m !! number of insertions, default 1
        integer, intent(out), optional :: status !! 0 success; 1 invalid input
        type(fitpack_bivariate_tck_type) :: refined
        refined = tck%insert_knot(axis,x,m,status)
    end function bispline_insert_knot

    ! The integral over a rectangle is the tensor product of two exact
    ! univariate B-spline integrations. Like FITPACK's dblint, out-of-domain
    ! parts of the requested rectangle are ignored, rather than extrapolated.
    function bivariate_integral(self, xa, xb, ya, yb, status) result(v)
        class(fitpack_bivariate_tck_type), intent(in) :: self !! fitted tensor B-spline
        real(dp), intent(in) :: xa !! first x integration bound
        real(dp), intent(in) :: xb !! second x integration bound
        real(dp), intent(in) :: ya !! first y integration bound
        real(dp), intent(in) :: yb !! second y integration bound
        integer, intent(out), optional :: status !! 0 success, 1 invalid
        real(dp) :: v
        real(dp), allocatable :: cx(:)
        type(bspline_interpolator_type) :: curve
        real(dp) :: xlo, xhi, ylo, yhi, sign_xy
        integer :: ix, nx, ny, code

        code = bispline_invalid_input
        v = ieee_value(0.0_dp,ieee_quiet_nan)
        if (.not. self%valid()) goto 900
        if (.not. ieee_is_finite(xa) .or. .not. ieee_is_finite(xb) .or. &
            .not. ieee_is_finite(ya) .or. .not. ieee_is_finite(yb)) goto 900
        nx = size(self%tx)-self%kx-1
        ny = size(self%ty)-self%ky-1
        xlo = max(min(xa,xb),self%tx(self%kx+1))
        xhi = min(max(xa,xb),self%tx(nx+1))
        ylo = max(min(ya,yb),self%ty(self%ky+1))
        yhi = min(max(ya,yb),self%ty(ny+1))
        if (xhi <= xlo .or. yhi <= ylo) then
            v = 0.0_dp
            code = bispline_success
            goto 900
        end if
        sign_xy = 1.0_dp
        if (xb < xa) sign_xy = -sign_xy
        if (yb < ya) sign_xy = -sign_xy
        allocate(cx(nx))
        do ix = 1,nx
            curve = bspline(self%ty,self%c((ix-1)*ny+1:ix*ny),self%ky)
            cx(ix) = curve%integrate(ylo,yhi)
        end do
        curve = bspline(self%tx,cx,self%kx)
        v = sign_xy*curve%integrate(xlo,xhi)
        code = bispline_success
900     if (present(status)) status = code
    end function bivariate_integral

    function bispline_integral(xa, xb, ya, yb, tck, status) result(v)
        real(dp), intent(in) :: xa !! first x integration bound
        real(dp), intent(in) :: xb !! second x integration bound
        real(dp), intent(in) :: ya !! first y integration bound
        real(dp), intent(in) :: yb !! second y integration bound
        type(fitpack_bivariate_tck_type), intent(in) :: tck !! bivariate spline tuple
        integer, intent(out), optional :: status !! 0 success, 1 invalid
        real(dp) :: v
        v = tck%integral(xa,xb,ya,yb,status)
    end function bispline_integral

    pure logical function bivariate_valid(self) result(ok)
        class(fitpack_bivariate_tck_type), intent(in) :: self !! TCK to validate
        integer :: nx, ny
        ok = allocated(self%tx) .and. allocated(self%ty) .and. allocated(self%c)
        if (.not. ok) return
        if (self%kx < 1 .or. self%kx > 5 .or. self%ky < 1 .or. self%ky > 5) then
            ok = .false.
            return
        end if
        nx = size(self%tx)-self%kx-1
        ny = size(self%ty)-self%ky-1
        ok = nx >= self%kx+1 .and. ny >= self%ky+1
        if (.not. ok) return
        ok = size(self%c) >= nx*ny
        if (.not. ok) return
        ok = all(ieee_is_finite(self%tx)) .and. all(ieee_is_finite(self%ty)) .and. &
             all(ieee_is_finite(self%c(:nx*ny)))
        if (.not. ok) return
        ok = all(self%tx(2:) >= self%tx(:size(self%tx)-1)) .and. &
             all(self%ty(2:) >= self%ty(:size(self%ty)-1))
        if (.not. ok) return
        ok = self%tx(self%kx+1) < self%tx(nx+1) .and. &
             self%ty(self%ky+1) < self%ty(ny+1)
    end function bivariate_valid

    pure integer function locate_span(t, k, x) result(span)
        real(dp), intent(in) :: t(:) !! nondecreasing knots
        integer, intent(in) :: k !! polynomial degree
        real(dp), intent(in) :: x !! coordinate, already clipped to domain
        integer :: left, right, mid, nc
        nc = size(t)-k-1
        if (x >= t(nc+1)) then
            span = nc
            return
        end if
        left = k+1
        right = nc+1
        do while (right-left > 1)
            mid = (right+left)/2
            if (x < t(mid)) then
                right = mid
            else
                left = mid
            end if
        end do
        span = left
    end function locate_span

    ! At most k+1 basis functions contribute. Reuse the tested 1D de Boor
    ! implementation, including its analytical derivatives; no global grid
    ! or dense tensor basis matrix is created.
    subroutine local_basis(t, k, x, nu, span, vals)
        real(dp), intent(in) :: t(:) !! knots
        integer, intent(in) :: k !! degree
        real(dp), intent(in) :: x !! query coordinate
        integer, intent(in) :: nu !! derivative order
        integer, intent(out) :: span !! active knot interval
        real(dp), intent(out) :: vals(0:) !! local B-spline basis/derivatives
        real(dp), allocatable :: coefficients(:)
        type(bspline_interpolator_type) :: one
        real(dp) :: xx
        integer :: i, n, j
        n = size(t)-k-1
        xx = max(t(k+1),min(x,t(n+1)))  ! FITPACK bispev boundary convention
        span = locate_span(t,k,xx)
        vals = 0.0_dp
        if (nu > k) return
        allocate(coefficients(n))
        do j=0,k
            i = span-k+j
            coefficients = 0.0_dp
            coefficients(i) = 1.0_dp
            one = bspline(t,coefficients,k)
            if (nu == 0) then
                vals(j) = one%evaluate(xx)
            else
                vals(j) = one%derivative(xx,nu)
            end if
        end do
    end subroutine local_basis

    function bisplev_scalar(x,y,tck,dx,dy,status) result(z)
        real(dp), intent(in) :: x !! evaluation x-coordinate
        real(dp), intent(in) :: y !! evaluation y-coordinate
        type(fitpack_bivariate_tck_type), intent(in) :: tck !! fitted tensor B-spline
        integer, intent(in), optional :: dx !! x partial derivative order, default 0
        integer, intent(in), optional :: dy !! y partial derivative order, default 0
        integer, intent(out), optional :: status !! 0 success, 1 invalid
        real(dp) :: z
        real(dp) :: bx(0:5), by(0:5)
        integer :: nu_x, nu_y, sx, sy, ix, iy, ny, code
        z = ieee_value(0.0_dp,ieee_quiet_nan)
        code = bispline_invalid_input
        if (.not.tck%valid()) goto 900
        if (.not.ieee_is_finite(x) .or. .not.ieee_is_finite(y)) goto 900
        nu_x=0
        nu_y=0
        if (present(dx)) nu_x=dx
        if (present(dy)) nu_y=dy
        if (nu_x < 0 .or. nu_y < 0 .or. nu_x >= tck%kx .or. nu_y >= tck%ky) goto 900
        call local_basis(tck%tx,tck%kx,x,nu_x,sx,bx(:tck%kx))
        call local_basis(tck%ty,tck%ky,y,nu_y,sy,by(:tck%ky))
        ny = size(tck%ty)-tck%ky-1
        z = 0.0_dp
        do ix=0,tck%kx
            do iy=0,tck%ky
                z = z + bx(ix)*by(iy)*tck%c((sx-tck%kx+ix-1)*ny+sy-tck%ky+iy)
            end do
        end do
        code = bispline_success
900     if (present(status)) status = code
    end function bisplev_scalar

    function bisplev_grid(x,y,tck,dx,dy,status) result(z)
        real(dp), intent(in) :: x(:) !! x evaluation coordinates
        real(dp), intent(in) :: y(:) !! y evaluation coordinates
        type(fitpack_bivariate_tck_type), intent(in) :: tck !! bivariate TCK
        integer, intent(in), optional :: dx !! x derivative order
        integer, intent(in), optional :: dy !! y derivative order
        integer, intent(out), optional :: status !! evaluation status
        real(dp), allocatable :: z(:,:)
        integer :: i, j, code, first_error
        allocate(z(size(x),size(y)))
        first_error = bispline_success
        do j=1,size(y)
            do i=1,size(x)
                z(i,j) = bisplev_scalar(x(i),y(j),tck,dx,dy,code)
                if (code /= bispline_success) first_error = code
            end do
        end do
        if (present(status)) status = first_error
    end function bisplev_grid

    function bisplev_points(x,y,tck,dx,dy,status) result(z)
        real(dp), intent(in) :: x(:) !! paired x coordinates
        real(dp), intent(in) :: y(:) !! paired y coordinates
        type(fitpack_bivariate_tck_type), intent(in) :: tck !! bivariate TCK
        integer, intent(in), optional :: dx !! x derivative order
        integer, intent(in), optional :: dy !! y derivative order
        integer, intent(out), optional :: status !! evaluation status
        real(dp), allocatable :: z(:)
        integer :: i, code, first_error
        if (size(x) /= size(y)) then
            allocate(z(0))
            if (present(status)) status = bispline_invalid_input
            return
        end if
        allocate(z(size(x)))
        first_error = bispline_success
        do i=1,size(x)
            z(i) = bisplev_scalar(x(i),y(i),tck,dx,dy,code)
            if (code /= bispline_success) first_error = code
        end do
        if (present(status)) status = first_error
    end function bisplev_points
end module scifort_interpolate_fitpack_bivariate
