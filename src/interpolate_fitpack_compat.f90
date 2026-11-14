! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Dependency-free FITPACK TCK compatibility layer: spline representation,
! Interpolation, native scalar/parametric smoothing, fixed-knot LSQ and calculus.
module scifort_interpolate_fitpack_compat
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_interpolate_curfit_core, only : fitpack_smooth_curve,fitpack_smooth_parametric
    use scifort_interpolate_curfit_core, only : fitpack_periodic_curve
    use scifort_interpolate_curfit_core, only : fitpack_periodic_parametric
    use scifort_interpolate_bspline, only : bspline_interpolator_type, bspline, &
        interpolated_univariate_spline, lsq_univariate_spline, &
        bspline_status_success
    implicit none
    private

    integer, parameter, public :: fitpack_tck_success = 0
    integer, parameter, public :: fitpack_tck_invalid_input = 1
    integer, parameter, public :: fitpack_tck_unsupported_smoothing = 2
    integer, parameter, public :: fitpack_tck_singular = 3
    integer, parameter, public :: fitpack_tck_out_of_bounds = 4

    type, public :: fitpack_tck_type
        ! FITPACK/SciPy-style (t,c,k), where splrep pads c to the length of t.
        real(dp), allocatable :: t(:)
        real(dp), allocatable :: c(:)
        integer :: k = -1
    contains
        procedure :: valid => tck_valid
        procedure :: ncoeff => tck_ncoeff
    end type

    ! A parametric FITPACK tuple keeps the independent parameter values and a
    ! shared knot vector for the coordinate channels.  Unlike splrep, SciPy's
    ! splprep stores *unpadded* coefficients for each coordinate.
    type, public :: fitpack_parametric_tck_type
        real(dp), allocatable :: t(:)
        real(dp), allocatable :: c(:,:)  ! (number_of_B_spline_coefficients, coordinates)
        real(dp), allocatable :: u(:)    ! input parameter, estimated if absent
        integer :: k = -1
    contains
        procedure :: valid => parametric_tck_valid
        procedure :: channels => parametric_tck_channels
    end type fitpack_parametric_tck_type

    public :: splrep, splprep, splev, splint, splder, splantider, spalde
    interface splev
        module procedure splev_scalar, splev_array
        module procedure splev_parametric_scalar, splev_parametric_array
    end interface

contains

    pure logical function parametric_tck_valid(self) result(ok)
        class(fitpack_parametric_tck_type), intent(in) :: self !! parametric spline
        integer :: nc
        ok = allocated(self%t) .and. allocated(self%c) .and. allocated(self%u)
        if (.not. ok) return
        if (self%k < 1 .or. self%k > 5) then
            ok = .false.
            return
        end if
        nc = size(self%t)-self%k-1
        ok = nc >= self%k+1 .and. size(self%c,1) == nc .and. &
             size(self%c,2) >= 1 .and. size(self%u) >= self%k+1
        if (.not. ok) return
        ok = all(ieee_is_finite(self%t)) .and. all(ieee_is_finite(self%c)) .and. &
             all(ieee_is_finite(self%u))
        if (.not. ok) return
        ok = all(self%t(2:) >= self%t(:size(self%t)-1)) .and. &
             all(self%u(2:) > self%u(:size(self%u)-1))
        if (.not. ok) return
        ok = self%t(self%k+1) < self%t(nc+1)
    end function parametric_tck_valid

    pure integer function parametric_tck_channels(self) result(n)
        class(fitpack_parametric_tck_type), intent(in) :: self !! parametric spline
        n = 0
        if (self%valid()) n = size(self%c,2)
    end function parametric_tck_channels

    pure elemental function chord_difference(a,b,scale) result(delta)
        implicit none
        real(dp), intent(in) :: a !! Finite current coordinate component.
        real(dp), intent(in) :: b !! Finite preceding coordinate component.
        real(dp), intent(in) :: scale !! Common positive finite coordinate scale.
        real(dp) :: delta
        ! Subtract first for nearby large values; divide first only to avoid overflow.
        if (a>=0.0_dp .and. b<0.0_dp) then
            if (a>huge(a)+b) then
                delta=a/scale-b/scale
                return
            end if
        else if (a<0.0_dp .and. b>=0.0_dp) then
            if (b>huge(a)+a) then
                delta=a/scale-b/scale
                return
            end if
        end if
        delta=(a-b)/scale
    end function chord_difference

    ! SciPy splprep: interpolation (s=0), automatic-knot smoothing (s>0),
    ! or supplied full knot vector (nonperiodic task=-1).
    ! Input ordinates are transposed relative to scipy.splprep's sequence:
    ! coordinates(sample,dimension), with up to arbitrary coordinate count.
    function splprep(x,u,w,k,s,task,t,ub,ue,per,status) result(tck)
        implicit none
        real(dp), intent(in) :: x(:,:) !! curve coordinates (samples, dimensions)
        real(dp), intent(in), optional :: u(:) !! strictly increasing parameter
        real(dp), intent(in), optional :: w(:) !! positive FITPACK residual multipliers
        integer, intent(in), optional :: k !! degree 1:5 (default 3)
        real(dp), intent(in), optional :: s !! total weighted squared-residual target, default zero
        integer, intent(in), optional :: task !! 0=interpolation/smoothing, -1=fixed knots
        real(dp), intent(in), optional :: t(:) !! *complete* FITPACK knots if task=-1
        real(dp), intent(in), optional :: ub !! left parameter domain bound
        real(dp), intent(in), optional :: ue !! right parameter domain bound
        logical, intent(in), optional :: per !! Closed periodic curve, default false.
        integer, intent(out), optional :: status !! 0=success, 1=invalid, 2=unsupported, 3=singular
        type(fitpack_parametric_tck_type) :: tck
        type(bspline_interpolator_type) :: fit
        real(dp), allocatable :: param(:),distances(:),coefficients(:),knots(:),interior(:)
        real(dp), allocatable :: delta(:)
        real(dp) :: left, right, total, target, coordinate_scale
        integer :: n, ndim, deg, mode, code, fitcode, j, nc, used
        logical :: periodic,scaled_chords
        code = fitpack_tck_invalid_input
        n = size(x,1)
        ndim = size(x,2)
        periodic=.false.
        if (present(per)) periodic=per
        deg = 3
        if (present(k)) deg = k
        mode = 0
        if (present(task)) mode = task
        target = 0.0_dp
        if (present(s)) target = s
        if (deg < 1 .or. deg > 5 .or. n <= deg .or. ndim < 1) goto 900
        if (mode /= 0 .and. mode /= -1) goto 900
        used=n
        if (periodic) used=n-1
        if (.not. all(ieee_is_finite(x(:used,:)))) goto 900
        if (.not. ieee_is_finite(target)) goto 900
        if (target < 0.0_dp) goto 900
        if (periodic .and. mode==-1) then
            if (present(t)) code=fitpack_tck_unsupported_smoothing
            goto 900
        end if
        if (present(w)) then
            if (size(w) /= n) goto 900
            if (.not. all(ieee_is_finite(w(:used)))) goto 900
            if (any(w(:used) <= 0.0_dp)) goto 900
        end if
        allocate(param(n),stat=fitcode)
        if (fitcode /= 0) then
            code=fitpack_tck_singular
            goto 900
        end if
        if (present(u)) then
            if (size(u) /= n) goto 900
            if (.not. all(ieee_is_finite(u))) goto 900
            param = u
            if (any(param(2:) <= param(:n-1))) goto 900
        else
            allocate(distances(n),delta(ndim),stat=fitcode)
            if (fitcode /= 0) then
                code=fitpack_tck_singular
                goto 900
            end if
            coordinate_scale=1.0_dp
            scaled_chords=periodic .or. (mode==0 .and. target>0.0_dp)
            if (scaled_chords) then
                coordinate_scale=maxval(abs(x(:used,:)))
                if (coordinate_scale == 0.0_dp) goto 900
            end if
            distances(1) = 0.0_dp
            do j=2,n
                ! Common scaling preserves Euclidean chord lengths up to a constant
                ! while preventing overflow before subtracting large coordinates.
                if (periodic .and. j==n) then
                    delta=chord_difference(x(1,:),x(n-1,:),coordinate_scale)
                else if (scaled_chords) then
                    delta=chord_difference(x(j,:),x(j-1,:),coordinate_scale)
                else
                    delta=x(j,:)-x(j-1,:)
                end if
                total = maxval(abs(delta))
                if (total == 0.0_dp) then
                    if (periodic .and. j==n) then
                        if (any(x(1,:)/=x(n-1,:))) code=fitpack_tck_singular
                    else
                        if (any(x(j,:)/=x(j-1,:))) code=fitpack_tck_singular
                    end if
                    goto 900
                end if
                if (total <= 0.0_dp .or. .not. ieee_is_finite(total)) goto 900
                distances(j)=distances(j-1)+total*sqrt(sum((delta/total)**2))
            end do
            total = distances(n)
            if (total <= 0.0_dp .or. .not. ieee_is_finite(total)) goto 900
            param = distances/total
        end if
        left = param(1)
        right = param(n)
        if (present(ub)) left = ub
        if (present(ue)) right = ue
        if (.not. ieee_is_finite(left) .or. .not. ieee_is_finite(right)) goto 900
        if (left > param(1) .or. right < param(n) .or. left >= right) goto 900
        if (periodic) then
            if (present(t)) goto 900
            if (left/=param(1) .or. right/=param(n)) goto 900
            call fitpack_periodic_parametric(param,x,w,deg,target,tck%t,tck%c,code)
            if (code==0) then
                call move_alloc(param,tck%u)
                tck%k=deg
                if (.not.tck%valid()) code=fitpack_tck_singular
            end if
            goto 900
        end if
        if (mode == 0) then
            if (present(t)) goto 900
            if (target > 0.0_dp) then
                call fitpack_smooth_parametric(param,x,w,deg,target,[left,right], &
                    tck%t,tck%c,code)
                if (code == 0) then
                    call move_alloc(param,tck%u)
                    tck%k=deg
                    if (.not.tck%valid()) code=fitpack_tck_singular
                end if
                goto 900
            end if
            if (left /= param(1) .or. right /= param(n)) then
                code = fitpack_tck_invalid_input
                goto 900
            end if
        else
            if (.not. present(t)) goto 900
            if (size(t) < 2*deg+2) goto 900
            if (.not. all(ieee_is_finite(t))) goto 900
            if (any(t(2:) < t(:size(t)-1))) goto 900
            if (any(t(:deg+1) /= left)) goto 900
            if (any(t(size(t)-deg:) /= right)) goto 900
            if (any(t(deg+2:size(t)-deg-1) <= left) .or. &
                any(t(deg+2:size(t)-deg-1) >= right)) goto 900
            interior = t(deg+2:size(t)-deg-1)
        end if
        do j=1,ndim
            if (mode == 0) then
                fit = interpolated_univariate_spline(param,x(:,j),k=deg,status=fitcode)
            else
                fit = lsq_univariate_spline(param,x(:,j),interior,w=w,k=deg, &
                    bbox=[left,right],status=fitcode)
            end if
            if (fitcode /= bspline_status_success .or. .not. fit%valid()) then
                code = fitpack_tck_singular
                goto 900
            end if
            if (j == 1) then
                knots = fit%knots()
                nc = size(knots)-deg-1
                allocate(tck%c(nc,ndim))
            end if
            coefficients = fit%coefficients()
            if (size(coefficients) /= nc) then
                code = fitpack_tck_singular
                goto 900
            end if
            tck%c(:,j) = coefficients
        end do
        tck%t = knots
        tck%u = param
        tck%k = deg
        if (.not. tck%valid()) then
            code = fitpack_tck_singular
            goto 900
        end if
        code = fitpack_tck_success
900     if (code /= fitpack_tck_success) then
            if (allocated(tck%c)) deallocate(tck%c)
            if (allocated(tck%t)) deallocate(tck%t)
            if (allocated(tck%u)) deallocate(tck%u)
            tck%k = -1
        end if
        if (present(status)) status = code
    end function splprep

    function splev_parametric_scalar(x,tck,der,ext,status) result(v)
        real(dp), intent(in) :: x !! parameter at which to evaluate curve
        type(fitpack_parametric_tck_type), intent(in) :: tck !! parametric spline
        integer, intent(in), optional :: der !! derivative order
        integer, intent(in), optional :: ext !! outside-domain policy 0:3
        integer, intent(out), optional :: status !! evaluation status
        real(dp), allocatable :: v(:)
        real(dp), allocatable :: matrix(:,:)
        integer :: code
        matrix = splev_parametric_array([x],tck,der,ext,code)
        v = matrix(1,:)
        if (present(status)) status = code
    end function splev_parametric_scalar

    function splev_parametric_array(x,tck,der,ext,status) result(v)
        real(dp), intent(in) :: x(:) !! parameter evaluation sites
        type(fitpack_parametric_tck_type), intent(in) :: tck !! parametric spline
        integer, intent(in), optional :: der !! derivative order
        integer, intent(in), optional :: ext !! FITPACK outside-domain policy
        integer, intent(out), optional :: status !! success, invalid or out-of-bounds
        real(dp), allocatable :: v(:,:) !! (evaluations, dimensions)
        type(fitpack_tck_type) :: single
        integer :: j, code
        if (.not. tck%valid()) then
            allocate(v(size(x),0))
            if (present(status)) status = fitpack_tck_invalid_input
            return
        end if
        allocate(v(size(x),size(tck%c,2)))
        single%t = tck%t
        single%k = tck%k
        do j=1,size(tck%c,2)
            single%c = tck%c(:,j)
            v(:,j) = splev_array(x,single,der,ext,code)
            if (code /= fitpack_tck_success) exit
        end do
        if (code /= fitpack_tck_success) v = ieee_value(0.0_dp,ieee_quiet_nan)
        if (present(status)) status = code
    end function splev_parametric_array

    pure logical function tck_valid(self) result(ok)
        class(fitpack_tck_type), intent(in) :: self !! tuple to validate
        integer :: nc
        ok = allocated(self%t) .and. allocated(self%c) .and. self%k >= 0
        if (.not. ok) return
        nc = size(self%t)-self%k-1
        ok = nc >= self%k+1 .and. size(self%c) >= nc
        if (.not. ok) return
        ok = all(ieee_is_finite(self%t)) .and. all(ieee_is_finite(self%c(:nc)))
        if (.not. ok) return
        ok = all(self%t(2:) >= self%t(:size(self%t)-1))
        if (.not. ok) return
        ok = self%t(self%k+1) < self%t(nc+1)
    end function tck_valid

    pure integer function tck_ncoeff(self) result(nc)
        class(fitpack_tck_type), intent(in) :: self !! tuple to query
        nc = 0
        if (self%valid()) nc = size(self%t)-self%k-1
    end function tck_ncoeff

    function pack_tck(s) result(tck)
        type(bspline_interpolator_type), intent(in) :: s !! fitted B-spline
        type(fitpack_tck_type) :: tck
        real(dp), allocatable :: values(:)
        if (.not. s%valid()) return
        tck%t = s%knots()
        values = s%coefficients()
        tck%c = [values,spread(0.0_dp,1,s%degree()+1)]
        tck%k = s%degree()
    end function pack_tck

    ! Task 0: interpolation at s=0 or automatic-knot smoothing at s>0.
    ! Task -1: specified knots and LSQ.
    function splrep(x,y,w,k,s,task,t,bbox,status,per) result(tck)
        implicit none
        real(dp), intent(in) :: x(:) !! increasing sample abscissae
        real(dp), intent(in) :: y(:) !! scalar sample values
        real(dp), intent(in), optional :: w(:) !! positive FITPACK residual multipliers
        integer, intent(in), optional :: k !! spline degree 1..5 (default 3)
        real(dp), intent(in), optional :: s !! weighted residual target (s=0 for interpolation)
        integer, intent(in), optional :: task !! 0 interpolation/smoothing, -1 fixed knots LSQ
        real(dp), intent(in), optional :: t(:) !! interior knots when task=-1
        real(dp), intent(in), optional :: bbox(:) !! left and right support bounds
        integer, intent(out), optional :: status !! success, invalid, unsupported, singular
        logical, intent(in), optional :: per !! Periodic fit over x(1):x(m), default false.
        type(fitpack_tck_type) :: tck
        type(bspline_interpolator_type) :: fit
        real(dp) :: target, bounds(2)
        integer :: deg, mode, code, used
        logical :: periodic
        periodic=.false.
        if (present(per)) periodic=per
        deg = 3
        if (present(k)) deg = k
        mode = 0
        if (present(task)) mode = task
        code = fitpack_tck_invalid_input
        if (deg < 1 .or. deg > 5) goto 900
        if (size(x)<deg+1 .or. size(y)/=size(x)) goto 900
        used=size(x)
        if (periodic) used=used-1
        if (.not.all(ieee_is_finite(x))) goto 900
        if (.not.all(ieee_is_finite(y(:used)))) goto 900
        if (any(x(2:) <= x(:size(x)-1))) goto 900
        if (present(w)) then
            if (size(w) /= size(x)) goto 900
            if (.not. all(ieee_is_finite(w(:used)))) goto 900
            if (any(w(:used) <= 0.0_dp)) goto 900
        end if
        target = 0.0_dp
        if (present(s)) target = s
        if (.not. ieee_is_finite(target)) goto 900
        if (target < 0.0_dp) goto 900
        bounds = [x(1),x(size(x))]
        if (present(bbox)) then
            if (size(bbox) /= 2) goto 900
            if (.not. all(ieee_is_finite(bbox))) goto 900
            if (bbox(1) > x(1) .or. bbox(2) < x(size(x))) goto 900
            bounds = bbox
        end if
        if (periodic) then
            if (mode/=0) then
                if (mode==-1 .and. present(t)) code=fitpack_tck_unsupported_smoothing
                goto 900
            end if
            if (present(t)) goto 900
            if (any(bounds/=[x(1),x(size(x))])) goto 900
            call fitpack_periodic_curve(x,y,w,deg,target,tck%t,tck%c,code)
            if (code==0) tck%k=deg
            goto 900
        end if
        if (mode == 0) then
            if (present(t)) goto 900
            if (target > 0.0_dp) then
                call fitpack_smooth_curve(x,y,w,deg,target,bounds,tck%t,tck%c,code)
                if (code == 0) tck%k = deg
                goto 900
            end if
            if (present(bbox)) then
                if (size(bbox) /= 2) goto 900
                ! Exact interpolation with enlarged support is not supported.
                if (any(abs(bbox-[x(1),x(size(x))]) > 0.0_dp)) goto 900
            end if
            fit = interpolated_univariate_spline(x,y,k=deg,status=code)
        else if (mode == -1) then
            if (.not. present(t)) goto 900
            fit = lsq_univariate_spline(x,y,t,w=w,k=deg,bbox=bbox,status=code)
        else
            goto 900
        end if
        if (code == bspline_status_success) then
            tck = pack_tck(fit)
            code = fitpack_tck_success
        else if (code /= fitpack_tck_invalid_input) then
            code = fitpack_tck_singular
        end if
900     if (present(status)) status = code
    end function splrep

    function make_spline(tck,code) result(s)
        type(fitpack_tck_type), intent(in) :: tck !! SciPy tuple
        integer, intent(out) :: code !! construction status
        type(bspline_interpolator_type) :: s
        code = fitpack_tck_invalid_input
        if (.not. tck%valid()) return
        s = bspline(tck%t,tck%c,tck%k,status=code)
    end function make_spline

    function splev_scalar(x,tck,der,ext,status) result(v)
        real(dp), intent(in) :: x !! evaluation coordinate
        type(fitpack_tck_type), intent(in) :: tck !! tuple
        integer, intent(in), optional :: der !! derivative order (default 0)
        integer, intent(in), optional :: ext !! 0 extrapolate,1 zero,2 error,3 constant
        integer, intent(out), optional :: status !! evaluation status
        real(dp) :: v
        real(dp) :: tmp(1)
        integer :: code
        tmp = splev_array([x],tck,der,ext,code)
        v = tmp(1)
        if (present(status)) status=code
    end function splev_scalar

    function splev_array(x,tck,der,ext,status) result(v)
        real(dp), intent(in) :: x(:) !! evaluation coordinates
        type(fitpack_tck_type), intent(in) :: tck !! tuple
        integer, intent(in), optional :: der !! derivative order (default 0)
        integer, intent(in), optional :: ext !! outside-domain policy
        integer, intent(out), optional :: status !! success or invalid/out-of-bounds
        real(dp) :: v(size(x))
        type(bspline_interpolator_type) :: spline_obj
        integer :: code, order, policy, i
        real(dp) :: lo, hi, xi
        v = ieee_value(0.0_dp,ieee_quiet_nan)
        code = fitpack_tck_invalid_input
        order = 0
        if (present(der)) order = der
        policy = 0
        if (present(ext)) policy=ext
        if (order < 0 .or. order > tck%k) goto 900
        if (policy < 0 .or. policy > 3) goto 900
        if (.not. all(ieee_is_finite(x))) goto 900
        spline_obj = make_spline(tck,code)
        if (code /= fitpack_tck_success) goto 900
        lo = tck%t(tck%k+1)
        hi = tck%t(tck%ncoeff()+1)
        if (policy == 2) then
            if (any(x < lo) .or. any(x > hi)) then
                code = fitpack_tck_out_of_bounds
                goto 900
            end if
        end if
        do i=1,size(x)
            xi = x(i)
            if (xi < lo .or. xi > hi) then
                select case(policy)
                case(1)
                    v(i) = 0.0_dp
                    cycle
                case(3)
                    xi = min(max(xi,lo),hi)
                end select
            end if
            if (order == 0) then
                v(i) = spline_obj%evaluate(xi)
            else
                v(i) = spline_obj%derivative(xi,nu=order)
            end if
        end do
900     if (present(status)) status=code
    end function splev_array

    function splint(a,b,tck,status) result(v)
        real(dp), intent(in) :: a !! lower bound
        real(dp), intent(in) :: b !! upper bound
        type(fitpack_tck_type), intent(in) :: tck !! tuple
        integer, intent(out), optional :: status !! success or invalid
        real(dp) :: v,lo,hi,aa,bb
        integer :: code
        type(bspline_interpolator_type) :: spline_obj
        v=ieee_value(0.0_dp,ieee_quiet_nan)
        spline_obj=make_spline(tck,code)
        if (code /= fitpack_tck_success) goto 900
        if (.not. ieee_is_finite(a) .or. .not. ieee_is_finite(b)) then
            code = fitpack_tck_invalid_input
            goto 900
        end if
        lo=tck%t(tck%k+1)
        hi=tck%t(tck%ncoeff()+1)
        aa=max(lo,min(hi,a))
        bb=max(lo,min(hi,b))
        v = spline_obj%integrate(aa,bb)
900     if (present(status)) status=code
    end function splint

    function splder(tck,n,status) result(out)
        type(fitpack_tck_type), intent(in) :: tck !! input tuple
        integer, intent(in), optional :: n !! derivative order, default 1
        integer, intent(out), optional :: status !! status code
        type(fitpack_tck_type) :: out
        type(bspline_interpolator_type) :: a,b
        integer :: code,order
        order = 1
        if (present(n)) order=n
        a=make_spline(tck,code)
        if (code /= fitpack_tck_success) goto 900
        b=a%derivative_spline(nu=order,status=code)
        if (code == bspline_status_success) out=pack_tck(b)
900     if (present(status)) status=code
    end function splder

    function splantider(tck,n,status) result(out)
        type(fitpack_tck_type), intent(in) :: tck !! input tuple
        integer, intent(in), optional :: n !! antiderivative order, default 1
        integer, intent(out), optional :: status !! status code
        type(fitpack_tck_type) :: out
        type(bspline_interpolator_type) :: a,b
        integer :: code,order
        order = 1
        if (present(n)) order=n
        a=make_spline(tck,code)
        if (code /= fitpack_tck_success) goto 900
        b=a%antiderivative(nu=order,status=code)
        if (code == bspline_status_success) out=pack_tck(b)
900     if (present(status)) status=code
    end function splantider

    function spalde(x,tck,status) result(v)
        real(dp), intent(in) :: x !! evaluation coordinate, within domain
        type(fitpack_tck_type), intent(in) :: tck !! tuple
        integer, intent(out), optional :: status !! result status
        real(dp), allocatable :: v(:)
        integer :: i,code
        allocate(v(0))
        code=fitpack_tck_invalid_input
        if (.not. tck%valid()) goto 900
        if (.not. ieee_is_finite(x)) goto 900
        if (x < tck%t(tck%k+1) .or. x > tck%t(tck%ncoeff()+1)) then
            code=fitpack_tck_out_of_bounds
            goto 900
        end if
        deallocate(v)
        allocate(v(0:tck%k))
        do i=0,tck%k
            v(i)=splev_scalar(x,tck,der=i,status=code)
            if (code /= fitpack_tck_success) exit
        end do
900     if (present(status)) status=code
    end function spalde
end module scifort_interpolate_fitpack_compat
