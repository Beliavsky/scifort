! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! General, nonuniform B-splines: de Boor evaluation, analytic derivatives and
! antiderivatives, interpolation by collocation, and knot insertion.
! Independent implementation; SciPy is an external numerical test oracle.
module scifort_interpolate_bspline
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
    use scifort_kinds, only: dp
    implicit none
    private

    integer, parameter, public :: bspline_status_success = 0
    integer, parameter, public :: bspline_status_invalid_input = 1
    integer, parameter, public :: bspline_status_singular = 2

    type, public :: bspline_interpolator_type
        private
        real(dp), allocatable :: t(:)
        real(dp), allocatable :: c(:)
        integer :: k = -1
        logical :: extrapolate = .true.
        logical :: periodic = .false.
    contains
        procedure :: valid => spline_valid
        procedure :: knots => spline_knots
        procedure :: coefficients => spline_coefficients
        procedure :: degree => spline_degree
        procedure, private :: spline_eval_scalar
        procedure, private :: spline_eval_vector
        generic :: evaluate => spline_eval_scalar, spline_eval_vector
        procedure, private :: spline_deriv_scalar
        procedure, private :: spline_deriv_vector
        generic :: derivative => spline_deriv_scalar, spline_deriv_vector
        procedure :: derivative_spline => spline_derivative_spline
        procedure :: antiderivative => spline_antiderivative
        procedure :: integrate => spline_integrate
        procedure :: insert_knot => spline_insert_knot
    end type bspline_interpolator_type

    type, public :: bspline_array_interpolator_type
        private
        type(bspline_interpolator_type), allocatable :: columns(:)
    contains
        procedure :: valid => array_valid
        procedure :: channels => array_channels
        procedure, private :: array_eval_scalar
        procedure, private :: array_eval_vector
        generic :: evaluate => array_eval_scalar, array_eval_vector
        procedure, private :: array_deriv_scalar
        procedure, private :: array_deriv_vector
        generic :: derivative => array_deriv_scalar, array_deriv_vector
        procedure :: integrate => array_integrate
        procedure :: knots => array_knots
        procedure :: coefficients => array_coefficients
        procedure :: insert_knot => array_insert_knot
    end type bspline_array_interpolator_type

    public :: bspline, make_interp_spline, make_lsq_spline, bspline_design_matrix
    public :: lsq_univariate_spline, interpolated_univariate_spline

    interface bspline
        module procedure create_bspline
        module procedure create_bspline_array
    end interface

    interface make_lsq_spline
        module procedure make_lsq_spline_scalar
        module procedure make_lsq_spline_channels
    end interface

    interface lsq_univariate_spline
        module procedure lsq_spline_scalar
        module procedure lsq_spline_matrix
    end interface

    interface interpolated_univariate_spline
        module procedure interpolated_spline_scalar
        module procedure interpolated_spline_matrix
    end interface

    interface make_interp_spline
        module procedure make_interp_spline_scalar
        module procedure make_interp_spline_array
    end interface

contains

    pure function nan_value() result(v)
        real(dp) :: v
        v = ieee_value(0.0_dp, ieee_quiet_nan)
    end function nan_value

    pure logical function spline_valid(self) result(ok)
        class(bspline_interpolator_type), intent(in) :: self !! spline to validate
        ok = allocated(self%t) .and. allocated(self%c) .and. self%k >= 0
        if (.not. ok) return
        ok = size(self%t) == size(self%c)+self%k+1
    end function spline_valid

    pure logical function array_valid(self) result(ok)
        class(bspline_array_interpolator_type), intent(in) :: self !! array spline to validate
        integer :: j
        ok = allocated(self%columns)
        if (.not. ok) return
        ok = size(self%columns) > 0
        if (.not. ok) return
        do j = 1, size(self%columns)
            if (.not. self%columns(j)%valid()) ok = .false.
        end do
    end function array_valid

    pure integer function array_channels(self) result(n)
        class(bspline_array_interpolator_type), intent(in) :: self !! possibly uninitialized spline
        n = 0
        if (allocated(self%columns)) n = size(self%columns)
    end function array_channels

    pure integer function spline_degree(self) result(k)
        class(bspline_interpolator_type), intent(in) :: self !! spline
        k = self%k
    end function spline_degree

    function spline_knots(self) result(t)
        class(bspline_interpolator_type), intent(in) :: self !! spline
        real(dp), allocatable :: t(:)
        if (allocated(self%t)) then
            t = self%t
        else
            allocate(t(0))
        end if
    end function spline_knots

    function spline_coefficients(self) result(c)
        class(bspline_interpolator_type), intent(in) :: self !! spline
        real(dp), allocatable :: c(:)
        if (allocated(self%c)) then
            c = self%c
        else
            allocate(c(0))
        end if
    end function spline_coefficients

    function array_knots(self) result(t)
        class(bspline_array_interpolator_type), intent(in) :: self !! multichannel spline
        real(dp), allocatable :: t(:)
        allocate(t(0))
        if (self%valid()) t = self%columns(1)%knots()
    end function array_knots

    function array_coefficients(self) result(c)
        class(bspline_array_interpolator_type), intent(in) :: self !! multichannel spline
        real(dp), allocatable :: c(:,:)
        integer :: j
        allocate(c(0,0))
        if (.not. self%valid()) return
        deallocate(c)
        allocate(c(size(self%columns(1)%c),size(self%columns)))
        do j = 1,size(self%columns)
            c(:,j) = self%columns(j)%c
        end do
    end function array_coefficients

    ! SciPy's BSpline(t,c,k): c may contain extra trailing entries, ignored.
    function create_bspline(t, c, k, status, extrapolate, periodic) result(s)
        real(dp), intent(in) :: t(:) !! nondecreasing finite knots, length >= 2*k+2
        real(dp), intent(in) :: c(:) !! real coefficients, at least size(t)-k-1
        integer, intent(in) :: k !! polynomial degree, nonnegative
        integer, intent(out), optional :: status !! 0 success, 1 invalid input
        logical, intent(in), optional :: extrapolate !! extend endpoint polynomial pieces (default true)
        logical, intent(in), optional :: periodic !! periodically wrap queries instead (default false)
        type(bspline_interpolator_type) :: s
        integer :: nc, code

        code = bspline_status_invalid_input
        if (k < 0) goto 900
        if (size(t) < 2*k+2) goto 900
        nc = size(t)-k-1
        if (size(c) < nc) goto 900
        if (.not. all(ieee_is_finite(t))) goto 900
        if (.not. all(ieee_is_finite(c(1:nc)))) goto 900
        if (any(t(2:) < t(:size(t)-1))) goto 900
        if (t(k+1) >= t(nc+1)) goto 900
        s%t = t
        s%c = c(1:nc)
        s%k = k
        if (present(extrapolate)) s%extrapolate = extrapolate
        if (present(periodic)) s%periodic = periodic
        code = bspline_status_success
900     if (present(status)) status = code
    end function create_bspline

    function create_bspline_array(t, c, k, status, extrapolate, periodic) result(s)
        real(dp), intent(in) :: t(:) !! shared nondecreasing knots
        real(dp), intent(in) :: c(:,:) !! coefficients: (coefficient, channel)
        integer, intent(in) :: k !! polynomial degree
        integer, intent(out), optional :: status !! 0 success, 1 invalid input
        logical, intent(in), optional :: extrapolate !! allow endpoint extrapolation
        logical, intent(in), optional :: periodic !! periodic evaluation
        type(bspline_array_interpolator_type) :: s
        integer :: j, code
        code = bspline_status_invalid_input
        if (size(c,2) < 1) goto 900
        allocate(s%columns(size(c,2)))
        do j=1,size(c,2)
            s%columns(j) = create_bspline(t,c(:,j),k,code,extrapolate,periodic)
            if (code /= bspline_status_success) then
                deallocate(s%columns)
                goto 900
            end if
        end do
900     if (present(status)) status = code
    end function create_bspline_array

    ! Right-continuous knot span; the right boundary belongs to the final span.
    pure integer function find_span(t, k, x) result(s)
        real(dp), intent(in) :: t(:) !! nondecreasing knots
        integer, intent(in) :: k !! spline degree
        real(dp), intent(in) :: x !! finite evaluation coordinate
        integer :: lo, hi, mid, nc
        nc = size(t)-k-1
        if (x < t(k+1)) then
            s = k+1
            return
        end if
        if (x >= t(nc+1)) then
            s = nc
            return
        end if
        lo = k+1
        hi = nc+1
        do while (hi-lo > 1)
            mid = (hi+lo)/2
            if (x < t(mid)) then
                hi = mid
            else
                lo = mid
            end if
        end do
        s = lo
    end function find_span

    ! Nonzero B-spline basis functions on the selected span, N(0:k).
    pure subroutine basis_values(t, k, span, x, vals)
        real(dp), intent(in) :: t(:) !! knots
        integer, intent(in) :: k !! degree
        integer, intent(in) :: span !! current span, 1-based
        real(dp), intent(in) :: x !! evaluation point
        real(dp), intent(out) :: vals(0:) !! basis values (0:k)
        real(dp) :: left(k), right(k), den, saved, temp
        integer :: j, r
        vals = 0.0_dp
        vals(0) = 1.0_dp
        do j = 1,k
            left(j) = x-t(span+1-j)
            right(j) = t(span+j)-x
            saved = 0.0_dp
            do r = 0,j-1
                den = right(r+1)+left(j-r)
                temp = 0.0_dp
                if (den /= 0.0_dp) temp = vals(r)/den
                vals(r) = saved+right(r+1)*temp
                saved = left(j-r)*temp
            end do
            vals(j) = saved
        end do
    end subroutine basis_values

    pure function eval_bspline(t, c, k, x, extrapolate, periodic) result(v)
        real(dp), intent(in) :: t(:) !! knots, length size(c)+k+1
        real(dp), intent(in) :: c(:) !! spline coefficients
        integer, intent(in) :: k !! degree
        real(dp), intent(in) :: x !! query point
        logical, intent(in) :: extrapolate !! endpoint polynomial extrapolation
        logical, intent(in) :: periodic !! wrap query to the base interval
        real(dp) :: v, xx, den, alpha
        real(dp) :: d(0:k)
        integer :: s, j, r, idx, nc
        nc = size(c)
        v = nan_value()
        if (.not. ieee_is_finite(x)) return
        xx = x
        if (periodic) xx = t(k+1)+modulo(x-t(k+1),t(nc+1)-t(k+1))
        if (.not. extrapolate .and. .not. periodic) then
            if (xx < t(k+1) .or. xx > t(nc+1)) return
        end if
        s = find_span(t,k,xx)
        do j=0,k
            d(j) = c(s-k+j)
        end do
        do r=1,k
            do j=k,r,-1
                idx = s-k+j
                den = t(idx+k-r+1)-t(idx)
                alpha = 0.0_dp
                if (den /= 0.0_dp) alpha = (xx-t(idx))/den
                d(j) = (1.0_dp-alpha)*d(j-1)+alpha*d(j)
            end do
        end do
        v = d(k)
    end function eval_bspline

    pure real(dp) function spline_eval_scalar(self,x) result(v)
        class(bspline_interpolator_type), intent(in) :: self !! spline
        real(dp), intent(in) :: x !! evaluation point
        v = nan_value()
        if (.not. self%valid()) return
        v = eval_bspline(self%t,self%c,self%k,x,self%extrapolate,self%periodic)
    end function spline_eval_scalar

    pure function spline_eval_vector(self,x) result(v)
        class(bspline_interpolator_type), intent(in) :: self !! spline
        real(dp), intent(in) :: x(:) !! evaluation points
        real(dp) :: v(size(x))
        integer :: i
        do i=1,size(x)
            v(i) = self%spline_eval_scalar(x(i))
        end do
    end function spline_eval_vector

    pure function array_eval_scalar(self,x) result(v)
        class(bspline_array_interpolator_type), intent(in) :: self !! multichannel spline
        real(dp), intent(in) :: x !! evaluation point
        real(dp), allocatable :: v(:)
        integer :: j
        allocate(v(self%channels()))
        do j=1,size(v)
            v(j) = self%columns(j)%evaluate(x)
        end do
    end function array_eval_scalar

    pure function array_eval_vector(self,x) result(v)
        class(bspline_array_interpolator_type), intent(in) :: self !! multichannel spline
        real(dp), intent(in) :: x(:) !! evaluation points
        real(dp), allocatable :: v(:,:)
        integer :: j
        allocate(v(size(x),self%channels()))
        do j=1,size(v,2)
            v(:,j) = self%columns(j)%evaluate(x)
        end do
    end function array_eval_vector

    ! Derivative in B-spline form: reduce degree and drop outermost knots.
    function spline_derivative_spline(self, nu, status) result(s)
        class(bspline_interpolator_type), intent(in) :: self !! spline to differentiate
        integer, intent(in), optional :: nu !! nonnegative derivative order, default 1
        integer, intent(out), optional :: status !! success / invalid input
        type(bspline_interpolator_type) :: s
        real(dp), allocatable :: t(:), c(:), d(:)
        integer :: r, i, deg, order, code
        real(dp) :: denom
        code = bspline_status_invalid_input
        if (.not. self%valid()) goto 900
        order = 1
        if (present(nu)) order = nu
        if (order < 0 .or. order > self%k) goto 900
        t = self%t
        c = self%c
        deg = self%k
        do r=1,order
            allocate(d(size(c)-1))
            do i=1,size(d)
                denom = t(i+deg+1)-t(i+1)
                ! SciPy rejects derivative objects at knots of multiplicity k+1.
                if (denom == 0.0_dp) goto 900
                d(i) = real(deg,dp)*(c(i+1)-c(i))/denom
            end do
            c = d
            deallocate(d)
            t = t(2:size(t)-1)
            deg = deg-1
        end do
        s = create_bspline(t,c,deg,status=code,extrapolate=self%extrapolate,periodic=self%periodic)
900     if (present(status)) status = code
    end function spline_derivative_spline

    pure real(dp) function spline_deriv_scalar(self,x,nu) result(v)
        class(bspline_interpolator_type), intent(in) :: self !! spline
        real(dp), intent(in) :: x !! evaluation point
        integer, intent(in), optional :: nu !! derivative order, default 1
        real(dp), allocatable :: t(:), c(:), d(:)
        integer :: order, r, deg, i
        real(dp) :: den
        v = nan_value()
        if (.not. self%valid()) return
        order = 1
        if (present(nu)) order = nu
        if (order < 0) return
        if (order == 0) then
            v = self%spline_eval_scalar(x)
            return
        end if
        if (order > self%k) then
            v = 0.0_dp
            if (order == self%k+1 .and. .not. self%extrapolate .and. .not. self%periodic) then
                if (x < self%t(self%k+1) .or. x > self%t(size(self%c)+1)) v = nan_value()
            end if
            if (.not. ieee_is_finite(x)) v = nan_value()
            return
        end if
        t = self%t
        c = self%c
        deg = self%k
        do r=1,order
            allocate(d(size(c)-1))
            do i=1,size(d)
                den = t(i+deg+1)-t(i+1)
                d(i) = 0.0_dp
                if (den /= 0.0_dp) d(i) = real(deg,dp)*(c(i+1)-c(i))/den
            end do
            c = d
            deallocate(d)
            t = t(2:size(t)-1)
            deg = deg-1
        end do
        v = eval_bspline(t,c,deg,x,self%extrapolate,self%periodic)
    end function spline_deriv_scalar

    pure function spline_deriv_vector(self,x,nu) result(v)
        class(bspline_interpolator_type), intent(in) :: self !! spline
        real(dp), intent(in) :: x(:) !! evaluation points
        integer, intent(in), optional :: nu !! derivative order
        real(dp) :: v(size(x))
        integer :: i
        do i=1,size(x)
            v(i) = self%spline_deriv_scalar(x(i),nu)
        end do
    end function spline_deriv_vector

    pure function array_deriv_scalar(self,x,nu) result(v)
        class(bspline_array_interpolator_type), intent(in) :: self !! spline
        real(dp), intent(in) :: x !! query point
        integer, intent(in), optional :: nu !! derivative order
        real(dp), allocatable :: v(:)
        integer :: j
        allocate(v(self%channels()))
        do j=1,size(v)
            v(j) = self%columns(j)%derivative(x,nu)
        end do
    end function array_deriv_scalar

    pure function array_deriv_vector(self,x,nu) result(v)
        class(bspline_array_interpolator_type), intent(in) :: self !! spline
        real(dp), intent(in) :: x(:) !! evaluation points
        integer, intent(in), optional :: nu !! derivative order
        real(dp), allocatable :: v(:,:)
        integer :: j
        allocate(v(size(x),self%channels()))
        do j=1,size(v,2)
            v(:,j) = self%columns(j)%derivative(x,nu)
        end do
    end function array_deriv_vector

    function spline_antiderivative(self, nu, status) result(s)
        class(bspline_interpolator_type), intent(in) :: self !! spline to integrate
        integer, intent(in), optional :: nu !! positive integral count (default 1)
        integer, intent(out), optional :: status !! status
        type(bspline_interpolator_type) :: s
        real(dp), allocatable :: t(:), c(:), d(:)
        integer :: i, r, order, deg, code
        code = bspline_status_invalid_input
        if (.not. self%valid()) goto 900
        order = 1
        if (present(nu)) order = nu
        if (order < 0 .or. order > 20) goto 900
        t = self%t
        c = self%c
        deg = self%k
        do r = 1,order
            allocate(d(size(c)+1))
            d(1) = 0.0_dp
            do i = 1,size(c)
                d(i+1) = d(i)+c(i)*(t(i+deg+1)-t(i))/real(deg+1,dp)
            end do
            c = d
            deallocate(d)
            t = [t(1),t,t(size(t))]
            deg = deg+1
        end do
        s = create_bspline(t,c,deg,status=code, &
            extrapolate=self%extrapolate .and. .not. self%periodic,periodic=.false.)
900     if (present(status)) status = code
    end function spline_antiderivative

    recursive pure real(dp) function spline_integrate(self, a, b) result(v)
        class(bspline_interpolator_type), intent(in) :: self !! spline
        real(dp), intent(in) :: a !! integration lower limit
        real(dp), intent(in) :: b !! integration upper limit
        real(dp), allocatable :: t(:), c(:), d(:)
        real(dp) :: left, right, lo, hi, length, period, aa, bb, total, fleft, fhi
        integer :: i, deg
        v = nan_value()
        if (.not. self%valid()) return
        if (.not. ieee_is_finite(a) .or. .not. ieee_is_finite(b)) return
        if (a == b) then
            v = 0.0_dp
            if (.not. self%extrapolate .and. .not. self%periodic) then
                if (a < self%t(self%k+1) .or. a > self%t(size(self%c)+1)) v = nan_value()
            end if
            return
        end if
        if (b < a) then
            v = -self%integrate(b,a)
            return
        end if
        lo = self%t(self%k+1)
        hi = self%t(size(self%c)+1)
        aa = a
        bb = b
        if (.not. self%extrapolate .and. .not. self%periodic) then
            aa = max(lo,a)
            bb = min(hi,b)
            if (aa >= bb) then
                v = 0.0_dp
                return
            end if
        end if
        deg = self%k
        t = self%t
        c = self%c
        allocate(d(size(c)+1))
        d(1) = 0.0_dp
        do i=1,size(c)
            d(i+1) = d(i)+c(i)*(t(i+deg+1)-t(i))/real(deg+1,dp)
        end do
        t = [t(1),t,t(size(t))]
        if (self%periodic) then
            length = hi-lo
            total = eval_bspline(t,d,deg+1,hi,.true.,.false.) &
                -eval_bspline(t,d,deg+1,lo,.true.,.false.)
            ! Integral of periodically repeated base polynomial pieces.
            period = floor((bb-lo)/length,kind=8)-floor((aa-lo)/length,kind=8)
            left = lo+modulo(aa-lo,length)
            right = lo+modulo(bb-lo,length)
            fleft = eval_bspline(t,d,deg+1,left,.true.,.false.)
            fhi = eval_bspline(t,d,deg+1,right,.true.,.false.)
            v = period*total+fhi-fleft
        else
            v = eval_bspline(t,d,deg+1,bb,.true.,.false.) &
                -eval_bspline(t,d,deg+1,aa,.true.,.false.)
        end if
    end function spline_integrate

    pure function array_integrate(self,a,b) result(v)
        class(bspline_array_interpolator_type), intent(in) :: self !! spline
        real(dp), intent(in) :: a !! lower integration limit
        real(dp), intent(in) :: b !! upper integration limit
        real(dp), allocatable :: v(:)
        integer :: j
        allocate(v(self%channels()))
        do j=1,size(v)
            v(j) = self%columns(j)%integrate(a,b)
        end do
    end function array_integrate

    ! Boehm's shape-preserving one-knot insertion, repeated m times.
    function spline_insert_knot(self, x, m, status) result(s)
        class(bspline_interpolator_type), intent(in) :: self !! original spline
        real(dp), intent(in) :: x !! new knot within the base interval
        integer, intent(in), optional :: m !! number of insertions, default 1
        integer, intent(out), optional :: status !! input validation status
        type(bspline_interpolator_type) :: s
        real(dp), allocatable :: t(:), c(:), tt(:), cc(:)
        real(dp) :: den, alpha
        integer :: reps, r, span, multiplicity, j, n, k, code
        code = bspline_status_invalid_input
        if (.not. self%valid()) goto 900
        if (.not. ieee_is_finite(x)) goto 900
        k = self%k
        n = size(self%c)
        if (x < self%t(k+1) .or. x > self%t(n+1)) goto 900
        reps = 1
        if (present(m)) reps = m
        if (reps < 0) goto 900
        t = self%t
        c = self%c
        do r=1,reps
            span = find_span(t,k,x)
            multiplicity = count(t == x)
            if (multiplicity >= k+1) goto 900
            n = size(c)
            allocate(tt(size(t)+1), cc(n+1))
            tt(1:span) = t(1:span)
            tt(span+1) = x
            tt(span+2:) = t(span+1:)
            do j=1,span-k
                cc(j) = c(j)
            end do
            do j=span-multiplicity+1,n+1
                cc(j) = c(j-1)
            end do
            do j=span-k+1,span-multiplicity
                den = t(j+k)-t(j)
                alpha = 0.0_dp
                if (den /= 0.0_dp) alpha = (x-t(j))/den
                cc(j) = alpha*c(j)+(1.0_dp-alpha)*c(j-1)
            end do
            t = tt
            c = cc
            deallocate(tt,cc)
        end do
        s = create_bspline(t,c,k,status=code,extrapolate=self%extrapolate,periodic=self%periodic)
900     if (present(status)) status = code
    end function spline_insert_knot

    function array_insert_knot(self,x,m,status) result(s)
        class(bspline_array_interpolator_type), intent(in) :: self !! original spline
        real(dp), intent(in) :: x !! insertion knot
        integer, intent(in), optional :: m !! repetitions, default 1
        integer, intent(out), optional :: status !! status
        type(bspline_array_interpolator_type) :: s
        integer :: j, code
        code = bspline_status_invalid_input
        if (.not. self%valid()) goto 900
        allocate(s%columns(size(self%columns)))
        do j=1,size(s%columns)
            s%columns(j) = self%columns(j)%insert_knot(x,m,status=code)
            if (code /= bspline_status_success) then
                deallocate(s%columns)
                goto 900
            end if
        end do
900     if (present(status)) status=code
    end function array_insert_knot

    ! Dense design matrix, unlike SciPy's CSR sparse result. Caller may pass
    ! this to a least-squares or banded linear system solver.
    function bspline_design_matrix(x, t, k, status, extrapolate) result(mat)
        real(dp), intent(in) :: x(:) !! evaluation coordinates
        real(dp), intent(in) :: t(:) !! knot vector
        integer, intent(in) :: k !! polynomial degree
        integer, intent(out), optional :: status !! 0 success; 1 invalid
        logical, intent(in), optional :: extrapolate !! evaluate outside the base interval (default false)
        real(dp), allocatable :: mat(:,:)
        type(bspline_interpolator_type) :: s
        real(dp), allocatable :: c(:)
        real(dp) :: v(0:max(0,k))
        integer :: i, j, span, nc, code
        logical :: allow
        code = bspline_status_invalid_input
        allocate(mat(0,0))
        if (k < 0 .or. size(t) < 2*k+2) goto 900
        nc = size(t)-k-1
        allocate(c(nc))
        c = 0.0_dp
        s = create_bspline(t,c,k,status=code)
        if (code /= bspline_status_success) goto 900
        code = bspline_status_invalid_input
        if (.not. all(ieee_is_finite(x))) goto 900
        allow = .false.
        if (present(extrapolate)) allow = extrapolate
        if (.not. allow) then
            if (any(x < t(k+1)) .or. any(x > t(nc+1))) goto 900
        end if
        deallocate(mat)
        allocate(mat(size(x),nc))
        mat = 0.0_dp
        do i=1,size(x)
            span = find_span(t,k,x(i))
            call basis_values(t,k,span,x(i),v)
            do j=0,k
                mat(i,span-k+j) = v(j)
            end do
        end do
        code = bspline_status_success
900     if (present(status)) status = code
    end function bspline_design_matrix

    ! Dense system with partial pivoting. No system BLAS or LAPACK dependency.
    ! This targets moderate 1D collocation systems; future tensor backends
    ! can use a banded solver for larger sample counts.
    subroutine solve_collocation(mat, rhs, code)
        real(dp), intent(inout) :: mat(:,:) !! square collocation matrix; overwritten by upper factor
        real(dp), intent(inout) :: rhs(:,:) !! right-hand sides; overwritten by solution
        integer, intent(out) :: code !! 0 success; 2 singular
        real(dp), allocatable :: tmp(:)
        real(dp) :: pivot, scale, mult
        integer :: n, nrhs, i, j, p
        n = size(mat,1)
        nrhs = size(rhs,2)
        code = bspline_status_singular
        if (size(mat,2) /= n .or. size(rhs,1) /= n) return
        allocate(tmp(max(n,nrhs)))
        scale = maxval(abs(mat))
        if (scale <= 0.0_dp) return
        do i=1,n
            p = i-1+maxloc(abs(mat(i:n,i)),dim=1)
            pivot = mat(p,i)
            if (abs(pivot) <= 32.0_dp*epsilon(1.0_dp)*scale) return
            if (p /= i) then
                tmp(:n) = mat(i,:)
                mat(i,:) = mat(p,:)
                mat(p,:) = tmp(:n)
                tmp(:nrhs) = rhs(i,:)
                rhs(i,:) = rhs(p,:)
                rhs(p,:) = tmp(:nrhs)
            end if
            do j=i+1,n
                mult = mat(j,i)/mat(i,i)
                mat(j,i) = 0.0_dp
                mat(j,i+1:n) = mat(j,i+1:n)-mult*mat(i,i+1:n)
                rhs(j,:) = rhs(j,:)-mult*rhs(i,:)
            end do
        end do
        do i=n,1,-1
            do j=i+1,n
                rhs(i,:) = rhs(i,:)-mat(i,j)*rhs(j,:)
            end do
            rhs(i,:) = rhs(i,:)/mat(i,i)
        end do
        code = bspline_status_success
    end subroutine solve_collocation

    ! Derivative weights on the entire coefficient vector, used only at the
    ! left/right edges of collocation systems with boundary conditions.
    subroutine boundary_weights(t,k,x,order,w)
        real(dp), intent(in) :: t(:) !! knot vector
        integer, intent(in) :: k !! polynomial degree
        real(dp), intent(in) :: x !! boundary coordinate
        integer, intent(in) :: order !! derivative order (1 or 2)
        real(dp), intent(out) :: w(:) !! weights for each coefficient
        type(bspline_interpolator_type) :: basis
        real(dp), allocatable :: c(:)
        integer :: i, nc
        nc = size(t)-k-1
        allocate(c(nc))
        do i=1,nc
            c = 0.0_dp
            c(i) = 1.0_dp
            basis = create_bspline(t,c,k)
            w(i) = basis%derivative(x,order)
        end do
    end subroutine boundary_weights

    ! Collocation with periodic endpoint derivatives.  The n data conditions
    ! and k-1 matching derivative conditions determine n+k-1 coefficients.
    ! Knot gaps are continued cyclically beyond the physical interval; the
    ! even-degree variant uses shifted interior sites. The knot-placement
    ! convention and endpoint constraints are adapted from SciPy 1.17.0's
    ! scipy/interpolate/_bsplines.py (BSD-3-Clause); see CODE_PROVENANCE.md
    ! and THIRD_PARTY_LICENSES.md. No FITPACK/LAPACK runtime dependency.
    subroutine solve_periodic_collocation(x, y, k, t, coeff, status)
        real(dp), intent(in) :: x(:) !! strictly increasing sample sites, including both endpoints
        real(dp), intent(in) :: y(:,:) !! sampled ordinates (sample,channel), same values at endpoints
        integer, intent(in) :: k !! polynomial degree in 1:15
        real(dp), allocatable, intent(out) :: t(:) !! extended periodic knot sequence
        real(dp), allocatable, intent(out) :: coeff(:,:) !! fitted B-spline coefficients
        integer, intent(out) :: status !! 0 successful, 1 invalid, 2 singular collocation
        real(dp), allocatable :: adjusted(:), spacing(:), lhs(:,:), rhs(:,:), left(:), right(:)
        real(dp) :: basis(0:max(0,k)), scale, row_scale
        integer :: n, nc, i, j, r, span

        status = bspline_status_invalid_input
        n = size(x)
        if (k < 1 .or. k > 15 .or. n < 2) return
        if (size(y,1) /= n .or. size(y,2) < 1) return
        if (.not. all(ieee_is_finite(x)) .or. .not. all(ieee_is_finite(y))) return
        if (any(x(2:) <= x(:n-1))) return
        scale = max(1.0_dp,maxval(abs(y)))
        if (any(abs(y(1,:)-y(n,:)) > 1.0e-12_dp*scale)) return

        allocate(adjusted(n),spacing(n-1))
        adjusted = x
        if (mod(k,2) == 0 .and. n > 2) then
            adjusted(2:n-1) = x(2:n-1)-0.5_dp*(x(2:n-1)-x(1:n-2))
        end if
        spacing = adjusted(2:n)-adjusted(1:n-1)
        nc = n+k-1
        allocate(t(n+2*k),lhs(nc,nc),rhs(nc,size(y,2)),left(nc),right(nc))
        t(k+1:k+n) = adjusted
        do i = 1,k
            j = mod(i-1,n-1)
            t(k+1-i) = t(k+2-i)-spacing(n-1-j)
            t(k+n+i) = t(k+n+i-1)+spacing(j+1)
        end do
        lhs = 0.0_dp
        rhs = 0.0_dp
        do r = 1,k-1
            call boundary_weights(t,k,x(1),r,left)
            call boundary_weights(t,k,x(n),r,right)
            lhs(r,:) = left-right
            ! Higher-order endpoint derivatives may be many orders of
            ! magnitude larger than value rows. Scale each homogeneous
            ! continuity condition before using the shared dense pivot
            ! threshold, or valid high-degree systems appear singular.
            row_scale = maxval(abs(lhs(r,:)))
            if (row_scale <= 0.0_dp) then
                status = bspline_status_singular
                return
            end if
            lhs(r,:) = lhs(r,:)/row_scale
        end do
        do i = 1,n
            span = find_span(t,k,x(i))
            call basis_values(t,k,span,x(i),basis)
            lhs(k-1+i,span-k:span) = basis
        end do
        rhs(k:k+n-1,:) = y
        call solve_collocation(lhs,rhs,status)
        if (status == bspline_status_success) coeff = rhs
    end subroutine solve_periodic_collocation

    function make_interp_spline_scalar(x, y, k, bc_type, status, extrapolate, knots) result(s)
        real(dp), intent(in) :: x(:) !! strictly increasing finite sample positions
        real(dp), intent(in) :: y(:) !! finite data ordinates, shape same as x
        integer, intent(in), optional :: k !! spline polynomial degree, default 3
        character(len=*), intent(in), optional :: bc_type !! not-a-knot (default), natural, clamped, periodic
        integer, intent(out), optional :: status !! success, invalid-input, or singular
        logical, intent(in), optional :: extrapolate !! extrapolate or wrap for periodic splines (default true)
        real(dp), intent(in), optional :: knots(:) !! explicit knot vector, if provided
        type(bspline_interpolator_type) :: s
        real(dp), allocatable :: t(:), mat(:,:), rhs(:,:), basis(:,:), row(:)
        integer :: deg, n, nc, code, i, j, span, nl
        character(len=:), allocatable :: bc
        code = bspline_status_invalid_input
        n = size(x)
        deg = 3
        if (present(k)) deg = k
        if (deg < 0 .or. deg > 15 .or. n < 2) goto 900
        if (size(y) /= n) goto 900
        if (.not. all(ieee_is_finite(x)) .or. .not. all(ieee_is_finite(y))) goto 900
        if (any(x(2:) <= x(:n-1))) goto 900
        bc = 'not-a-knot'
        if (present(bc_type)) bc = trim(bc_type)
        if (bc /= 'not-a-knot' .and. bc /= 'natural' .and. bc /= 'clamped' &
            .and. bc /= 'periodic') goto 900
        if (bc == 'periodic') then
            if (present(knots)) goto 900
            ! SciPy treats degrees zero and one as the ordinary interpolant,
            ! even when periodic boundary conditions are selected.
            if (deg <= 1) then
                if (abs(y(1)-y(n)) > 1.0e-12_dp*max(1.0_dp,maxval(abs(y)))) goto 900
                bc = 'not-a-knot'
            else
                call solve_periodic_collocation(x,reshape(y,[n,1]),deg,t,rhs,code)
                if (code /= bspline_status_success) goto 900
                if (present(extrapolate)) then
                    s = create_bspline(t,rhs(:,1),deg,status=code, &
                        extrapolate=.false.,periodic=extrapolate)
                else
                    s = create_bspline(t,rhs(:,1),deg,status=code, &
                        extrapolate=.false.,periodic=.true.)
                end if
                goto 900
            end if
        end if
        if (bc /= 'not-a-knot' .and. deg /= 3) goto 900
        if (deg <= 1) then
            if (bc /= 'not-a-knot' .or. present(knots)) goto 900
            if (deg == 0) then
                t = [x,x(n)]
            else
                t = [x(1),x,x(n)]
            end if
            s = create_bspline(t,y,deg,status=code,extrapolate=extrapolate)
            goto 900
        end if
        if (present(knots)) then
            t = knots
        else if (bc /= 'not-a-knot') then
            ! Augmented knots: k repeated outer positions in addition to x.
            t = [spread(x(1),1,deg),x,spread(x(n),1,deg)]
        else
            if (n < deg+1) goto 900
            ! SciPy's not-a-knot construction for odd/even degrees.
            allocate(row(n-1))
            row = 0.5_dp*(x(1:n-1)+x(2:n))
            if (mod(deg,2) == 1) then
                i = (deg+1)/2
                t = [spread(x(1),1,deg+1),x(i+1:n-i),spread(x(n),1,deg+1)]
            else
                i = deg/2
                t = [spread(x(1),1,deg+1),row(i+1:n-1-i),spread(x(n),1,deg+1)]
            end if
        end if
        if (size(t) < 2*deg+2) goto 900
        nc = size(t)-deg-1
        nl = 0
        if (bc /= 'not-a-knot') nl = 1
        if (nc /= n+2*nl) goto 900
        if (.not. all(ieee_is_finite(t))) goto 900
        if (any(t(2:) < t(:size(t)-1))) goto 900
        if (t(deg+1) > x(1) .or. t(nc+1) < x(n)) goto 900
        allocate(mat(nc,nc),rhs(nc,1))
        mat = 0.0_dp
        rhs = 0.0_dp
        do i=1,n
            span = find_span(t,deg,x(i))
            if (span > nc .or. span < deg+1) goto 900
            allocate(basis(0:deg,1))
            call basis_values(t,deg,span,x(i),basis(:,1))
            do j=0,deg
                mat(i+nl,span-deg+j) = basis(j,1)
            end do
            deallocate(basis)
        end do
        rhs(nl+1:nl+n,1) = y
        if (nl == 1) then
            i = 1
            if (bc == 'natural') i = 2
            call boundary_weights(t,deg,x(1),i,mat(1,:))
            call boundary_weights(t,deg,x(n),i,mat(nc,:))
        end if
        call solve_collocation(mat,rhs,code)
        if (code /= bspline_status_success) goto 900
        s = create_bspline(t,rhs(:,1),deg,status=code,extrapolate=extrapolate)
900     if (present(status)) status=code
    end function make_interp_spline_scalar

    function make_interp_spline_array(x, y, k, bc_type, status, extrapolate, knots) result(s)
        real(dp), intent(in) :: x(:) !! sample abscissas
        real(dp), intent(in) :: y(:,:) !! ordinates (sample, channel), SciPy axis=0
        integer, intent(in), optional :: k !! spline polynomial degree, default 3
        character(len=*), intent(in), optional :: bc_type !! not-a-knot, natural, clamped, periodic
        integer, intent(out), optional :: status !! construction status
        logical, intent(in), optional :: extrapolate !! extrapolate or wrap for periodic splines (default true)
        real(dp), intent(in), optional :: knots(:) !! optional knot vector
        type(bspline_array_interpolator_type) :: s
        real(dp), allocatable :: t(:), coeff(:,:)
        integer :: j, code, deg
        code = bspline_status_invalid_input
        if (size(y,1) /= size(x) .or. size(y,2) < 1) goto 900
        if (present(bc_type)) then
            if (trim(bc_type) == 'periodic') then
                if (present(knots)) goto 900
                deg = 3
                if (present(k)) deg = k
                if (deg <= 1) goto 100
                call solve_periodic_collocation(x,y,deg,t,coeff,code)
                if (code /= bspline_status_success) goto 900
                if (present(extrapolate)) then
                    s = create_bspline_array(t,coeff,deg,status=code, &
                        extrapolate=.false.,periodic=extrapolate)
                else
                    s = create_bspline_array(t,coeff,deg,status=code, &
                        extrapolate=.false.,periodic=.true.)
                end if
                goto 900
            end if
        end if
100     allocate(s%columns(size(y,2)))
        do j=1,size(y,2)
            s%columns(j) = make_interp_spline_scalar(x,y(:,j),k,bc_type,code,extrapolate,knots)
            if (code /= bspline_status_success) then
                deallocate(s%columns)
                goto 900
            end if
        end do
900     if (present(status)) status = code
    end function make_interp_spline_array


    ! SciPy's make_lsq_spline with an explicit complete knot vector.
    ! The compact streaming Givens QR algorithm stores O((k+1)*ncoeff)
    ! numbers, avoids the loss of accuracy from normal equations, and
    ! factors the weighted basis only once for all ordinate channels.
    ! Unlike LSQUnivariateSpline, this allows nondecreasing sample sites,
    ! zero observation weights, repeated internal knots, and degrees 0:15.
    subroutine fit_full_knots_lsq(x, y, t, k, w_in, c, status)
        real(dp), intent(in) :: x(:) !! finite nondecreasing observation positions (m)
        real(dp), intent(in) :: y(:,:) !! finite observations (m, number of channels)
        real(dp), intent(in) :: t(:) !! complete nondecreasing knot vector (n+k+1)
        integer, intent(in) :: k !! polynomial degree, 0:15
        real(dp), intent(in), optional :: w_in(:) !! nonnegative residual weights, length m; default 1
        real(dp), allocatable, intent(out) :: c(:,:) !! fitted coefficients (n, channels), only if successful
        integer, intent(out) :: status !! 0 success, 1 invalid input, 2 rank deficient
        real(dp), allocatable :: weights(:), upper(:,:), rhs(:,:), z(:), basis(:), row(:)
        real(dp) :: a, b, rho, scale, cosine, sine, old, tmp, tol
        integer :: m, n, nch, i, j, l, span, first, ch

        status = bspline_status_invalid_input
        m = size(x)
        nch = size(y,2)
        if (k < 0 .or. k > 15 .or. m < k+1 .or. nch < 1) return
        if (size(y,1) /= m .or. size(t) < 2*k+2) return
        n = size(t)-k-1
        if (n > m) return
        if (.not. all(ieee_is_finite(x)) .or. .not. all(ieee_is_finite(y)) &
            .or. .not. all(ieee_is_finite(t))) return
        if (any(x(2:) < x(:m-1))) return
        if (any(t(2:) < t(:size(t)-1))) return
        if (t(k+1) >= t(n+1)) return
        if (x(1) < t(k+1) .or. x(m) > t(n+1)) return
        if (present(w_in)) then
            if (size(w_in) /= m) return
            if (.not. all(ieee_is_finite(w_in))) return
            if (any(w_in < 0.0_dp)) return
            weights = w_in
        else
            allocate(weights(m))
            weights = 1.0_dp
        end if
        allocate(upper(0:k,n), rhs(n,nch), basis(0:k), row(0:k), z(nch))
        upper = 0.0_dp
        rhs = 0.0_dp
        do i = 1,m
            if (weights(i) == 0.0_dp) cycle
            span = find_span(t,k,x(i))
            first = span-k
            call basis_values(t,k,span,x(i),basis)
            row = basis*weights(i)
            z = y(i,:)*weights(i)
            ! Monotone x guarantees the support does not move left. At
            ! each step the nonzero row entries start at column first;
            ! earlier rows have already been eliminated by Givens rotations.
            do j = first,span
                a = upper(0,j)
                b = row(j-first)
                scale = max(abs(a),abs(b))
                if (scale == 0.0_dp) cycle
                rho = scale*sqrt((a/scale)**2+(b/scale)**2)
                cosine = a/rho
                sine = b/rho
                upper(0,j) = rho
                do l = j+1,span
                    old = upper(l-j,j)
                    tmp = row(l-first)
                    upper(l-j,j) = cosine*old+sine*tmp
                    row(l-first) = -sine*old+cosine*tmp
                end do
                do ch=1,nch
                    old = rhs(j,ch)
                    tmp = z(ch)
                    rhs(j,ch) = cosine*old+sine*tmp
                    z(ch) = -sine*old+cosine*tmp
                end do
            end do
        end do
        status = bspline_status_singular
        if (.not. all(ieee_is_finite(upper)) .or. .not. all(ieee_is_finite(rhs))) return
        tol = 64.0_dp*epsilon(1.0_dp)*maxval(abs(upper(0,:)))
        if (any(abs(upper(0,:)) <= tol)) return
        allocate(c(n,nch))
        c = rhs
        do j = n,1,-1
            do l = j+1,min(n,j+k)
                c(j,:) = c(j,:)-upper(l-j,j)*c(l,:)
            end do
            c(j,:) = c(j,:)/upper(0,j)
        end do
        if (.not. all(ieee_is_finite(c))) then
            deallocate(c)
            return
        end if
        status = bspline_status_success
    end subroutine fit_full_knots_lsq

    function make_lsq_spline_scalar(x, y, t, k, w, status, extrapolate, method) result(s)
        real(dp), intent(in) :: x(:) !! nondecreasing observation coordinates (m)
        real(dp), intent(in) :: y(:) !! scalar observations (m)
        real(dp), intent(in) :: t(:) !! *complete* B-spline knots, unlike LSQUnivariateSpline
        integer, intent(in), optional :: k !! degree 0:15, default 3
        real(dp), intent(in), optional :: w(:) !! nonnegative residual weights, default one
        integer, intent(out), optional :: status !! 0 success, 1 invalid input, 2 singular
        logical, intent(in), optional :: extrapolate !! polynomial extrapolation (default true)
        character(len=*), intent(in), optional :: method !! 'qr' only; other methods rejected
        type(bspline_interpolator_type) :: s
        real(dp), allocatable :: coeff(:,:), ordinates(:,:)
        integer :: deg, code

        code = bspline_status_invalid_input
        deg = 3
        if (present(k)) deg = k
        if (present(method)) then
            if (trim(method) /= 'qr') goto 900
        end if
        ordinates = reshape(y,[size(y),1])
        call fit_full_knots_lsq(x,ordinates,t,deg,w,coeff,code)
        if (code == bspline_status_success) then
            s = create_bspline(t,coeff(:,1),deg,status=code,extrapolate=extrapolate)
        end if
900     if (present(status)) status = code
    end function make_lsq_spline_scalar

    function make_lsq_spline_channels(x, y, t, k, w, status, extrapolate, method) result(s)
        real(dp), intent(in) :: x(:) !! nondecreasing sample coordinates (m)
        real(dp), intent(in) :: y(:,:) !! observations (m, channels)
        real(dp), intent(in) :: t(:) !! complete knot vector (n+k+1)
        integer, intent(in), optional :: k !! degree 0:15, default 3
        real(dp), intent(in), optional :: w(:) !! nonnegative residual weights, default one
        integer, intent(out), optional :: status !! 0 success, 1 invalid input, 2 singular
        logical, intent(in), optional :: extrapolate !! polynomial extrapolation (default true)
        character(len=*), intent(in), optional :: method !! 'qr' only
        type(bspline_array_interpolator_type) :: s
        real(dp), allocatable :: coeff(:,:)
        integer :: deg, code

        code = bspline_status_invalid_input
        deg = 3
        if (present(k)) deg = k
        if (present(method)) then
            if (trim(method) /= 'qr') goto 900
        end if
        call fit_full_knots_lsq(x,y,t,deg,w,coeff,code)
        if (code == bspline_status_success) then
            s = create_bspline_array(t,coeff,deg,status=code,extrapolate=extrapolate)
        end if
900     if (present(status)) status = code
    end function make_lsq_spline_channels

    ! FITPACK-style least-squares B-spline with user-supplied interior knots.
    ! SciPy's LSQUnivariateSpline minimizes sum((w_i * (y_i - f(x_i)))**2).
    ! This entry point preserves its historical degree 1:5, strict x,
    ! strictly increasing interior knots and positive-weight constraints;
    ! the shared banded QR solver also supports the more general factory.
    subroutine fit_lsq_spline(x, y, interior, w_in, k_in, bbox, t, c, status)
        real(dp), intent(in) :: x(:) !! strictly increasing sample sites, length m
        real(dp), intent(in) :: y(:,:) !! observations, shape (m, channels)
        real(dp), intent(in) :: interior(:) !! strictly increasing internal knots
        real(dp), intent(in), optional :: w_in(:) !! positive residual weights
        integer, intent(in), optional :: k_in !! degree, default 3, supported 1:5
        real(dp), intent(in), optional :: bbox(:) !! support bounds [left,right]
        real(dp), allocatable, intent(out) :: t(:) !! complete knot vector
        real(dp), allocatable, intent(out) :: c(:,:) !! fitted coefficients
        integer, intent(out) :: status !! success, invalid input or singular
        real(dp) :: left, right
        integer :: m, deg

        status = bspline_status_invalid_input
        m = size(x)
        deg = 3
        if (present(k_in)) deg = k_in
        if (deg < 1 .or. deg > 5 .or. m < deg+1 .or. size(y,1) /= m &
            .or. size(y,2) < 1 .or. size(interior) > m-deg-1) return
        if (.not. all(ieee_is_finite(x)) .or. .not. all(ieee_is_finite(y)) &
            .or. .not. all(ieee_is_finite(interior))) return
        if (any(x(2:) <= x(:m-1))) return
        if (size(interior) > 1) then
            if (any(interior(2:) <= interior(:size(interior)-1))) return
        end if
        left = x(1)
        right = x(m)
        if (present(bbox)) then
            if (size(bbox) /= 2) return
            if (.not. all(ieee_is_finite(bbox))) return
            left = bbox(1)
            right = bbox(2)
        end if
        if (.not. (left < right .and. left <= x(1) .and. right >= x(m))) return
        if (any(interior <= left) .or. any(interior >= right)) return
        if (present(w_in)) then
            if (size(w_in) /= m) return
            if (.not. all(ieee_is_finite(w_in))) return
            if (any(w_in <= 0.0_dp)) return
        end if
        t = [spread(left,1,deg+1),interior,spread(right,1,deg+1)]
        call fit_full_knots_lsq(x,y,t,deg,w_in,c,status)
    end subroutine fit_lsq_spline

    function lsq_spline_scalar(x, y, t, w, k, bbox, status, extrapolate) result(s)
        real(dp), intent(in) :: x(:) !! sample coordinates, strictly increasing
        real(dp), intent(in) :: y(:) !! scalar data values
        real(dp), intent(in) :: t(:) !! interior knots (not repeated end knots)
        real(dp), intent(in), optional :: w(:) !! positive residual weights
        integer, intent(in), optional :: k !! degree 1:5, default 3
        real(dp), intent(in), optional :: bbox(:) !! support interval [left,right]
        integer, intent(out), optional :: status !! result code 0/1/2
        logical, intent(in), optional :: extrapolate !! allow extrapolation (default true)
        type(bspline_interpolator_type) :: s
        real(dp), allocatable :: knots(:), coef(:,:), data(:,:)
        integer :: deg, code
        data = reshape(y, [size(y),1])
        call fit_lsq_spline(x,data,t,w,k,bbox,knots,coef,code)
        if (code == bspline_status_success) then
            deg = 3
            if (present(k)) deg = k
            s = create_bspline(knots,coef(:,1),deg,status=code,extrapolate=extrapolate)
        end if
        if (present(status)) status = code
    end function lsq_spline_scalar

    function lsq_spline_matrix(x, y, t, w, k, bbox, status, extrapolate) result(s)
        real(dp), intent(in) :: x(:) !! sample coordinates, strictly increasing
        real(dp), intent(in) :: y(:,:) !! output channels, shape (n,channels)
        real(dp), intent(in) :: t(:) !! interior knots
        real(dp), intent(in), optional :: w(:) !! positive residual weights
        integer, intent(in), optional :: k !! spline degree 1:5, default 3
        real(dp), intent(in), optional :: bbox(:) !! domain bounds [left,right]
        integer, intent(out), optional :: status !! result code 0/1/2
        logical, intent(in), optional :: extrapolate !! allow extrapolation (default true)
        type(bspline_array_interpolator_type) :: s
        real(dp), allocatable :: knots(:), coef(:,:)
        integer :: deg, code
        call fit_lsq_spline(x,y,t,w,k,bbox,knots,coef,code)
        if (code == bspline_status_success) then
            deg = 3
            if (present(k)) deg = k
            s = create_bspline_array(knots,coef,deg,status=code,extrapolate=extrapolate)
        end if
        if (present(status)) status = code
    end function lsq_spline_matrix

    ! Exact-data counterpart to LSQUnivariateSpline. For s=0, the FITPACK
    ! interpolating spline has the not-a-knot knot placement used here.
    function interpolated_spline_scalar(x, y, k, status, extrapolate) result(s)
        real(dp), intent(in) :: x(:) !! strictly increasing observations
        real(dp), intent(in) :: y(:) !! scalar ordinates
        integer, intent(in), optional :: k !! polynomial degree, default 3
        integer, intent(out), optional :: status !! construction status
        logical, intent(in), optional :: extrapolate !! endpoint extrapolation
        type(bspline_interpolator_type) :: s
        if (present(k)) then
            if (k < 1 .or. k > 5) then
                if (present(status)) status=bspline_status_invalid_input
                return
            end if
        end if
        s = make_interp_spline_scalar(x,y,k,status=status,extrapolate=extrapolate)
    end function interpolated_spline_scalar

    function interpolated_spline_matrix(x, y, k, status, extrapolate) result(s)
        real(dp), intent(in) :: x(:) !! strictly increasing observations
        real(dp), intent(in) :: y(:,:) !! observations (sample,channel)
        integer, intent(in), optional :: k !! degree, default 3
        integer, intent(out), optional :: status !! construction status
        logical, intent(in), optional :: extrapolate !! endpoint extrapolation
        type(bspline_array_interpolator_type) :: s
        if (present(k)) then
            if (k < 1 .or. k > 5) then
                if (present(status)) status=bspline_status_invalid_input
                return
            end if
        end if
        s = make_interp_spline_array(x,y,k,status=status,extrapolate=extrapolate)
    end function interpolated_spline_matrix

end module scifort_interpolate_bspline
