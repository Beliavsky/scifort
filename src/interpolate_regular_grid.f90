! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Rectilinear nearest, multilinear, monotone PCHIP, and tensor-product B-splines.
! Independent implementation of standard interpolation and tensor-product
! collocation. Reuses SciFort's own B-spline basis matrix; no upstream source
! has been imported or copied.
module scifort_interpolate_regular_grid
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan, ieee_value, ieee_quiet_nan
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_kinds, only : dp
    implicit none
    private

    integer, parameter :: max_grid_dim = 6
    integer, parameter, public :: grid_status_success = 0
    integer, parameter, public :: grid_status_invalid_input = 1
    integer, parameter, public :: grid_status_out_of_bounds = 2
    integer, parameter, public :: grid_status_singular = 3

    type, public :: grid_axis
        real(dp), allocatable :: points(:) !! One strictly monotonic, finite coordinate vector
    end type grid_axis

    type, public :: regular_grid_interpolator_type
        private
        type(grid_axis), allocatable :: axes(:)
        integer, allocatable :: lengths(:)
        integer, allocatable :: strides(:)
        real(dp), allocatable :: values(:,:)  ! (grid cells in Fortran order, channels)
        type(grid_axis), allocatable :: spline_knots(:) !! not-a-knot B-spline knot vectors
        real(dp), allocatable :: spline_coefficients(:,:) !! fitted coefficients, cell x channel
        integer :: spline_degree = -1 !! degree 1, 3 or 5; -1 for nearest/linear/PCHIP
        character(len=7) :: method = 'linear'
        logical :: bounds_error = .true.
        logical :: use_fill = .true.
        real(dp) :: fill_value = 0.0_dp
    contains
        procedure :: valid => valid_grid
        procedure :: dimensions => grid_dimensions
        procedure :: channels => grid_channels
        procedure, private :: evaluate_one
        procedure, private :: evaluate_many
        generic :: evaluate => evaluate_one, evaluate_many
        procedure :: evaluate_channels => evaluate_all_channels
        procedure :: evaluate_values => evaluate_scalar_batch
    end type regular_grid_interpolator_type

    type, public :: complex_regular_grid_interpolator_type
        private
        type(regular_grid_interpolator_type) :: real_part
        type(regular_grid_interpolator_type) :: imag_part
    contains
        procedure :: valid => valid_complex_grid
        procedure :: dimensions => complex_grid_dimensions
        procedure :: channels => complex_grid_channels
        procedure, private :: evaluate_complex_one
        procedure, private :: evaluate_complex_many
        generic :: evaluate => evaluate_complex_one, evaluate_complex_many
        procedure :: evaluate_channels => evaluate_complex_channels
        procedure :: evaluate_values => evaluate_complex_values
    end type complex_regular_grid_interpolator_type

    public :: regular_grid_interpolator, complex_regular_grid_interpolator
    public :: interpn, complex_interpn

    interface interpn
        module procedure interpn_real
    end interface
    interface complex_interpn
        module procedure interpn_complex
    end interface

contains

    function regular_grid_interpolator(axes, values, method, bounds_error, fill_value, extrapolate, status) result(grid)
        type(grid_axis), intent(in) :: axes(:) !! N strictly increasing or decreasing coordinate axes
        real(dp), intent(in) :: values(..) !! rank N grid data, or rank N+1 with trailing channels
        character(len=*), intent(in), optional :: method !! 'linear', 'nearest', 'pchip', 'slinear', 'cubic', 'quintic'
        logical, intent(in), optional :: bounds_error !! flag out-of-range queries (default true)
        real(dp), intent(in), optional :: fill_value !! out-of-domain value (default quiet NaN)
        logical, intent(in), optional :: extrapolate !! if true, use edge cells outside instead of a fill value
        integer, intent(out), optional :: status !! zero on success, invalid input otherwise
        type(regular_grid_interpolator_type) :: grid
        real(dp), allocatable :: flat(:)
        integer :: npoint, nchannel, code
        integer, allocatable :: shape_values(:)

        code = grid_status_invalid_input
        if (.not. valid_axes(axes)) goto 900
        npoint = grid_cell_count(axes)
        shape_values = shape(values)
        if (rank(values) == size(axes)) then
            if (any(shape(values) /= axis_sizes(axes))) goto 900
            nchannel = 1
        else if (rank(values) == size(axes)+1) then
            if (any(shape_values(:size(axes)) /= axis_sizes(axes))) goto 900
            nchannel = size(values)/npoint
            if (nchannel < 1) goto 900
        else
            goto 900
        end if
        select rank (values)
        rank (1)
            flat = values
        rank (2)
            flat = reshape(values,[size(values)])
        rank (3)
            flat = reshape(values,[size(values)])
        rank (4)
            flat = reshape(values,[size(values)])
        rank (5)
            flat = reshape(values,[size(values)])
        rank (6)
            flat = reshape(values,[size(values)])
        rank (7)
            flat = reshape(values,[size(values)])
        rank default
            goto 900
        end select
        call initialize_grid(grid, axes, reshape(flat,[npoint,nchannel]), method, &
            bounds_error, fill_value, extrapolate, code)
900     if (present(status)) status = code
    end function regular_grid_interpolator

    function complex_regular_grid_interpolator(axes, values, method, bounds_error, fill_value, extrapolate, &
        status) result(grid)
        type(grid_axis), intent(in) :: axes(:) !! N strictly monotone coordinate vectors
        complex(dp), intent(in) :: values(..) !! complex rank N data or rank N+1 channels
        character(len=*), intent(in), optional :: method !! 'linear', 'nearest' or B-spline; PCHIP rejects complex data
        logical, intent(in), optional :: bounds_error !! whether out-of-range queries report failure
        complex(dp), intent(in), optional :: fill_value !! complex out-of-domain fill value
        logical, intent(in), optional :: extrapolate !! edge-cell extrapolation if true
        integer, intent(out), optional :: status !! success or invalid-input status
        type(complex_regular_grid_interpolator_type) :: grid
        complex(dp), allocatable :: flat(:)
        real(dp), allocatable :: yr(:,:), yi(:,:)
        real(dp) :: fillr, filli
        integer :: npoint, nchannel, code, code_im
        integer, allocatable :: shape_values(:)

        code = grid_status_invalid_input
        if (.not. valid_axes(axes)) goto 900
        ! SciPy >= 1.15 rejects complex-valued PCHIP interpolation.
        if (present(method)) then
            if (method == 'pchip') goto 900
        end if
        npoint = grid_cell_count(axes)
        shape_values = shape(values)
        if (rank(values) == size(axes)) then
            if (any(shape(values) /= axis_sizes(axes))) goto 900
            nchannel = 1
        else if (rank(values) == size(axes)+1) then
            if (any(shape_values(:size(axes)) /= axis_sizes(axes))) goto 900
            nchannel = size(values)/npoint
            if (nchannel < 1) goto 900
        else
            goto 900
        end if
        select rank (values)
        rank (1)
            flat = values
        rank (2)
            flat = reshape(values,[size(values)])
        rank (3)
            flat = reshape(values,[size(values)])
        rank (4)
            flat = reshape(values,[size(values)])
        rank (5)
            flat = reshape(values,[size(values)])
        rank (6)
            flat = reshape(values,[size(values)])
        rank (7)
            flat = reshape(values,[size(values)])
        rank default
            goto 900
        end select
        yr = reshape(real(flat, dp),[npoint,nchannel])
        yi = reshape(aimag(flat),[npoint,nchannel])
        fillr = ieee_value(0.0_dp,ieee_quiet_nan)
        filli = 0.0_dp
        if (present(fill_value)) then
            fillr = real(fill_value,dp)
            filli = aimag(fill_value)
        end if
        call initialize_grid(grid%real_part, axes, yr, method, bounds_error, fillr, extrapolate, code)
        if (code /= grid_status_success) goto 900
        call initialize_grid(grid%imag_part, axes, yi, method, bounds_error, filli, extrapolate, code_im)
        code = code_im
900     if (present(status)) status = code
    end function complex_regular_grid_interpolator

    pure logical function is_spline_method(method) result(ok)
        character(len=*), intent(in) :: method !! interpolation method name
        ok = method == 'slinear' .or. method == 'cubic' .or. method == 'quintic'
    end function is_spline_method

    pure logical function supported_grid_method(method) result(ok)
        character(len=*), intent(in) :: method !! interpolation method name
        ok = method == 'nearest' .or. method == 'linear' .or. method == 'pchip' .or. is_spline_method(method)
    end function supported_grid_method

    ! Fit exact interpolating B-splines along each axis. Collocation has
    ! half-bandwidth <= degree, so the specialized no-pivot banded LU uses
    ! O(n*degree) storage and O(n*degree**2) work, instead of the previous
    ! O(n**2) dense matrix and O(n**3) LU factorization.
    ! Collocation of B-spline bases at Schoenberg-Whitney points is totally
    ! nonnegative. For valid strictly interior knots its leading pivots are
    ! positive; a small-pivot check catches ill-conditioned input rather
    ! than quietly producing an inaccurate answer.
    subroutine initialize_tensor_spline(grid, values, status)
        type(regular_grid_interpolator_type), intent(inout) :: grid !! validated grid and axis strides
        real(dp), intent(in) :: values(:,:) !! original grid samples, cells x channels
        integer, intent(out) :: status !! success or singular/invalid input
        real(dp), allocatable :: coords(:), knots(:), band(:,:), rhs(:)
        real(dp) :: basis(0:5)
        integer :: deg, n, d, mid, ncell, stride, first, rest, ch, j, i, col, code
        logical :: descend

        status = grid_status_invalid_input
        deg = 1
        if (grid%method == 'cubic') deg = 3
        if (grid%method == 'quintic') deg = 5
        if (any(grid%lengths < deg+1)) return
        if (.not. all(ieee_is_finite(values))) return
        grid%spline_degree = deg
        allocate(grid%spline_knots(size(grid%axes)))
        grid%spline_coefficients = values
        ncell = size(values,1)
        do d = 1, size(grid%axes)
            n = grid%lengths(d)
            descend = grid%axes(d)%points(n) < grid%axes(d)%points(1)
            coords = grid%axes(d)%points
            if (descend) coords = -coords
            mid = (deg+1)/2
            knots = [spread(coords(1),1,deg+1), coords(mid+1:n-mid), &
                spread(coords(n),1,deg+1)]
            grid%spline_knots(d)%points = knots
            allocate(band(-deg:deg,1:n), rhs(n))
            band = 0.0_dp
            do i = 1,n
                call grid_basis(knots,deg,coords(i),first,basis)
                do j = 0,deg
                    col = first+j
                    if (col >= 1 .and. col <= n) then
                        if (abs(i-col) <= deg) band(i-col,col) = basis(j)
                    end if
                end do
            end do
            call factor_grid_banded(band,deg,code)
            if (code /= grid_status_success) then
                status = code
                return
            end if
            stride = grid%strides(d)
            do ch = 1, size(values,2)
                do first = 1, ncell, n*stride
                    do rest = 0, stride-1
                        do j = 1, n
                            rhs(j) = grid%spline_coefficients(first+rest+(j-1)*stride,ch)
                        end do
                        call solve_grid_banded(band,deg,rhs)
                        do j = 1, n
                            grid%spline_coefficients(first+rest+(j-1)*stride,ch) = rhs(j)
                        end do
                    end do
                end do
            end do
            deallocate(band, rhs)
        end do
        status = grid_status_success
    end subroutine initialize_tensor_spline

    subroutine factor_grid_banded(a, degree, status)
        integer, intent(in) :: degree !! half-bandwidth of the collocation matrix, <= 5
        real(dp), intent(inout) :: a(-degree:,:) !! a(i-j,j), overwritten by band LU factors
        integer, intent(out) :: status !! zero on success, singular otherwise
        real(dp) :: threshold, factor, pivot
        integer :: n, i, row, col

        n = size(a,2)
        status = grid_status_singular
        threshold = 32.0_dp*epsilon(1.0_dp)*maxval(abs(a))
        if (threshold <= 0.0_dp) return
        do i = 1,n
            pivot = a(0,i)
            if (.not. ieee_is_finite(pivot) .or. abs(pivot) <= threshold) return
            do row = i+1,min(n,i+degree)
                factor = a(row-i,i)/pivot
                a(row-i,i) = factor
                do col = i+1,min(n,i+degree)
                    if (abs(row-col) <= degree) &
                        a(row-col,col) = a(row-col,col)-factor*a(i-col,col)
                end do
            end do
        end do
        status = grid_status_success
    end subroutine factor_grid_banded

    subroutine solve_grid_banded(a, degree, rhs)
        integer, intent(in) :: degree !! half-bandwidth
        real(dp), intent(in) :: a(-degree:,:) !! factored band matrix
        real(dp), intent(inout) :: rhs(:) !! RHS overwritten by solution
        integer :: i, j, n

        n = size(rhs)
        do i = 2,n
            do j = max(1,i-degree),i-1
                rhs(i) = rhs(i)-a(i-j,j)*rhs(j)
            end do
        end do
        do i = n,1,-1
            do j = i+1,min(n,i+degree)
                rhs(i) = rhs(i)-a(i-j,j)*rhs(j)
            end do
            rhs(i) = rhs(i)/a(0,i)
        end do
    end subroutine solve_grid_banded

    ! SciPy's recursive PCHIP evaluator applies 1D PCHIP from the LAST
    ! interpolation axis to the FIRST. PCHIP is nonlinear in the ordinates;
    ! reversing this order changes the answer. The current stage is laid out
    ! in Fortran order, with the axis being eliminated as its slowest index.
    subroutine evaluate_tensor_pchip(grid, x, y)
        class(regular_grid_interpolator_type), intent(in) :: grid !! initialized PCHIP plan
        real(dp), intent(in) :: x(:) !! query coordinates, one per axis
        real(dp), intent(out) :: y(:) !! computed output channels
        real(dp), allocatable :: stage(:,:), folded(:,:), knots(:), samples(:)
        real(dp) :: query
        integer :: d, stride, n, ch, i, j

        stage = grid%values
        do d = size(grid%axes),1,-1
            stride = product(grid%lengths(:d-1))
            n = grid%lengths(d)
            allocate(folded(stride,size(y)), samples(n))
            knots = grid%axes(d)%points
            query = x(d)
            if (knots(n) < knots(1)) then
                knots = -knots
                query = -query
            end if
            do ch = 1,size(y)
                do i = 1,stride
                    do j = 1,n
                        samples(j) = stage(i+(j-1)*stride,ch)
                    end do
                    folded(i,ch) = evaluate_pchip_fiber(knots,samples,query)
                end do
            end do
            call move_alloc(folded,stage)
            deallocate(samples)
        end do
        y = stage(1,:)
    end subroutine evaluate_tensor_pchip

    pure real(dp) function pchip_endpoint(h1,h2,s1,s2) result(d)
        real(dp), intent(in) :: h1 !! distance of adjacent endpoint interval
        real(dp), intent(in) :: h2 !! distance of following interval
        real(dp), intent(in) :: s1 !! endpoint secant slope
        real(dp), intent(in) :: s2 !! next secant slope

        d = ((2.0_dp*h1+h2)*s1-h1*s2)/(h1+h2)
        if (s1 == 0.0_dp .or. ((d > 0.0_dp) .neqv. (s1 > 0.0_dp))) then
            d = 0.0_dp
        else if ((s1 > 0.0_dp) .neqv. (s2 > 0.0_dp)) then
            if (abs(d) > 3.0_dp*abs(s1)) d = 3.0_dp*s1
        end if
    end function pchip_endpoint

    ! Fritsch-Butland harmonic interior slopes and shape-preserving Hermite
    ! evaluation. The same formulas are used by SciFort's 1D PCHIP code.
    pure real(dp) function evaluate_pchip_fiber(knots, samples, query) result(val)
        real(dp), intent(in) :: knots(:) !! finite, ascending coordinates; length >= 4
        real(dp), intent(in) :: samples(:) !! finite ordinates, same size
        real(dp), intent(in) :: query !! interpolation or extrapolation point
        real(dp) :: h(size(knots)-1), delta(size(knots)-1), slope(size(knots))
        real(dp) :: dx, hi, d, u, s, weight1, weight2
        integer :: n, i, seg, lo, up, mid

        n = size(knots)
        h = knots(2:n)-knots(1:n-1)
        delta = (samples(2:n)-samples(1:n-1))/h
        slope(1) = pchip_endpoint(h(1),h(2),delta(1),delta(2))
        slope(n) = pchip_endpoint(h(n-1),h(n-2),delta(n-1),delta(n-2))
        do i = 2,n-1
            if (delta(i-1) == 0.0_dp .or. delta(i) == 0.0_dp) then
                slope(i) = 0.0_dp
            else if ((delta(i-1) > 0.0_dp) .neqv. (delta(i) > 0.0_dp)) then
                slope(i) = 0.0_dp
            else
                weight1 = (2.0_dp*h(i)+h(i-1))/(3.0_dp*(h(i-1)+h(i)))
                weight2 = (h(i)+2.0_dp*h(i-1))/(3.0_dp*(h(i-1)+h(i)))
                slope(i) = 1.0_dp/(weight1/delta(i-1)+weight2/delta(i))
            end if
        end do
        if (query <= knots(1)) then
            seg = 1
        else if (query >= knots(n)) then
            seg = n-1
        else
            lo = 1
            up = n
            do while (up-lo > 1)
                mid = (up+lo)/2
                if (query >= knots(mid)) then
                    lo = mid
                else
                    up = mid
                end if
            end do
            seg = lo
        end if
        dx = query-knots(seg)
        hi = h(seg)
        d = delta(seg)
        ! c2 and c3 are local Hermite coefficients; Horner minimizes work.
        u = (3.0_dp*d-2.0_dp*slope(seg)-slope(seg+1))/hi
        s = (slope(seg)+slope(seg+1)-2.0_dp*d)/(hi*hi)
        val = samples(seg)+dx*(slope(seg)+dx*(u+dx*s))
    end function evaluate_pchip_fiber

    ! Local B-spline support and basis values (Cox-de Boor recurrence).
    ! In particular, the endpoint span also applies to polynomial
    ! extrapolation outside the grid without a separate special-case fit.
    subroutine grid_basis(t, degree, x, first, basis, derivative_order, left_limit)
        real(dp), intent(in) :: t(:) !! nondecreasing knot vector
        integer, intent(in) :: degree !! spline polynomial degree, in {1,3,5}
        real(dp), intent(in) :: x !! finite coordinate in transformed axis
        integer, intent(out) :: first !! first nonzero coefficient index
        real(dp), intent(out) :: basis(0:5) !! local basis weights or derivatives, 0:degree valid
        integer, intent(in), optional :: derivative_order !! nonnegative derivative order
        logical, intent(in), optional :: left_limit !! at interior knots select left-side polynomial
        real(dp) :: left(5), right(5), saved, temp, denominator
        real(dp) :: levels(0:5,0:6,0:5), a, b
        integer :: derivative, q, i, index
        integer :: n, span, lo, hi, mid, j, r
        logical :: take_left

        take_left = .false.
        if (present(left_limit)) take_left = left_limit
        n = size(t)-degree-1
        if (x <= t(degree+1)) then
            span = degree+1
        else if (x >= t(n+1)) then
            span = n
        else
            lo = degree+1
            hi = n+1
            do while (hi-lo > 1)
                mid = (lo+hi)/2
                if (x < t(mid) .or. (take_left .and. x == t(mid))) then
                    hi = mid
                else
                    lo = mid
                end if
            end do
            span = lo
        end if
        first = span-degree
        basis = 0.0_dp
        derivative = 0
        if (present(derivative_order)) derivative = derivative_order
        if (derivative > degree) return
        if (derivative > 0) then
            levels = 0.0_dp
            levels(0,degree,0) = 1.0_dp
            do q = 1, degree
                do i = span-q, span
                    index = i-first
                    a = 0.0_dp
                    b = 0.0_dp
                    denominator = t(i+q)-t(i)
                    if (denominator /= 0.0_dp) a = (x-t(i))/denominator
                    denominator = t(i+q+1)-t(i+1)
                    if (denominator /= 0.0_dp) b = (t(i+q+1)-x)/denominator
                    if (a /= 0.0_dp) levels(q,index,0) = levels(q,index,0) + a*levels(q-1,index,0)
                    if (b /= 0.0_dp) levels(q,index,0) = levels(q,index,0) + b*levels(q-1,index+1,0)
                end do
            end do
            do r = 1,derivative
                do q = r,degree
                    do i = span-q,span
                        index = i-first
                        a = 0.0_dp
                        b = 0.0_dp
                        denominator = t(i+q)-t(i)
                        if (denominator /= 0.0_dp) a = real(q,dp)/denominator
                        denominator = t(i+q+1)-t(i+1)
                        if (denominator /= 0.0_dp) b = real(q,dp)/denominator
                        levels(q,index,r) = a*levels(q-1,index,r-1) - b*levels(q-1,index+1,r-1)
                    end do
                end do
            end do
            basis(0:degree) = levels(degree,0:degree,derivative)
            return
        end if
        basis(0) = 1.0_dp
        do j = 1, degree
            left(j) = x-t(span+1-j)
            right(j) = t(span+j)-x
            saved = 0.0_dp
            do r = 0,j-1
                denominator = right(r+1)+left(j-r)
                temp = 0.0_dp
                if (denominator /= 0.0_dp) temp = basis(r)/denominator
                basis(r) = saved+right(r+1)*temp
                saved = left(j-r)*temp
            end do
            basis(j) = saved
        end do
    end subroutine grid_basis

    subroutine evaluate_tensor_spline(grid, x, y, nu)
        class(regular_grid_interpolator_type), intent(in) :: grid !! initialized spline plan
        real(dp), intent(in) :: x(:) !! one physical point, ndim coordinates
        real(dp), intent(out) :: y(:) !! interpolated output channels
        integer, intent(in), optional :: nu(:) !! nonnegative per-axis derivative orders
        real(dp) :: basis(0:5,max_grid_dim), coord, weight
        integer :: first(max_grid_dim), deg, nd, k, j, digit, offset, index, combos

        deg = grid%spline_degree
        nd = size(grid%axes)
        do k = 1, nd
            coord = x(k)
            if (grid%axes(k)%points(grid%lengths(k)) < grid%axes(k)%points(1)) coord = -coord
            if (present(nu)) then
                call grid_basis(grid%spline_knots(k)%points,deg,coord,first(k),basis(:,k),nu(k), &
                    grid%axes(k)%points(grid%lengths(k)) < grid%axes(k)%points(1))
            else
                call grid_basis(grid%spline_knots(k)%points,deg,coord,first(k),basis(:,k))
            end if
            ! The transformed descending axis uses q=-x. The chain rule
            ! contributes a factor (-1)**nu for odd derivatives.
            if (present(nu)) then
                if (grid%axes(k)%points(grid%lengths(k)) < grid%axes(k)%points(1)) then
                    if (mod(nu(k),2) == 1) basis(:,k) = -basis(:,k)
                end if
            end if
        end do
        y = 0.0_dp
        combos = (deg+1)**nd
        do index = 0,combos-1
            offset = 1
            weight = 1.0_dp
            j = index
            do k = 1,nd
                digit = mod(j,deg+1)
                j = j/(deg+1)
                weight = weight*basis(digit,k)
                offset = offset+(first(k)+digit-1)*grid%strides(k)
            end do
            y = y+weight*grid%spline_coefficients(offset,:)
        end do
    end subroutine evaluate_tensor_spline

    subroutine initialize_grid(grid, axes, values, method, bounds_error, fill_value, extrapolate, status)
        type(regular_grid_interpolator_type), intent(out) :: grid !! initialized output grid
        type(grid_axis), intent(in) :: axes(:) !! validated monotone grid axes
        real(dp), intent(in) :: values(:,:) !! gridpoint-major columns, Fortran order
        character(len=*), intent(in), optional :: method !! nearest, linear, pchip, slinear, cubic or quintic
        logical, intent(in), optional :: bounds_error !! whether to reject out-of-domain queries
        real(dp), intent(in), optional :: fill_value !! value for out-of-bounds queries
        logical, intent(in), optional :: extrapolate !! whether to extrapolate
        integer, intent(out) :: status !! zero for success
        integer :: k

        status = grid_status_invalid_input
        if (.not. valid_axes(axes)) return
        if (size(values,1) /= grid_cell_count(axes) .or. size(values,2) < 1) return
        if (present(method)) then
            if (.not. supported_grid_method(method)) return
            grid%method = method
        end if
        grid%bounds_error = .true.
        if (present(bounds_error)) grid%bounds_error = bounds_error
        grid%fill_value = ieee_value(0.0_dp,ieee_quiet_nan)
        if (present(fill_value)) grid%fill_value = fill_value
        grid%use_fill = .true.
        if (present(extrapolate)) grid%use_fill = .not. extrapolate
        grid%axes = axes
        grid%lengths = axis_sizes(axes)
        allocate(grid%strides(size(axes)))
        grid%strides(1) = 1
        do k = 2, size(axes)
            grid%strides(k) = grid%strides(k-1)*grid%lengths(k-1)
        end do
        if (grid%method == 'pchip') then
            if (any(grid%lengths < 4)) return
            if (.not. all(ieee_is_finite(values))) return
        end if
        if (is_spline_method(grid%method)) then
            call initialize_tensor_spline(grid, values, status)
            if (status /= grid_status_success) return
        end if
        grid%values = values
        status = grid_status_success
    end subroutine initialize_grid

    logical function valid_axes(axes) result(ok)
        type(grid_axis), intent(in) :: axes(:) !! axes to validate
        integer :: i, n
        real(dp), allocatable :: deltas(:)
        integer(int64) :: cells

        ok = .false.
        if (size(axes) < 1 .or. size(axes) > max_grid_dim) return
        cells = 1_int64
        do i = 1, size(axes)
            if (.not. allocated(axes(i)%points)) return
            n = size(axes(i)%points)
            if (n < 1) return
            if (.not. all(ieee_is_finite(axes(i)%points))) return
            if (n > 1) then
                deltas = axes(i)%points(2:n)-axes(i)%points(1:n-1)
                if (.not. all(ieee_is_finite(deltas))) return
                if (.not. (all(deltas > 0.0_dp) .or. all(deltas < 0.0_dp))) return
            end if
            if (cells > int(huge(0),int64)/int(n,int64)) return
            cells = cells*int(n,int64)
        end do
        ok = .true.
    end function valid_axes

    function axis_sizes(axes) result(lengths)
        type(grid_axis), intent(in) :: axes(:) !! initialized axis array
        integer :: lengths(size(axes))
        integer :: k
        do k = 1, size(axes)
            lengths(k) = size(axes(k)%points)
        end do
    end function axis_sizes

    integer function grid_cell_count(axes) result(count)
        type(grid_axis), intent(in) :: axes(:) !! validated axis collection
        count = product(axis_sizes(axes))
    end function grid_cell_count

    logical function valid_grid(self)
        class(regular_grid_interpolator_type), intent(in) :: self !! grid to test
        valid_grid = allocated(self%values) .and. allocated(self%axes)
    end function valid_grid

    integer function grid_dimensions(self)
        class(regular_grid_interpolator_type), intent(in) :: self !! interpolation object
        grid_dimensions = 0
        if (self%valid()) grid_dimensions = size(self%axes)
    end function grid_dimensions

    integer function grid_channels(self)
        class(regular_grid_interpolator_type), intent(in) :: self !! interpolation object
        grid_channels = 0
        if (self%valid()) grid_channels = size(self%values,2)
    end function grid_channels

    function evaluate_one(self, x, status, method, nu) result(y)
        class(regular_grid_interpolator_type), intent(in) :: self !! configured grid
        real(dp), intent(in) :: x(:) !! one query of length N
        integer, intent(out), optional :: status !! success, invalid, or out of bounds
        character(len=*), intent(in), optional :: method !! override nearest, linear or B-spline
        integer, intent(in), optional :: nu(:) !! partial derivative order along each grid axis
        real(dp) :: y
        real(dp), allocatable :: channels(:)
        integer :: code
        y = ieee_value(0.0_dp,ieee_quiet_nan)
        code = grid_status_invalid_input
        if (self%channels() /= 1) goto 900
        channels = self%evaluate_channels(x, status=code, method=method, nu=nu)
        if (code == grid_status_success) y = channels(1)
900     if (present(status)) status = code
    end function evaluate_one

    recursive function evaluate_many(self, x, status, method, nu) result(y)
        class(regular_grid_interpolator_type), intent(in) :: self !! configured grid
        real(dp), intent(in) :: x(:,:) !! row-oriented queries, shape (nquery,ndim)
        integer, intent(out), optional :: status !! first failure code, or zero
        character(len=*), intent(in), optional :: method !! override interpolation method
        integer, intent(in), optional :: nu(:) !! partial derivative order along each grid axis
        real(dp), allocatable :: y(:,:)
        real(dp), allocatable :: yi(:)
        integer :: i, code, first_error
        type(regular_grid_interpolator_type) :: override_grid
        allocate(y(size(x,1),self%channels()))
        y = ieee_value(0.0_dp,ieee_quiet_nan)
        first_error = grid_status_success
        if (present(method) .and. self%valid()) then
            if (is_spline_method(method) .and. method /= self%method) then
                call initialize_grid(override_grid, self%axes, self%values, method, self%bounds_error, &
                    self%fill_value, .not. self%use_fill, code)
                if (code == grid_status_success) then
                    y = override_grid%evaluate(x, status=code, nu=nu)
                end if
                if (present(status)) status = code
                return
            end if
        end if
        do i = 1, size(x,1)
            yi = self%evaluate_channels(x(i,:),status=code,method=method,nu=nu)
            if (code == grid_status_success) y(i,:) = yi
            if (code /= grid_status_success .and. first_error == grid_status_success) first_error = code
        end do
        if (.not. self%valid() .or. size(x,2) /= self%dimensions()) first_error = grid_status_invalid_input
        if (present(status)) status = first_error
    end function evaluate_many

    function evaluate_scalar_batch(self, x, status, method, nu) result(y)
        class(regular_grid_interpolator_type), intent(in) :: self !! configured single-channel grid
        real(dp), intent(in) :: x(:,:) !! row-oriented queries, shape (nquery,ndim)
        integer, intent(out), optional :: status !! zero if all requests succeed
        character(len=*), intent(in), optional :: method !! nearest, linear, pchip, slinear, cubic or quintic
        integer, intent(in), optional :: nu(:) !! partial derivative order along each grid axis
        real(dp), allocatable :: y(:)
        integer :: code
        real(dp), allocatable :: matrix(:,:)
        matrix = self%evaluate(x,status=code,method=method,nu=nu)
        allocate(y(size(x,1)))
        y = ieee_value(0.0_dp,ieee_quiet_nan)
        if (self%channels() == 1) then
            y = matrix(:,1)
        else
            code = grid_status_invalid_input
        end if
        if (present(status)) status = code
    end function evaluate_scalar_batch

    recursive function evaluate_all_channels(self, x, status, method, nu) result(y)
        class(regular_grid_interpolator_type), intent(in) :: self !! configured grid
        real(dp), intent(in) :: x(:) !! coordinates of one interpolation point
        integer, intent(out), optional :: status !! zero on success
        character(len=*), intent(in), optional :: method !! interpolation method override
        integer, intent(in), optional :: nu(:) !! partial derivative order along each grid axis
        real(dp), allocatable :: y(:)
        integer :: idx(max_grid_dim), k, mask, offset, nch, nd, code, ncorners, loc
        type(regular_grid_interpolator_type) :: override_grid
        real(dp) :: frac(max_grid_dim), weight, tk, q, lo, hi
        logical :: outside, nearest
        character(len=7) :: how

        code = grid_status_invalid_input
        nch = self%channels()
        allocate(y(nch))
        y = ieee_value(0.0_dp,ieee_quiet_nan)
        if (.not. self%valid()) goto 900
        nd = self%dimensions()
        if (size(x) /= nd) goto 900
        how = self%method
        if (present(method)) then
            if (.not. supported_grid_method(method)) goto 900
            how = method
        end if
        ! Method changes from a previously constructed interpolation plan may
        ! require a new spline coefficient fit. It is transient rather than
        ! modifying an intent(in) object, matching SciPy's per-call override.
        if (present(nu)) then
            if (.not. is_spline_method(how)) goto 900
            if (size(nu) /= nd) goto 900
            if (any(nu < 0)) goto 900
        end if
        if (is_spline_method(how) .and. how /= self%method) then
            call initialize_grid(override_grid, self%axes, self%values, how, self%bounds_error, &
                self%fill_value, .not. self%use_fill, code)
            if (code == grid_status_success) y = override_grid%evaluate_channels(x, status=code, nu=nu)
            goto 900
        end if
        code = grid_status_success
        if (how == 'pchip' .and. self%method /= 'pchip') then
            if (any(self%lengths < 4)) then
                code = grid_status_invalid_input
                goto 900
            end if
            if (.not. all(ieee_is_finite(self%values))) then
                code = grid_status_invalid_input
                goto 900
            end if
        end if
        if (any(ieee_is_nan(x))) goto 900
        outside = .false.
        do k = 1, nd
            loc = axis_lower_cell(self%axes(k)%points, x(k))
            idx(k) = loc
            if (self%lengths(k) == 1) then
                frac(k) = 0.0_dp
            else
                lo = self%axes(k)%points(loc)
                hi = self%axes(k)%points(loc+1)
                frac(k) = (x(k)-lo)/(hi-lo)
            end if
            outside = outside .or. x(k) < min(self%axes(k)%points(1), &
                self%axes(k)%points(self%lengths(k))) .or. &
                x(k) > max(self%axes(k)%points(1), &
                self%axes(k)%points(self%lengths(k)))
        end do
        if (outside) then
            if (self%bounds_error) then
                code = grid_status_out_of_bounds
                goto 900
            end if
            if (self%use_fill) then
                y = self%fill_value
                goto 900
            end if
        end if
        if (is_spline_method(how)) then
            call evaluate_tensor_spline(self, x, y, nu)
            goto 900
        end if
        if (how == 'pchip') then
            call evaluate_tensor_pchip(self, x, y)
            goto 900
        end if
        nearest = how == 'nearest'
        y = 0.0_dp
        if (nearest) then
            offset = 1
            do k = 1, nd
                loc = idx(k)
                if (self%lengths(k) > 1) then
                    q = x(k)
                    lo = self%axes(k)%points(loc)
                    hi = self%axes(k)%points(loc+1)
                    if (lo < hi) then
                        if ((q-lo)/(hi-lo) > 0.5_dp) loc = loc+1
                    else
                        ! SciPy reverses descending axes before its nearest
                        ! tie-break: ties select the smaller physical knot.
                        if ((q-hi)/(lo-hi) <= 0.5_dp) loc = loc+1
                    end if
                end if
                offset = offset + (loc-1)*self%strides(k)
            end do
            y = self%values(offset,:)
        else
            ncorners = 2**nd
            do mask = 0, ncorners-1
                offset = 1
                weight = 1.0_dp
                do k = 1, nd
                    loc = idx(k)
                    tk = frac(k)
                    if (btest(mask,k-1) .and. self%lengths(k) > 1) then
                        loc = loc+1
                        weight = weight*tk
                    else if (self%lengths(k) > 1) then
                        weight = weight*(1.0_dp-tk)
                    else if (btest(mask,k-1)) then
                        weight = 0.0_dp
                    end if
                    offset = offset+(loc-1)*self%strides(k)
                end do
                y = y+weight*self%values(offset,:)
            end do
        end if
900     if (present(status)) status = code
    end function evaluate_all_channels

    integer function axis_lower_cell(axis, x) result(idx)
        real(dp), intent(in) :: axis(:) !! increasing or decreasing coordinate vector
        real(dp), intent(in) :: x !! query coordinate
        integer :: n, lo, hi, mid
        logical :: ascending
        n = size(axis)
        idx = 1
        if (n <= 1) return
        ascending = axis(n) > axis(1)
        if (ascending) then
            if (x <= axis(1)) return
            if (x >= axis(n)) then
                idx = n-1
                return
            end if
        else
            if (x >= axis(1)) return
            if (x <= axis(n)) then
                idx = n-1
                return
            end if
        end if
        lo = 1
        hi = n
        do while (hi-lo > 1)
            mid = (lo+hi)/2
            if (ascending .eqv. (x >= axis(mid))) then
                lo = mid
            else
                hi = mid
            end if
        end do
        idx = lo
    end function axis_lower_cell

    logical function valid_complex_grid(self)
        class(complex_regular_grid_interpolator_type), intent(in) :: self !! grid to validate
        valid_complex_grid = self%real_part%valid() .and. self%imag_part%valid()
    end function valid_complex_grid

    integer function complex_grid_dimensions(self)
        class(complex_regular_grid_interpolator_type), intent(in) :: self !! grid object
        complex_grid_dimensions = self%real_part%dimensions()
    end function complex_grid_dimensions

    integer function complex_grid_channels(self)
        class(complex_regular_grid_interpolator_type), intent(in) :: self !! grid object
        complex_grid_channels = self%real_part%channels()
    end function complex_grid_channels

    function evaluate_complex_channels(self, x, status, method, nu) result(y)
        class(complex_regular_grid_interpolator_type), intent(in) :: self !! complex-valued interpolator
        real(dp), intent(in) :: x(:) !! query coordinates
        integer, intent(out), optional :: status !! zero on success
        character(len=*), intent(in), optional :: method !! nearest, linear, pchip, slinear, cubic or quintic
        integer, intent(in), optional :: nu(:) !! partial derivative order along each grid axis
        complex(dp), allocatable :: y(:)
        real(dp), allocatable :: a(:), b(:)
        integer :: code, code2
        if (present(method)) then
            if (method == 'pchip') then
                allocate(y(self%channels()))
                y = cmplx(ieee_value(0.0_dp,ieee_quiet_nan),0.0_dp,kind=dp)
                if (present(status)) status = grid_status_invalid_input
                return
            end if
        end if
        a = self%real_part%evaluate_channels(x,status=code,method=method,nu=nu)
        b = self%imag_part%evaluate_channels(x,status=code2,method=method,nu=nu)
        y = cmplx(a,b,kind=dp)
        if (any(ieee_is_nan(x))) then
            y = cmplx(ieee_value(0.0_dp,ieee_quiet_nan),0.0_dp,kind=dp)
        end if
        if (code == grid_status_success) code = code2
        if (present(status)) status = code
    end function evaluate_complex_channels

    function evaluate_complex_one(self, x, status, method, nu) result(y)
        class(complex_regular_grid_interpolator_type), intent(in) :: self !! complex-valued interpolator
        real(dp), intent(in) :: x(:) !! query coordinates
        integer, intent(out), optional :: status !! zero on success
        character(len=*), intent(in), optional :: method !! nearest, linear, pchip, slinear, cubic or quintic
        integer, intent(in), optional :: nu(:) !! partial derivative order along each grid axis
        complex(dp) :: y
        complex(dp), allocatable :: z(:)
        integer :: code
        y = cmplx(ieee_value(0.0_dp,ieee_quiet_nan), &
            ieee_value(0.0_dp,ieee_quiet_nan),kind=dp)
        code = grid_status_invalid_input
        if (self%channels() /= 1) goto 900
        z = self%evaluate_channels(x,status=code,method=method,nu=nu)
        if (code == grid_status_success) y = z(1)
900     if (present(status)) status = code
    end function evaluate_complex_one

    function evaluate_complex_many(self, x, status, method, nu) result(y)
        class(complex_regular_grid_interpolator_type), intent(in) :: self !! complex grid
        real(dp), intent(in) :: x(:,:) !! queries shaped (nqueries,ndim)
        integer, intent(out), optional :: status !! first nonzero query status
        character(len=*), intent(in), optional :: method !! nearest, linear, pchip, slinear, cubic or quintic
        integer, intent(in), optional :: nu(:) !! partial derivative order along each grid axis
        complex(dp), allocatable :: y(:,:)
        real(dp), allocatable :: a(:,:), b(:,:)
        integer :: code, code2, i
        if (present(method)) then
            if (method == 'pchip') then
                allocate(y(size(x,1),self%channels()))
                y = cmplx(ieee_value(0.0_dp,ieee_quiet_nan),0.0_dp,kind=dp)
                if (present(status)) status = grid_status_invalid_input
                return
            end if
        end if
        a = self%real_part%evaluate(x,status=code,method=method,nu=nu)
        b = self%imag_part%evaluate(x,status=code2,method=method,nu=nu)
        y = cmplx(a,b,kind=dp)
        do i = 1, size(x,1)
            if (any(ieee_is_nan(x(i,:)))) &
                y(i,:) = cmplx(ieee_value(0.0_dp,ieee_quiet_nan),0.0_dp,kind=dp)
        end do
        if (code == grid_status_success) code = code2
        if (present(status)) status = code
    end function evaluate_complex_many

    function evaluate_complex_values(self, x, status, method, nu) result(y)
        class(complex_regular_grid_interpolator_type), intent(in) :: self !! complex grid
        real(dp), intent(in) :: x(:,:) !! queries shaped (nqueries,ndim)
        integer, intent(out), optional :: status !! zero on success
        character(len=*), intent(in), optional :: method !! nearest, linear, pchip, slinear, cubic or quintic
        integer, intent(in), optional :: nu(:) !! partial derivative order along each grid axis
        complex(dp), allocatable :: y(:)
        complex(dp), allocatable :: m(:,:)
        integer :: code
        m = self%evaluate(x,status=code,method=method,nu=nu)
        allocate(y(size(x,1)))
        y = cmplx(ieee_value(0.0_dp,ieee_quiet_nan), &
            ieee_value(0.0_dp,ieee_quiet_nan),kind=dp)
        if (self%channels() == 1) then
            y = m(:,1)
        else
            code = grid_status_invalid_input
        end if
        if (present(status)) status = code
    end function evaluate_complex_values

    function interpn_real(axes, values, x, method, bounds_error, fill_value, extrapolate, status) result(y)
        type(grid_axis), intent(in) :: axes(:) !! monotone coordinates for each dimension
        real(dp), intent(in) :: values(..) !! N-dimensional data or rank N+1 with channels
        real(dp), intent(in) :: x(:,:) !! interpolation points, shape (nquery,ndim)
        character(len=*), intent(in), optional :: method !! nearest, linear, pchip, slinear, cubic or quintic
        logical, intent(in), optional :: bounds_error !! flag out-of-domain coordinates
        real(dp), intent(in), optional :: fill_value !! out-of-domain value
        logical, intent(in), optional :: extrapolate !! permit extrapolation instead of fill
        integer, intent(out), optional :: status !! zero if all operations succeed
        real(dp), allocatable :: y(:,:)
        type(regular_grid_interpolator_type) :: grid
        integer :: code
        grid = regular_grid_interpolator(axes,values,method,bounds_error,fill_value,extrapolate,code)
        y = grid%evaluate(x,status=code)
        if (present(status)) status = code
    end function interpn_real

    function interpn_complex(axes, values, x, method, bounds_error, fill_value, extrapolate, status) result(y)
        type(grid_axis), intent(in) :: axes(:) !! monotone coordinates for each dimension
        complex(dp), intent(in) :: values(..) !! N-dimensional complex samples, optional channels
        real(dp), intent(in) :: x(:,:) !! interpolation points, shape (nquery,ndim)
        character(len=*), intent(in), optional :: method !! nearest, linear, pchip, slinear, cubic or quintic
        logical, intent(in), optional :: bounds_error !! flag out-of-domain queries
        complex(dp), intent(in), optional :: fill_value !! out-of-domain value
        logical, intent(in), optional :: extrapolate !! extrapolate from edge cells if true
        integer, intent(out), optional :: status !! zero if all operations succeed
        complex(dp), allocatable :: y(:,:)
        type(complex_regular_grid_interpolator_type) :: grid
        integer :: code
        grid = complex_regular_grid_interpolator(axes,values,method,bounds_error,fill_value,extrapolate,code)
        y = grid%evaluate(x,status=code)
        if (present(status)) status = code
    end function interpn_complex

end module scifort_interpolate_regular_grid
