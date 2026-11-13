! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! One-dimensional cubic interpolation, including global cubic splines.
! Independent implementation from the formulas of Fritsch & Butland (1984),
! Akima (1970), and the modified Akima slope definition (Moler & Ionita).
! SciPy is used as a numerical test reference, not as implementation source.
module scifort_interpolate
    use scifort_interpolate_sphere, only : lsq_sphere_bivariate_spline, &
        sphere_bivariate_spline_type, sphere_spline_success, sphere_spline_invalid_input, &
        sphere_spline_rank_deficient, sphere_spline_numerical_failure, sphere_spline_size_limit
    use scifort_interpolate_aaa, only : aaa,aaa_type,complex_aaa_type,aaa_success, &
        aaa_invalid_input,aaa_max_terms,aaa_numerical_failure,aaa_spectral_size_limit, &
        aaa_setup_failure,aaa_undefined
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_interpolate_regular_grid, only : grid_axis, regular_grid_interpolator_type, &
        complex_regular_grid_interpolator_type, regular_grid_interpolator, complex_regular_grid_interpolator, &
        interpn, complex_interpn, grid_status_success, grid_status_invalid_input, &
        grid_status_out_of_bounds, grid_status_singular
    use scifort_interpolate_smoothing, only : make_smoothing_spline, smoothing_status_success, &
        smoothing_status_invalid_input, smoothing_status_singular
    use scifort_interpolate_fitpack_bivariate, only : fitpack_bivariate_tck_type, &
        bisplev, bisplev_points, bispline_integral, rect_bivariate_spline, &
        bispline_partial_derivative, bispline_insert_knot, &
        bispline_success, bispline_invalid_input, bispline_singular
    use scifort_interpolate_fitpack_compat, only : fitpack_tck_type, fitpack_parametric_tck_type, splrep, splprep, splev, splint, splder, &
        splantider, spalde, fitpack_tck_success, fitpack_tck_invalid_input, &
        fitpack_tck_unsupported_smoothing, fitpack_tck_singular, fitpack_tck_out_of_bounds
    use scifort_interpolate_bspline, only : bspline_interpolator_type, bspline_array_interpolator_type, &
        bspline, make_interp_spline, make_lsq_spline, bspline_design_matrix, lsq_univariate_spline, &
        interpolated_univariate_spline, bspline_status_success, &
        bspline_status_invalid_input, bspline_status_singular
    use scifort_interpolate_legacy, only : interp1d_type, complex_interp1d_type, &
        interp1d, interp1d_success, interp1d_invalid_input, interp1d_out_of_bounds
    use scifort_interpolate_rbf, only : rbf_interpolator_type, complex_rbf_interpolator_type, &
        rbf_interpolator, rbf_status_success, rbf_status_invalid_input, rbf_status_singular
    use scifort_interpolate_global, only : krogh_interpolator_type, complex_krogh_interpolator_type, &
        floater_hormann_interpolator_type, complex_floater_hormann_interpolator_type, &
        krogh_interpolator, krogh_interpolate, floater_hormann_interpolator, &
        global_interp_success, global_interp_invalid_input
    use scifort_interpolate_barycentric, only : barycentric_interpolator_type, &
        complex_barycentric_interpolator_type, barycentric_interpolator, &
        barycentric_interpolate, barycentric_success, barycentric_invalid_input
    use scifort_interpolate_piecewise, only : ppoly_type, bpoly_type, &
        complex_ppoly_type, complex_bpoly_type, ppoly, bpoly, &
        ppoly_from_bernstein, bpoly_from_power, bpoly_from_derivatives, &
        complex_ppoly_from_bernstein, complex_bpoly_from_power, &
        piecewise_native_real_roots, piecewise_status_success, piecewise_status_invalid_input
    use scifort_interpolate_kdtree, only : kd_tree_type, kd_tree, kdtree_success, kdtree_invalid_input
    use scifort_interpolate_clough_tocher, only : clough_tocher_2d_interpolator_type, &
        complex_clough_tocher_2d_interpolator_type, clough_tocher_2d_interpolator, &
        ct_status_success, ct_status_invalid_input, ct_status_degenerate, ct_status_nonconverged
    use scifort_interpolate_scattered, only : nearest_nd_interpolator_type, &
        complex_nearest_nd_interpolator_type, linear_nd_interpolator_type, &
        complex_linear_nd_interpolator_type, nearest_nd_interpolator, linear_nd_interpolator, scattered_griddata => griddata, &
        scattered_status_success, scattered_status_invalid_input, scattered_status_degenerate, delaunay_nd
    implicit none
    private

    public :: lsq_sphere_bivariate_spline, sphere_bivariate_spline_type
    public :: sphere_spline_success, sphere_spline_invalid_input, sphere_spline_rank_deficient
    public :: sphere_spline_numerical_failure, sphere_spline_size_limit

    public :: aaa,aaa_type,complex_aaa_type,aaa_success,aaa_invalid_input,aaa_max_terms
    public :: aaa_numerical_failure,aaa_spectral_size_limit,aaa_setup_failure,aaa_undefined

    integer, parameter, public :: interpolate_status_success = 0
    integer, parameter, public :: interpolate_status_invalid_input = 1

    type, public :: piecewise_cubic_interpolator
        private
        real(dp), allocatable :: x(:)
        real(dp), allocatable :: c(:,:)  ! c(1:4,j): constant through cubic coefficient
        real(dp), allocatable :: area(:) ! cumulative exact integral at each knot
        logical :: allow_extrapolation = .false.
        logical :: periodic_extension = .false.
    contains
        procedure, private :: evaluate_scalar
        procedure, private :: evaluate_array
        generic :: evaluate => evaluate_scalar, evaluate_array
        procedure, private :: derivative_scalar
        procedure, private :: derivative_array
        generic :: derivative => derivative_scalar, derivative_array
        procedure :: integrate => integrate_piecewise
        procedure :: valid => valid_piecewise
        procedure :: knots => get_piecewise_knots
        procedure :: coefficients => get_piecewise_coefficients
    end type piecewise_cubic_interpolator

    type, public :: piecewise_cubic_array_interpolator
        private
        type(piecewise_cubic_interpolator), allocatable :: columns(:)
    contains
        procedure, private :: evaluate_matrix_scalar
        procedure, private :: evaluate_matrix_array
        generic :: evaluate => evaluate_matrix_scalar, evaluate_matrix_array
        procedure, private :: derivative_matrix_scalar
        procedure, private :: derivative_matrix_array
        generic :: derivative => derivative_matrix_scalar, derivative_matrix_array
        procedure :: integrate => integrate_matrix
        procedure :: valid => valid_matrix
        procedure :: coefficients => get_matrix_coefficients
    end type piecewise_cubic_array_interpolator

    public :: kd_tree_type, kd_tree, kdtree_success, kdtree_invalid_input
    public :: nearest_nd_interpolator_type, complex_nearest_nd_interpolator_type
    public :: linear_nd_interpolator_type, complex_linear_nd_interpolator_type
    public :: nearest_nd_interpolator, linear_nd_interpolator, griddata, delaunay_nd
    public :: clough_tocher_2d_interpolator_type, complex_clough_tocher_2d_interpolator_type
    public :: clough_tocher_2d_interpolator
    public :: ct_status_success, ct_status_invalid_input, ct_status_degenerate, ct_status_nonconverged
    public :: scattered_status_success, scattered_status_invalid_input, scattered_status_degenerate
    public :: pchip_interpolator, akima1d_interpolator, cubic_hermite_spline, cubic_spline, pchip_interpolate
    public :: make_smoothing_spline, smoothing_status_success, smoothing_status_invalid_input, smoothing_status_singular
    public :: interp1d_type, complex_interp1d_type, interp1d
    public :: interp1d_success, interp1d_invalid_input, interp1d_out_of_bounds
    public :: rbf_interpolator_type, complex_rbf_interpolator_type, rbf_interpolator
    public :: rbf_status_success, rbf_status_invalid_input, rbf_status_singular
    public :: fitpack_bivariate_tck_type, bisplev, bisplev_points, bispline_integral, rect_bivariate_spline, &
        bispline_partial_derivative, bispline_insert_knot, &
        bispline_success, bispline_invalid_input, bispline_singular
    public :: fitpack_tck_type, fitpack_parametric_tck_type, splrep, splprep, &
        splev, splint, splder, splantider, spalde
    public :: fitpack_tck_success, fitpack_tck_invalid_input, fitpack_tck_unsupported_smoothing, &
        fitpack_tck_singular, fitpack_tck_out_of_bounds
    public :: bspline_interpolator_type, bspline_array_interpolator_type, bspline, make_interp_spline, make_lsq_spline, &
        lsq_univariate_spline, interpolated_univariate_spline, &
        bspline_design_matrix, bspline_status_success, bspline_status_invalid_input, bspline_status_singular
    public :: barycentric_interpolator_type, complex_barycentric_interpolator_type
    public :: barycentric_interpolator, barycentric_interpolate, barycentric_success, barycentric_invalid_input
    public :: krogh_interpolator_type, complex_krogh_interpolator_type, krogh_interpolator, krogh_interpolate
    public :: floater_hormann_interpolator_type, complex_floater_hormann_interpolator_type, &
        floater_hormann_interpolator, global_interp_success, global_interp_invalid_input
    public :: ppoly_type, bpoly_type, complex_ppoly_type, complex_bpoly_type
    public :: ppoly, bpoly, ppoly_from_bernstein, bpoly_from_power, bpoly_from_derivatives
    public :: complex_ppoly_from_bernstein, complex_bpoly_from_power
    public :: piecewise_status_success, piecewise_status_invalid_input, piecewise_native_real_roots
    public :: grid_axis, regular_grid_interpolator_type, complex_regular_grid_interpolator_type
    public :: regular_grid_interpolator, complex_regular_grid_interpolator, interpn, complex_interpn
    public :: grid_status_success, grid_status_invalid_input, grid_status_out_of_bounds, grid_status_singular

    public :: ppoly_from_cubic, ppoly_from_spline
    interface griddata
        module procedure griddata_real_vector_ct, griddata_real_matrix_ct
        module procedure griddata_complex_vector_ct, griddata_complex_matrix_ct
    end interface

    interface ppoly_from_spline
        module procedure ppoly_from_spline_scalar, ppoly_from_spline_matrix
    end interface

    interface ppoly_from_cubic
        module procedure ppoly_from_cubic_scalar, ppoly_from_cubic_matrix
    end interface

    interface cubic_spline
        module procedure cubic_spline_1d
        module procedure cubic_spline_matrix
    end interface

    interface pchip_interpolate
        module procedure pchip_interpolate_scalar
        module procedure pchip_interpolate_array
    end interface

contains

    ! Extend the existing nearest/linear griddata wrappers with SciPy's
    ! 2-D cubic (Clough-Tocher) method without introducing a circular module
    ! dependency between the triangulation and patch implementations.
    function griddata_real_vector_ct(x,y,xi,method,fill_value,rescale,status) result(v)
        real(dp), intent(in) :: x(:,:) !! observation sites, (sites,dimensions)
        real(dp), intent(in) :: y(:) !! scalar ordinate at each site
        real(dp), intent(in) :: xi(:,:) !! query points, (queries,dimensions)
        character(len=*), intent(in), optional :: method !! cubic, linear (default), or nearest
        real(dp), intent(in), optional :: fill_value !! out-of-hull value for cubic/linear
        logical, intent(in), optional :: rescale !! rescale coordinate ranges
        integer, intent(out), optional :: status !! status result
        real(dp), allocatable :: v(:),tmp(:,:)
        tmp=griddata_real_matrix_ct(x,reshape(y,[size(y),1]),xi,method,fill_value,rescale,status)
        if (size(tmp,2)==0) then
            allocate(v(0))
        else
            v=tmp(:,1)
        end if
    end function griddata_real_vector_ct

    function griddata_real_matrix_ct(x,y,xi,method,fill_value,rescale,status) result(v)
        real(dp), intent(in) :: x(:,:) !! observation sites, (sites,dimensions)
        real(dp), intent(in) :: y(:,:) !! multichannel ordinates, (sites,channels)
        real(dp), intent(in) :: xi(:,:) !! query points, (queries,dimensions)
        character(len=*), intent(in), optional :: method !! cubic, linear (default), nearest
        real(dp), intent(in), optional :: fill_value !! convex-hull fill value
        logical, intent(in), optional :: rescale !! rescale coordinate ranges
        integer, intent(out), optional :: status !! success, invalid, or degenerate
        real(dp), allocatable :: v(:,:)
        type(clough_tocher_2d_interpolator_type) :: cubic
        integer :: code,initial_code
        character(len=:), allocatable :: method_name
        method_name='linear'
        if (present(method)) method_name=trim(method)
        if (method_name=='cubic') then
            cubic=clough_tocher_2d_interpolator(x,y,fill_value=fill_value,rescale=rescale,status=initial_code)
            allocate(v(0,0))
            code=initial_code
            if (initial_code==ct_status_success .or. initial_code==ct_status_nonconverged) then
                v=cubic%evaluate(xi,status=code)
                if (code==ct_status_success .and. initial_code==ct_status_nonconverged) code=initial_code
            end if
        else
            v=scattered_griddata(x,y,xi,method,fill_value,rescale,status=code)
        end if
        if (present(status)) status=code
    end function griddata_real_matrix_ct

    function griddata_complex_vector_ct(x,y,xi,method,fill_value,rescale,status) result(v)
        real(dp), intent(in) :: x(:,:) !! observation sites, (sites,dimensions)
        complex(dp), intent(in) :: y(:) !! complex ordinate at each site
        real(dp), intent(in) :: xi(:,:) !! query points, (queries,dimensions)
        character(len=*), intent(in), optional :: method !! cubic, linear, nearest
        complex(dp), intent(in), optional :: fill_value !! out-of-hull fill
        logical, intent(in), optional :: rescale !! rescale coordinate ranges
        integer, intent(out), optional :: status !! result status
        complex(dp), allocatable :: v(:),tmp(:,:)
        tmp=griddata_complex_matrix_ct(x,reshape(y,[size(y),1]),xi,method,fill_value,rescale,status)
        if (size(tmp,2)==0) then
            allocate(v(0))
        else
            v=tmp(:,1)
        end if
    end function griddata_complex_vector_ct

    function griddata_complex_matrix_ct(x,y,xi,method,fill_value,rescale,status) result(v)
        real(dp), intent(in) :: x(:,:) !! scattered data positions
        complex(dp), intent(in) :: y(:,:) !! multichannel complex ordinates
        real(dp), intent(in) :: xi(:,:) !! query points
        character(len=*), intent(in), optional :: method !! cubic, linear, or nearest
        complex(dp), intent(in), optional :: fill_value !! out-of-hull fill
        logical, intent(in), optional :: rescale !! scale spatial axes
        integer, intent(out), optional :: status !! status code
        complex(dp), allocatable :: v(:,:)
        type(complex_clough_tocher_2d_interpolator_type) :: cubic
        integer :: code,initial_code
        character(len=:), allocatable :: method_name
        method_name='linear'
        if (present(method)) method_name=trim(method)
        if (method_name=='cubic') then
            cubic=clough_tocher_2d_interpolator(x,y,fill_value=fill_value,rescale=rescale,status=initial_code)
            allocate(v(0,0))
            code=initial_code
            if (initial_code==ct_status_success .or. initial_code==ct_status_nonconverged) then
                v=cubic%evaluate(xi,status=code)
                if (code==ct_status_success .and. initial_code==ct_status_nonconverged) code=initial_code
            end if
        else
            v=scattered_griddata(x,y,xi,method,fill_value,rescale,status=code)
        end if
        if (present(status)) status=code
    end function griddata_complex_matrix_ct


    ! SciPy's PPoly.from_spline: construct Taylor/power coefficients at
    ! every B-spline knot, retaining repeated knots and the full knot array.
    function ppoly_from_spline_scalar(spline,extrapolate,periodic,status) result(poly)
        type(bspline_interpolator_type), intent(in) :: spline !! source 1D B-spline
        logical, intent(in), optional :: extrapolate !! allow endpoint polynomial extension (default true)
        logical, intent(in), optional :: periodic !! periodic PPoly evaluation (default false)
        integer, intent(out), optional :: status !! zero success, 1 invalid input
        type(ppoly_type) :: poly
        real(dp), allocatable :: c(:,:), knots(:)
        real(dp) :: factorial
        integer :: k, j, d, code
        logical :: ext, per
        code = piecewise_status_invalid_input
        if (.not. spline%valid()) goto 900
        k = spline%degree()
        knots = spline%knots()
        allocate(c(k+1,size(knots)-1))
        do j = 1,size(knots)-1
            factorial = 1.0_dp
            do d = 0,k
                if (d>0) factorial = factorial*real(d,dp)
                c(k+1-d,j) = spline%derivative(knots(j),nu=d)/factorial
            end do
        end do
        ext = .true.
        per = .false.
        if (present(extrapolate)) ext = extrapolate
        if (present(periodic)) per = periodic
        poly = ppoly(c,knots,extrapolate=ext,periodic=per,status=code)
900     if (present(status)) status = code
    end function ppoly_from_spline_scalar

    function ppoly_from_spline_matrix(spline,extrapolate,periodic,status) result(poly)
        type(bspline_array_interpolator_type), intent(in) :: spline !! multichannel B-spline
        logical, intent(in), optional :: extrapolate !! allow endpoint extrapolation (default true)
        logical, intent(in), optional :: periodic !! use periodic PPoly extrapolation (default false)
        integer, intent(out), optional :: status !! zero success, 1 invalid input
        type(ppoly_type) :: poly
        real(dp), allocatable :: c(:,:,:), knots(:), bs_coeff(:,:)
        real(dp) :: factorial
        integer :: k, j, d, channels, code
        logical :: ext, per
        code = piecewise_status_invalid_input
        if (.not. spline%valid()) goto 900
        knots = spline%knots()
        channels = spline%channels()
        ! The multichannel B-spline type uses the same knots and degree
        ! in every channel. Find the common degree via coefficient counts.
        bs_coeff = spline%coefficients()
        k = size(knots)-size(bs_coeff,1)-1
        allocate(c(k+1,size(knots)-1,channels))
        do j = 1,size(knots)-1
            factorial = 1.0_dp
            do d = 0,k
                if (d>0) factorial = factorial*real(d,dp)
                c(k+1-d,j,:) = spline%derivative(knots(j),nu=d)/factorial
            end do
        end do
        ext = .true.
        per = .false.
        if (present(extrapolate)) ext = extrapolate
        if (present(periodic)) per = periodic
        poly = ppoly(c,knots,extrapolate=ext,periodic=per,status=code)
900     if (present(status)) status = code
    end function ppoly_from_spline_matrix

    ! Bridge the existing SciFort cubic interpolators into SciPy-style PPoly.
    function ppoly_from_cubic_scalar(spline) result(poly)
        type(piecewise_cubic_interpolator), intent(in) :: spline !! original cubic interpolator
        type(ppoly_type) :: poly
        real(dp), allocatable :: coef(:,:), knots(:)
        if (.not. spline%valid()) return
        coef=spline%coefficients()
        knots=spline%knots()
        poly=ppoly(coef(4:1:-1,:),knots,extrapolate=spline%allow_extrapolation,&
            periodic=spline%periodic_extension)
    end function ppoly_from_cubic_scalar

    function ppoly_from_cubic_matrix(spline) result(poly)
        type(piecewise_cubic_array_interpolator), intent(in) :: spline !! multiple cubic data channels
        type(ppoly_type) :: poly
        real(dp), allocatable :: coef(:,:,:), knots(:)
        if (.not. spline%valid()) return
        coef=spline%coefficients()
        knots=spline%columns(1)%knots()
        poly=ppoly(coef(4:1:-1,:,:),knots,&
            extrapolate=spline%columns(1)%allow_extrapolation,&
            periodic=spline%columns(1)%periodic_extension)
    end function ppoly_from_cubic_matrix

    ! CubicSpline (scipy.interpolate.CubicSpline) for real, one-dimensional data.
    ! Endpoint values in natural/clamped/first/second conditions refer to
    ! derivatives of the appropriate order. Unlike PCHIP, this is globally C2.
    function cubic_spline_1d(x, y, bc_type, status, extrapolate, bc_left, bc_right, &
        left_value, right_value) result(poly)
        real(dp), intent(in) :: x(:) !! strictly increasing, finite knot positions, at least two
        real(dp), intent(in) :: y(:) !! finite ordinate values, same size as x
        character(len=*), intent(in), optional :: bc_type !! both endpoints: not-a-knot (default), natural, clamped, periodic
        integer, intent(out), optional :: status !! success or invalid-input status
        logical, intent(in), optional :: extrapolate !! override extrapolation policy; periodic is default for periodic splines
        character(len=*), intent(in), optional :: bc_left !! override left boundary: not-a-knot, natural, clamped, first, second
        character(len=*), intent(in), optional :: bc_right !! override right boundary condition
        real(dp), intent(in), optional :: left_value !! left derivative value for first, second, or clamped (default 0)
        real(dp), intent(in), optional :: right_value !! right derivative value for first, second, or clamped (default 0)
        type(piecewise_cubic_interpolator) :: poly
        real(dp), allocatable :: h(:), delta(:), slopes(:), lower(:), diag(:), upper(:), rhs(:)
        real(dp) :: left_d, right_d, h1, h2, total, curvature
        character(len=:), allocatable :: left_bc, right_bc
        integer :: n, i, code
        logical :: periodic

        code = interpolate_status_invalid_input
        left_bc = 'not-a-knot'
        right_bc = 'not-a-knot'
        if (present(bc_type)) then
            left_bc = trim(bc_type)
            right_bc = left_bc
        end if
        if (present(bc_left)) left_bc = trim(bc_left)
        if (present(bc_right)) right_bc = trim(bc_right)
        periodic = (left_bc == 'periodic' .and. right_bc == 'periodic')
        if ((left_bc == 'periodic' .or. right_bc == 'periodic') .and. .not. periodic) goto 900
        if (.not. valid_cubic_boundary(left_bc) .or. .not. valid_cubic_boundary(right_bc)) goto 900
        left_d = 0.0_dp
        right_d = 0.0_dp
        if (present(left_value)) left_d = left_value
        if (present(right_value)) right_d = right_value
        if (.not. ieee_is_finite(left_d) .or. .not. ieee_is_finite(right_d)) goto 900
        if (periodic .and. (present(left_value) .or. present(right_value))) goto 900
        call validate_samples(x, y, h, delta, code)
        if (code /= interpolate_status_success) goto 900
        code = interpolate_status_invalid_input
        n = size(x)
        allocate(slopes(n))
        if (periodic) then
            if (abs(y(1)-y(n)) > 1.0e-15_dp * max(1.0_dp,abs(y(1)),abs(y(n)))) goto 900
            call periodic_cubic_slopes(h, delta, slopes, code)
            if (code /= interpolate_status_success) goto 900
        else
            ! For two knots, SciPy replaces not-a-knot by a first derivative
            ! equal to the secant slope; for three double not-a-knot, use the
            ! unique quadratic passing through the three knots.
            if (n == 2) then
                if (left_bc == 'not-a-knot') then
                    left_bc = 'first'
                    left_d = delta(1)
                end if
                if (right_bc == 'not-a-knot') then
                    right_bc = 'first'
                    right_d = delta(1)
                end if
            end if
            if (n == 3 .and. left_bc == 'not-a-knot' .and. right_bc == 'not-a-knot') then
                curvature = (delta(2) - delta(1))/(h(1)+h(2))
                slopes(1) = delta(1)-h(1)*curvature
                slopes(2) = delta(1)+h(1)*curvature
                slopes(3) = delta(2)+h(2)*curvature
                code = interpolate_status_success
            else
                allocate(lower(n), diag(n), upper(n), rhs(n))
                lower = 0.0_dp
                diag = 0.0_dp
                upper = 0.0_dp
                rhs = 0.0_dp
                do i = 2, n-1
                    lower(i) = h(i)
                    diag(i) = 2.0_dp*(h(i-1)+h(i))
                    upper(i) = h(i-1)
                    rhs(i) = 3.0_dp*(h(i)*delta(i-1)+h(i-1)*delta(i))
                end do
                select case (left_bc)
                case ('not-a-knot')
                    h1 = h(1)
                    h2 = h(2)
                    total = h1+h2
                    diag(1) = h2
                    upper(1) = total
                    rhs(1) = ((h1+2.0_dp*total)*h2*delta(1) + h1*h1*delta(2))/total
                case ('natural', 'second')
                    diag(1) = 2.0_dp
                    upper(1) = 1.0_dp
                    rhs(1) = 3.0_dp*delta(1) - 0.5_dp*h(1)*left_d
                    if (left_bc == 'natural') rhs(1) = 3.0_dp*delta(1)
                case ('clamped', 'first')
                    diag(1) = 1.0_dp
                    rhs(1) = left_d
                end select
                select case (right_bc)
                case ('not-a-knot')
                    h1 = h(n-2)
                    h2 = h(n-1)
                    total = h1+h2
                    lower(n) = total
                    diag(n) = h1
                    rhs(n) = (h2*h2*delta(n-2)+(2.0_dp*total+h2)*h1*delta(n-1))/total
                case ('natural', 'second')
                    lower(n) = 1.0_dp
                    diag(n) = 2.0_dp
                    rhs(n) = 3.0_dp*delta(n-1) + 0.5_dp*h(n-1)*right_d
                    if (right_bc == 'natural') rhs(n) = 3.0_dp*delta(n-1)
                case ('clamped', 'first')
                    diag(n) = 1.0_dp
                    rhs(n) = right_d
                end select
                call solve_cubic_tridiagonal(lower, diag, upper, rhs, slopes, code)
                if (code /= interpolate_status_success) goto 900
            end if
        end if
        call make_piecewise(poly, x, y, h, delta, slopes, code)
        if (code == interpolate_status_success) then
            poly%allow_extrapolation = .true.
            poly%periodic_extension = periodic
            if (present(extrapolate)) then
                poly%allow_extrapolation = extrapolate
                poly%periodic_extension = .false.
            end if
        end if
900     if (present(status)) status = code
    end function cubic_spline_1d

    ! Matrix ordinate overload: y(i,k) is sample i of output channel k.
    ! This matches scipy.interpolate.CubicSpline(x, y, axis=0).
    function cubic_spline_matrix(x, y, bc_type, status, extrapolate, bc_left, bc_right, &
        left_values, right_values) result(poly)
        real(dp), intent(in) :: x(:) !! strictly increasing knot positions, length n >= 2
        real(dp), intent(in) :: y(:,:) !! shape (n, m): m independent real-valued data series
        character(len=*), intent(in), optional :: bc_type !! shared endpoint condition name, default not-a-knot
        integer, intent(out), optional :: status !! zero on success; nonzero if any series fails
        logical, intent(in), optional :: extrapolate !! shared extrapolation setting, periodic default when periodic
        character(len=*), intent(in), optional :: bc_left !! override shared left boundary condition
        character(len=*), intent(in), optional :: bc_right !! override shared right boundary condition
        real(dp), intent(in), optional :: left_values(:) !! length m left boundary derivative values (if applicable)
        real(dp), intent(in), optional :: right_values(:) !! length m right boundary derivative values (if applicable)
        type(piecewise_cubic_array_interpolator) :: poly
        real(dp) :: lv, rv
        integer :: j, m, code
        logical :: has_values

        code = interpolate_status_invalid_input
        m = size(y,2)
        if (size(y,1) /= size(x) .or. m < 1) goto 900
        if (present(left_values)) then
            if (size(left_values) /= m) goto 900
        end if
        if (present(right_values)) then
            if (size(right_values) /= m) goto 900
        end if
        has_values = present(left_values) .or. present(right_values)
        allocate(poly%columns(m))
        do j = 1,m
            lv = 0.0_dp
            rv = 0.0_dp
            if (present(left_values)) lv = left_values(j)
            if (present(right_values)) rv = right_values(j)
            if (has_values) then
                poly%columns(j) = cubic_spline_1d(x,y(:,j),bc_type=bc_type,status=code, &
                    extrapolate=extrapolate,bc_left=bc_left,bc_right=bc_right, &
                    left_value=lv,right_value=rv)
            else
                poly%columns(j) = cubic_spline_1d(x,y(:,j),bc_type=bc_type,status=code, &
                    extrapolate=extrapolate,bc_left=bc_left,bc_right=bc_right)
            end if
            if (code /= interpolate_status_success) then
                deallocate(poly%columns)
                goto 900
            end if
        end do
900     if (present(status)) status = code
    end function cubic_spline_matrix

    pure logical function valid_matrix(self) result(ok)
        class(piecewise_cubic_array_interpolator), intent(in) :: self !! potentially uninitialized matrix spline
        integer :: j
        ok = allocated(self%columns)
        if (.not. ok) return
        ok = size(self%columns) > 0
        if (.not. ok) return
        do j = 1,size(self%columns)
            if (.not. self%columns(j)%valid()) ok = .false.
        end do
    end function valid_matrix

    pure function evaluate_matrix_scalar(self, z, extrapolate) result(vals)
        class(piecewise_cubic_array_interpolator), intent(in) :: self !! fitted multichannel spline
        real(dp), intent(in) :: z !! single query coordinate
        logical, intent(in), optional :: extrapolate !! optional extension override
        real(dp), allocatable :: vals(:)
        integer :: j

        if (.not. self%valid()) then
            allocate(vals(0))
            return
        end if
        allocate(vals(size(self%columns)))
        do j = 1,size(self%columns)
            vals(j) = self%columns(j)%evaluate(z,extrapolate)
        end do
    end function evaluate_matrix_scalar

    pure function evaluate_matrix_array(self, z, extrapolate) result(vals)
        class(piecewise_cubic_array_interpolator), intent(in) :: self !! fitted multichannel spline
        real(dp), intent(in) :: z(:) !! query coordinates, first dimension of result
        logical, intent(in), optional :: extrapolate !! extension override
        real(dp), allocatable :: vals(:,:)
        integer :: j

        if (.not. self%valid()) then
            allocate(vals(size(z),0))
            return
        end if
        allocate(vals(size(z),size(self%columns)))
        do j = 1,size(self%columns)
            vals(:,j) = self%columns(j)%evaluate(z,extrapolate)
        end do
    end function evaluate_matrix_array

    pure function derivative_matrix_scalar(self, z, order, extrapolate) result(vals)
        class(piecewise_cubic_array_interpolator), intent(in) :: self !! fitted multichannel spline
        real(dp), intent(in) :: z !! query coordinate
        integer, intent(in), optional :: order !! derivative order, default 1
        logical, intent(in), optional :: extrapolate !! extension override
        real(dp), allocatable :: vals(:)
        integer :: j

        if (.not. self%valid()) then
            allocate(vals(0))
            return
        end if
        allocate(vals(size(self%columns)))
        do j = 1,size(self%columns)
            vals(j) = self%columns(j)%derivative(z,order,extrapolate)
        end do
    end function derivative_matrix_scalar

    pure function derivative_matrix_array(self, z, order, extrapolate) result(vals)
        class(piecewise_cubic_array_interpolator), intent(in) :: self !! fitted multichannel spline
        real(dp), intent(in) :: z(:) !! first dimension of result
        integer, intent(in), optional :: order !! derivative order, default 1
        logical, intent(in), optional :: extrapolate !! extension override
        real(dp), allocatable :: vals(:,:)
        integer :: j

        if (.not. self%valid()) then
            allocate(vals(size(z),0))
            return
        end if
        allocate(vals(size(z),size(self%columns)))
        do j = 1,size(self%columns)
            vals(:,j) = self%columns(j)%derivative(z,order,extrapolate)
        end do
    end function derivative_matrix_array

    pure function integrate_matrix(self, a, b, extrapolate) result(vals)
        class(piecewise_cubic_array_interpolator), intent(in) :: self !! fitted multichannel spline
        real(dp), intent(in) :: a !! integration starting point
        real(dp), intent(in) :: b !! integration ending point
        logical, intent(in), optional :: extrapolate !! extension override
        real(dp), allocatable :: vals(:)
        integer :: j

        if (.not. self%valid()) then
            allocate(vals(0))
            return
        end if
        allocate(vals(size(self%columns)))
        do j = 1,size(self%columns)
            vals(j) = self%columns(j)%integrate(a,b,extrapolate)
        end do
    end function integrate_matrix

    pure function get_matrix_coefficients(self) result(c)
        class(piecewise_cubic_array_interpolator), intent(in) :: self !! fitted multichannel spline
        real(dp), allocatable :: c(:,:,:) !! ascending powers, intervals, channels
        integer :: j, n, m

        if (.not. self%valid()) then
            allocate(c(4,0,0))
            return
        end if
        m = size(self%columns)
        n = size(self%columns(1)%x)-1
        allocate(c(4,n,m))
        do j = 1,m
            c(:,:,j) = self%columns(j)%coefficients()
        end do
    end function get_matrix_coefficients

    pure logical function valid_cubic_boundary(bc) result(ok)
        character(len=*), intent(in) :: bc !! normalized cubic boundary-condition name
        select case (bc)
        case ('not-a-knot','natural','clamped','first','second','periodic')
            ok = .true.
        case default
            ok = .false.
        end select
    end function valid_cubic_boundary

    subroutine solve_cubic_tridiagonal(lower, diag, upper, rhs, sol, status)
        real(dp), intent(in) :: lower(:) !! subdiagonal coefficients; lower(1) unused
        real(dp), intent(in) :: diag(:) !! diagonal coefficients
        real(dp), intent(in) :: upper(:) !! superdiagonal coefficients; upper(n) unused
        real(dp), intent(in) :: rhs(:) !! right-hand side
        real(dp), intent(out) :: sol(:) !! tridiagonal linear-system solution
        integer, intent(out) :: status !! zero on successful finite solve
        real(dp), allocatable :: b(:), v(:)
        real(dp) :: factor
        integer :: n, i

        status = interpolate_status_invalid_input
        n = size(rhs)
        if (n < 1 .or. size(lower) /= n .or. size(diag) /= n .or. &
            size(upper) /= n .or. size(sol) /= n) return
        if (any(.not. ieee_is_finite(diag)) .or. any(.not. ieee_is_finite(rhs))) return
        allocate(b(n),v(n))
        b = diag
        v = rhs
        if (.not. (abs(b(1)) > 0.0_dp)) return
        do i = 2,n
            factor = lower(i)/b(i-1)
            b(i) = b(i)-factor*upper(i-1)
            v(i) = v(i)-factor*v(i-1)
            if (.not. (abs(b(i)) > 0.0_dp)) return
        end do
        sol(n) = v(n)/b(n)
        do i = n-1,1,-1
            sol(i) = (v(i)-upper(i)*sol(i+1))/b(i)
        end do
        if (any(.not. ieee_is_finite(sol))) return
        status = interpolate_status_success
    end subroutine solve_cubic_tridiagonal

    subroutine periodic_cubic_slopes(h, delta, slopes, status)
        real(dp), intent(in) :: h(:) !! validated positive interval widths
        real(dp), intent(in) :: delta(:) !! validated interval secants with periodic endpoint values
        real(dp), intent(out) :: slopes(:) !! resulting first derivative at every knot
        integer, intent(out) :: status !! zero for valid periodic spline
        real(dp), allocatable :: lower(:), diag(:), upper(:), rhs(:), &
            work_diag(:), u(:), v(:), moment(:), all_moments(:)
        real(dp) :: length, gamma, alpha, beta, denom, factor
        integer :: n, m, i, code

        status = interpolate_status_invalid_input
        n = size(slopes)
        m = n-1
        if (n == 2) then
            slopes = 0.0_dp
            status = interpolate_status_success
            return
        end if
        allocate(all_moments(n))
        if (n == 3) then
            length = h(1)+h(2)
            all_moments(1) = 6.0_dp*(delta(1)-delta(2))/length
            all_moments(2) = -all_moments(1)
            all_moments(3) = all_moments(1)
        else
            ! Symmetric cyclic tridiagonal system on the unique knots.
            allocate(lower(m),diag(m),upper(m),rhs(m),work_diag(m), &
                u(m),v(m),moment(m))
            lower = 0.0_dp
            upper = 0.0_dp
            diag(1) = 2.0_dp*(h(1)+h(m))
            upper(1) = h(1)
            rhs(1) = 6.0_dp*(delta(1)-delta(m))
            do i = 2,m
                lower(i) = h(i-1)
                diag(i) = 2.0_dp*(h(i-1)+h(i))
                rhs(i) = 6.0_dp*(delta(i)-delta(i-1))
                if (i < m) upper(i) = h(i)
            end do
            alpha = h(m)
            beta = h(m)
            gamma = -diag(1)
            work_diag = diag
            work_diag(1) = diag(1)-gamma
            work_diag(m) = diag(m)-alpha*beta/gamma
            call solve_cubic_tridiagonal(lower,work_diag,upper,rhs,moment,code)
            if (code /= interpolate_status_success) return
            u = 0.0_dp
            u(1) = gamma
            u(m) = beta
            call solve_cubic_tridiagonal(lower,work_diag,upper,u,v,code)
            if (code /= interpolate_status_success) return
            denom = 1.0_dp+v(1)+alpha*v(m)/gamma
            if (.not. ieee_is_finite(denom)) return
            if (.not. (abs(denom) > 0.0_dp)) return
            factor = (moment(1)+alpha*moment(m)/gamma)/denom
            all_moments(1:m) = moment-factor*v
            all_moments(n) = all_moments(1)
        end if
        slopes(1) = delta(1)-h(1)*(2.0_dp*all_moments(1)+all_moments(2))/6.0_dp
        do i = 2,n-1
            slopes(i) = delta(i-1)+h(i-1)*(all_moments(i-1)+2.0_dp*all_moments(i))/6.0_dp
        end do
        slopes(n) = slopes(1)
        if (any(.not. ieee_is_finite(slopes))) return
        status = interpolate_status_success
    end subroutine periodic_cubic_slopes

    function cubic_hermite_spline(x, y, dydx, status, extrapolate) result(poly)
        real(dp), intent(in) :: x(:) !! strictly increasing finite knot coordinates, length >= 2
        real(dp), intent(in) :: y(:) !! finite ordinates at the knots
        real(dp), intent(in) :: dydx(:) !! caller-provided first derivatives at the knots
        integer, intent(out), optional :: status !! 0 on success, 1 on invalid input
        logical, intent(in), optional :: extrapolate !! extrapolate beyond knot range, default true
        type(piecewise_cubic_interpolator) :: poly
        real(dp), allocatable :: h(:), delta(:)
        integer :: code

        call validate_samples(x, y, h, delta, code)
        if (code == interpolate_status_success) then
            if (size(dydx) /= size(x)) then
                code = interpolate_status_invalid_input
            else if (any(.not. ieee_is_finite(dydx))) then
                code = interpolate_status_invalid_input
            else
                call make_piecewise(poly, x, y, h, delta, dydx, code)
                poly%allow_extrapolation = .true.
                if (present(extrapolate)) poly%allow_extrapolation = extrapolate
            end if
        end if
        if (present(status)) status = code
    end function cubic_hermite_spline

    function pchip_interpolator(x, y, status, extrapolate) result(poly)
        real(dp), intent(in) :: x(:) !! strictly increasing finite knot coordinates, length >= 2
        real(dp), intent(in) :: y(:) !! finite sample ordinates, same length as x
        integer, intent(out), optional :: status !! 0 on success, 1 on invalid input
        logical, intent(in), optional :: extrapolate !! extrapolate outside knots, default true
        type(piecewise_cubic_interpolator) :: poly
        real(dp), allocatable :: h(:), delta(:), slopes(:)
        integer :: n, i, code

        call validate_samples(x, y, h, delta, code)
        if (code == interpolate_status_success) then
            n = size(x)
            allocate(slopes(n))
            if (n == 2) then
                slopes = delta(1)
            else
                slopes(1) = pchip_end_slope(h(1), h(2), delta(1), delta(2))
                slopes(n) = pchip_end_slope(h(n-1), h(n-2), delta(n-1), delta(n-2))
                do i = 2, n-1
                    if (abs(delta(i-1)) <= 0.0_dp .or. abs(delta(i)) <= 0.0_dp) then
                        slopes(i) = 0.0_dp
                    else if ((delta(i-1) > 0.0_dp) .neqv. (delta(i) > 0.0_dp)) then
                        slopes(i) = 0.0_dp
                    else
                        ! Normalize the harmonic weights before combining reciprocals.
                        slopes(i) = 1.0_dp / &
                            (((2.0_dp*h(i) + h(i-1)) / (3.0_dp*(h(i-1)+h(i)))) / delta(i-1) + &
                            ((h(i) + 2.0_dp*h(i-1)) / (3.0_dp*(h(i-1)+h(i)))) / delta(i))
                    end if
                end do
            end if
            call make_piecewise(poly, x, y, h, delta, slopes, code)
            poly%allow_extrapolation = .true.
            if (present(extrapolate)) poly%allow_extrapolation = extrapolate
        end if
        if (present(status)) status = code
    end function pchip_interpolator

    function akima1d_interpolator(x, y, method, status, extrapolate) result(poly)
        real(dp), intent(in) :: x(:) !! strictly increasing finite knot coordinates, length >= 2
        real(dp), intent(in) :: y(:) !! finite sample ordinates at each knot
        character(len=*), intent(in), optional :: method !! 'akima' (default) or 'makima'
        integer, intent(out), optional :: status !! 0 on success, 1 on invalid input
        logical, intent(in), optional :: extrapolate !! extrapolate outside knots, default false
        type(piecewise_cubic_interpolator) :: poly
        real(dp), allocatable :: h(:), delta(:), slopes(:), m(:)
        real(dp) :: w1, w2, total, threshold
        logical :: modified
        integer :: n, i, code

        code = interpolate_status_invalid_input
        modified = .false.
        if (present(method)) then
            select case (trim(method))
            case ('akima')
                modified = .false.
            case ('makima')
                modified = .true.
            case default
                if (present(status)) status = code
                return
            end select
        end if
        call validate_samples(x, y, h, delta, code)
        if (code /= interpolate_status_success) then
            if (present(status)) status = code
            return
        end if
        n = size(x)
        allocate(slopes(n))
        if (n == 2) then
            ! Akima with only two samples reduces to the unique straight line.
            slopes = delta(1)
        else
            allocate(m(n+3))
            m(3:n+1) = delta
            m(2) = 2.0_dp*m(3) - m(4)
            m(1) = 2.0_dp*m(2) - m(3)
            m(n+2) = 2.0_dp*m(n+1) - m(n)
            m(n+3) = 2.0_dp*m(n+2) - m(n+1)
            ! In the ambiguous case m(i+1) != m(i+2) with zero weights,
            ! the reference default is the mean of the outer two slopes.
            slopes = 0.5_dp*(m(1:n)+m(4:n+3))
            threshold = 0.0_dp
            do i = 1, n
                w1 = abs(m(i+3)-m(i+2))
                w2 = abs(m(i+1)-m(i))
                if (modified) then
                    w1 = w1 + 0.5_dp*abs(m(i+3)+m(i+2))
                    w2 = w2 + 0.5_dp*abs(m(i+1)+m(i))
                end if
                threshold = max(threshold, w1+w2)
            end do
            threshold = 1.0e-9_dp*threshold
            do i = 1, n
                w1 = abs(m(i+3)-m(i+2))
                w2 = abs(m(i+1)-m(i))
                if (modified) then
                    w1 = w1 + 0.5_dp*abs(m(i+3)+m(i+2))
                    w2 = w2 + 0.5_dp*abs(m(i+1)+m(i))
                end if
                total = w1+w2
                if (total > threshold .and. ieee_is_finite(total)) then
                    slopes(i) = m(i+1) + (w2/total)*(m(i+2)-m(i+1))
                end if
            end do
        end if
        call make_piecewise(poly, x, y, h, delta, slopes, code)
        poly%allow_extrapolation = .false.
        if (present(extrapolate)) poly%allow_extrapolation = extrapolate
        if (present(status)) status = code
    end function akima1d_interpolator

    subroutine validate_samples(x, y, h, delta, status)
        real(dp), intent(in) :: x(:) !! candidate ordered abscissas
        real(dp), intent(in) :: y(:) !! candidate ordinate samples
        real(dp), allocatable, intent(out) :: h(:) !! positive interval widths on success
        real(dp), allocatable, intent(out) :: delta(:) !! secant slopes on success
        integer, intent(out) :: status !! 0 for valid finite data
        integer :: n, i

        status = interpolate_status_invalid_input
        n = size(x)
        if (n < 2 .or. size(y) /= n) return
        if (any(.not. ieee_is_finite(x)) .or. any(.not. ieee_is_finite(y))) return
        allocate(h(n-1), delta(n-1))
        do i = 1, n-1
            if (x(i+1) <= x(i)) return
            h(i) = x(i+1)-x(i)
            if (.not. ieee_is_finite(h(i)) .or. h(i) <= 0.0_dp) return
            delta(i) = (y(i+1)-y(i))/h(i)
        end do
        if (any(.not. ieee_is_finite(delta))) return
        status = interpolate_status_success
    end subroutine validate_samples

    pure function pchip_end_slope(h1, h2, s1, s2) result(d)
        real(dp), intent(in) :: h1 !! adjacent interval width at endpoint
        real(dp), intent(in) :: h2 !! next interior interval width
        real(dp), intent(in) :: s1 !! adjacent secant slope
        real(dp), intent(in) :: s2 !! next interior secant slope
        real(dp) :: d

        d = ((2.0_dp*h1+h2)*s1 - h1*s2)/(h1+h2)
        if (abs(s1) <= 0.0_dp .or. (d > 0.0_dp) .neqv. (s1 > 0.0_dp)) then
            d = 0.0_dp
        else if ((s1 > 0.0_dp) .neqv. (s2 > 0.0_dp)) then
            if (abs(d) > 3.0_dp*abs(s1)) d = 3.0_dp*s1
        end if
    end function pchip_end_slope

    subroutine make_piecewise(poly, x, y, h, delta, slopes, status)
        type(piecewise_cubic_interpolator), intent(inout) :: poly !! newly fitted result
        real(dp), intent(in) :: x(:) !! validated strictly increasing knot positions
        real(dp), intent(in) :: y(:) !! validated ordinates
        real(dp), intent(in) :: h(:) !! interval widths
        real(dp), intent(in) :: delta(:) !! secant slopes
        real(dp), intent(in) :: slopes(:) !! interpolated first derivatives
        integer, intent(out) :: status !! 0 on success, else invalid numerical values
        integer :: j, n

        status = interpolate_status_invalid_input
        if (any(.not. ieee_is_finite(slopes))) return
        n = size(x)
        allocate(poly%x(n), poly%c(4,n-1), poly%area(n))
        poly%x = x
        poly%area(1) = 0.0_dp
        do j = 1, n-1
            poly%c(1,j) = y(j)
            poly%c(2,j) = slopes(j)
            poly%c(3,j) = (3.0_dp*delta(j) - 2.0_dp*slopes(j) - slopes(j+1))/h(j)
            poly%c(4,j) = (slopes(j) + slopes(j+1) - 2.0_dp*delta(j))/(h(j)*h(j))
            poly%area(j+1) = poly%area(j) + primitive_segment(poly%c(:,j), h(j))
        end do
        if (any(.not. ieee_is_finite(poly%c)) .or. &
            any(.not. ieee_is_finite(poly%area))) then
            deallocate(poly%x, poly%c, poly%area)
            return
        end if
        status = interpolate_status_success
    end subroutine make_piecewise

    pure function primitive_segment(c, t) result(val)
        real(dp), intent(in) :: c(4) !! Hermite polynomial coefficients in ascending order
        real(dp), intent(in) :: t !! coordinate relative to left knot
        real(dp) :: val
        val = t*(c(1) + t*(c(2)/2.0_dp + t*(c(3)/3.0_dp + t*c(4)/4.0_dp)))
    end function primitive_segment

    pure integer function locate_interval(x, z) result(j)
        real(dp), intent(in) :: x(:) !! increasing fitted knot positions
        real(dp), intent(in) :: z !! finite query point
        integer :: lo, hi, mid

        lo = 1
        hi = size(x)
        if (z <= x(1)) then
            j = 1
            return
        else if (z >= x(hi)) then
            j = hi-1
            return
        end if
        do while (hi-lo > 1)
            mid = lo + (hi-lo)/2
            if (z < x(mid)) then
                hi = mid
            else
                lo = mid
            end if
        end do
        j = lo
    end function locate_interval

    pure function evaluate_scalar(self, z, extrapolate) result(val)
        class(piecewise_cubic_interpolator), intent(in) :: self !! fitted interpolator
        real(dp), intent(in) :: z !! coordinate at which to evaluate the cubic
        logical, intent(in), optional :: extrapolate !! override the fitted extension policy
        real(dp) :: val

        val = self%derivative_scalar(z, 0, extrapolate)
    end function evaluate_scalar

    pure function evaluate_array(self, z, extrapolate) result(vals)
        class(piecewise_cubic_interpolator), intent(in) :: self !! fitted interpolator
        real(dp), intent(in) :: z(:) !! vector of query coordinates
        logical, intent(in), optional :: extrapolate !! override extension policy
        real(dp) :: vals(size(z))
        integer :: i

        do i = 1, size(z)
            vals(i) = self%evaluate_scalar(z(i), extrapolate)
        end do
    end function evaluate_array

    pure function derivative_scalar(self, z, order, extrapolate) result(val)
        class(piecewise_cubic_interpolator), intent(in) :: self !! fitted interpolator
        real(dp), intent(in) :: z !! coordinate at which to evaluate a derivative
        integer, intent(in), optional :: order !! nonnegative derivative order, default 1
        logical, intent(in), optional :: extrapolate !! override extension policy
        real(dp) :: val, t, actual_z, period
        integer :: j, nu
        logical :: extend

        val = ieee_value(0.0_dp, ieee_quiet_nan)
        if (.not. self%valid()) return
        if (.not. ieee_is_finite(z)) return
        nu = 1
        if (present(order)) nu = order
        if (nu < 0) return
        extend = self%allow_extrapolation
        if (present(extrapolate)) extend = extrapolate
        if (.not. extend .and. (z < self%x(1) .or. z > self%x(size(self%x)))) return
        actual_z = z
        if (self%periodic_extension .and. .not. present(extrapolate)) then
            period = self%x(size(self%x))-self%x(1)
            actual_z = self%x(1)+modulo(z-self%x(1),period)
        end if
        j = locate_interval(self%x, actual_z)
        t = actual_z - self%x(j)
        select case (nu)
        case (0)
            val = self%c(1,j) + t*(self%c(2,j) + t*(self%c(3,j) + t*self%c(4,j)))
        case (1)
            val = self%c(2,j) + t*(2.0_dp*self%c(3,j) + t*3.0_dp*self%c(4,j))
        case (2)
            val = 2.0_dp*self%c(3,j) + t*6.0_dp*self%c(4,j)
        case (3)
            val = 6.0_dp*self%c(4,j)
        case default
            val = 0.0_dp
        end select
    end function derivative_scalar

    pure function derivative_array(self, z, order, extrapolate) result(vals)
        class(piecewise_cubic_interpolator), intent(in) :: self !! fitted interpolator
        real(dp), intent(in) :: z(:) !! coordinates at which to evaluate derivatives
        integer, intent(in), optional :: order !! nonnegative derivative order, default 1
        logical, intent(in), optional :: extrapolate !! override extension policy
        real(dp) :: vals(size(z))
        integer :: i

        do i = 1, size(z)
            vals(i) = self%derivative_scalar(z(i), order, extrapolate)
        end do
    end function derivative_array

    pure function integrate_piecewise(self, a, b, extrapolate) result(val)
        class(piecewise_cubic_interpolator), intent(in) :: self !! fitted interpolator
        real(dp), intent(in) :: a !! left integration endpoint; may exceed b
        real(dp), intent(in) :: b !! right integration endpoint
        logical, intent(in), optional :: extrapolate !! override extension policy
        real(dp) :: val
        logical :: extend

        val = ieee_value(0.0_dp, ieee_quiet_nan)
        if (.not. self%valid()) return
        if (.not. ieee_is_finite(a) .or. .not. ieee_is_finite(b)) return
        extend = self%allow_extrapolation
        if (present(extrapolate)) extend = extrapolate
        if (.not. extend) then
            if (min(a,b) < self%x(1) .or. max(a,b) > self%x(size(self%x))) return
        end if
        if (self%periodic_extension .and. .not. present(extrapolate)) then
            val = periodic_cumulative_integral(self, b) - &
                periodic_cumulative_integral(self, a)
        else
            val = cumulative_integral(self, b) - cumulative_integral(self, a)
        end if
    end function integrate_piecewise

    pure function periodic_cumulative_integral(self, z) result(val)
        class(piecewise_cubic_interpolator), intent(in) :: self !! periodic fitted interpolator
        real(dp), intent(in) :: z !! finite integration endpoint, possibly outside base period
        real(dp) :: val, period, r, cycles
        integer :: n

        n = size(self%x)
        period = self%x(n)-self%x(1)
        r = modulo(z-self%x(1), period)
        cycles = anint(((z-self%x(1))-r)/period)
        val = cycles*self%area(n)+cumulative_integral(self,self%x(1)+r)
    end function periodic_cumulative_integral

    pure function cumulative_integral(self, z) result(val)
        class(piecewise_cubic_interpolator), intent(in) :: self !! valid interpolator
        real(dp), intent(in) :: z !! finite position; extrapolation allowed
        real(dp) :: val
        integer :: j

        j = locate_interval(self%x, z)
        val = self%area(j) + primitive_segment(self%c(:,j), z-self%x(j))
    end function cumulative_integral

    pure logical function valid_piecewise(self) result(ok)
        class(piecewise_cubic_interpolator), intent(in) :: self !! potentially uninitialized object
        ok = allocated(self%x) .and. allocated(self%c) .and. allocated(self%area)
    end function valid_piecewise

    pure function get_piecewise_knots(self) result(x)
        class(piecewise_cubic_interpolator), intent(in) :: self !! fitted interpolator
        real(dp), allocatable :: x(:)
        if (allocated(self%x)) then
            x = self%x
        else
            allocate(x(0))
        end if
    end function get_piecewise_knots

    pure function get_piecewise_coefficients(self) result(c)
        class(piecewise_cubic_interpolator), intent(in) :: self !! fitted interpolator
        real(dp), allocatable :: c(:,:)
        ! In contrast to scipy.interpolate.PPoly.c, this array is ascending
        ! powers: row 1 is constant, row 4 is cubic.
        if (allocated(self%c)) then
            c = self%c
        else
            allocate(c(4,0))
        end if
    end function get_piecewise_coefficients

    function pchip_interpolate_scalar(x, y, z, order, status, extrapolate) result(val)
        real(dp), intent(in) :: x(:) !! strictly increasing interpolation knots
        real(dp), intent(in) :: y(:) !! sampled ordinates
        real(dp), intent(in) :: z !! query coordinate
        integer, intent(in), optional :: order !! derivative order, default zero
        integer, intent(out), optional :: status !! fitting status, zero on success
        logical, intent(in), optional :: extrapolate !! default true
        real(dp) :: val
        type(piecewise_cubic_interpolator) :: poly
        integer :: code, nu

        poly = pchip_interpolator(x, y, code, extrapolate)
        nu = 0
        if (present(order)) nu = order
        val = poly%derivative(z, nu)
        if (present(status)) status = code
    end function pchip_interpolate_scalar

    function pchip_interpolate_array(x, y, z, order, status, extrapolate) result(vals)
        real(dp), intent(in) :: x(:) !! strictly increasing interpolation knots
        real(dp), intent(in) :: y(:) !! sampled ordinates
        real(dp), intent(in) :: z(:) !! vector of query coordinates
        integer, intent(in), optional :: order !! derivative order, default zero
        integer, intent(out), optional :: status !! fitting status, zero on success
        logical, intent(in), optional :: extrapolate !! default true
        real(dp) :: vals(size(z))
        type(piecewise_cubic_interpolator) :: poly
        integer :: code, nu

        poly = pchip_interpolator(x, y, code, extrapolate)
        nu = 0
        if (present(order)) nu = order
        vals = poly%derivative(z, nu)
        if (present(status)) status = code
    end function pchip_interpolate_array

end module scifort_interpolate
