! SPDX-License-Identifier: MIT AND BSD-3-Clause
! Copyright (c) 2026 SciFort contributors
! Copyright (c) 2001-2002 Enthought, Inc. 2003, SciPy Developers.
!
! Adaptive multidimensional cubature over hyperrectangles. The Genz-Malik
! degree-7/degree-5 embedded rule and Gauss-Kronrod 15/21 constants are adapted
! from SciPy 1.17.0 scipy/integrate/_rules, under SciPy's BSD-3-Clause license.
! The adaptive storage, Fortran callbacks, point splitting, and corrected
! infinite-limit transformations are SciFort code. See CODE_PROVENANCE.md.

module scifort_integrate_cubature
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_integrate_quadrature, only : integrate_status_invalid_argument, &
        integrate_status_invalid_shape, integrate_status_numerical_failure, &
        integrate_status_success
    use scifort_kinds, only : dp
    use scifort_math, only : quiet_nan
    use scifort_orthogonal_roots, only : orthogonal_roots_status_success, roots_legendre
    implicit none
    private

    integer, parameter, public :: cubature_status_max_subdivisions = 4

    type :: cubature_region
        real(dp), allocatable :: lower(:)
        real(dp), allocatable :: upper(:)
        real(dp), allocatable :: estimate(:)
        real(dp), allocatable :: error(:)
        real(dp) :: score = 0.0_dp
    end type cubature_region

    abstract interface
        function cubature_scalar_integrand(x) result(value)
            import :: dp
            implicit none
            real(dp), intent(in) :: x(:) !! one multidimensional integration point
            real(dp) :: value
        end function cubature_scalar_integrand

        subroutine cubature_vector_integrand(x, values)
            import :: dp
            implicit none
            real(dp), intent(in) :: x(:) !! one multidimensional integration point
            real(dp), intent(out) :: values(:) !! integrand output at x
        end subroutine cubature_vector_integrand
    end interface

    interface cubature
        module procedure cubature_scalar
        module procedure cubature_vector
    end interface cubature

    public :: cubature

    real(dp), parameter :: gk21_nodes(21) = [ &
        0.995657163025808080735527280689003_dp, &
        0.973906528517171720077964012084452_dp, &
        0.930157491355708226001207180059508_dp, &
        0.865063366688984510732096688423493_dp, &
        0.780817726586416897063717578345042_dp, &
        0.679409568299024406234327365114874_dp, &
        0.562757134668604683339000099272694_dp, &
        0.433395394129247190799265943165784_dp, &
        0.294392862701460198131126603103866_dp, &
        0.148874338981631210884826001129720_dp, &
        0.0_dp, &
        -0.148874338981631210884826001129720_dp, &
        -0.294392862701460198131126603103866_dp, &
        -0.433395394129247190799265943165784_dp, &
        -0.562757134668604683339000099272694_dp, &
        -0.679409568299024406234327365114874_dp, &
        -0.780817726586416897063717578345042_dp, &
        -0.865063366688984510732096688423493_dp, &
        -0.930157491355708226001207180059508_dp, &
        -0.973906528517171720077964012084452_dp, &
        -0.995657163025808080735527280689003_dp]

    real(dp), parameter :: gk21_weights(21) = [ &
        0.011694638867371874278064396062192_dp, &
        0.032558162307964727478818972459390_dp, &
        0.054755896574351996031381300244580_dp, &
        0.075039674810919952767043140916190_dp, &
        0.093125454583697605535065465083366_dp, &
        0.109387158802297641899210590325805_dp, &
        0.123491976262065851077958109831074_dp, &
        0.134709217311473325928054001771707_dp, &
        0.142775938577060080797094273138717_dp, &
        0.147739104901338491374841515972068_dp, &
        0.149445554002916905664936468389821_dp, &
        0.147739104901338491374841515972068_dp, &
        0.142775938577060080797094273138717_dp, &
        0.134709217311473325928054001771707_dp, &
        0.123491976262065851077958109831074_dp, &
        0.109387158802297641899210590325805_dp, &
        0.093125454583697605535065465083366_dp, &
        0.075039674810919952767043140916190_dp, &
        0.054755896574351996031381300244580_dp, &
        0.032558162307964727478818972459390_dp, &
        0.011694638867371874278064396062192_dp]

    real(dp), parameter :: gk15_nodes(15) = [ &
        0.991455371120812639206854697526329_dp, &
        0.949107912342758524526189684047851_dp, &
        0.864864423359769072789712788640926_dp, &
        0.741531185599394439863864773280788_dp, &
        0.586087235467691130294144838258730_dp, &
        0.405845151377397166906606412076961_dp, &
        0.207784955007898467600689403773245_dp, &
        0.0_dp, &
        -0.207784955007898467600689403773245_dp, &
        -0.405845151377397166906606412076961_dp, &
        -0.586087235467691130294144838258730_dp, &
        -0.741531185599394439863864773280788_dp, &
        -0.864864423359769072789712788640926_dp, &
        -0.949107912342758524526189684047851_dp, &
        -0.991455371120812639206854697526329_dp]

    real(dp), parameter :: gk15_weights(15) = [ &
        0.022935322010529224963732008058970_dp, &
        0.063092092629978553290700663189204_dp, &
        0.104790010322250183839876322541518_dp, &
        0.140653259715525918745189590510238_dp, &
        0.169004726639267902826583426598550_dp, &
        0.190350578064785409913256402421014_dp, &
        0.204432940075298892414161999234649_dp, &
        0.209482141084727828012999174891714_dp, &
        0.204432940075298892414161999234649_dp, &
        0.190350578064785409913256402421014_dp, &
        0.169004726639267902826583426598550_dp, &
        0.140653259715525918745189590510238_dp, &
        0.104790010322250183839876322541518_dp, &
        0.063092092629978553290700663189204_dp, &
        0.022935322010529224963732008058970_dp]

contains

    subroutine cubature_scalar(f, a, b, result, abserr, rule, rtol, atol, &
            max_subdivisions, points, neval, subdivisions, status)
        procedure(cubature_scalar_integrand) :: f !! scalar integrand callback
        real(dp), intent(in) :: a(:) !! lower bounds; components may be infinite or reversed
        real(dp), intent(in) :: b(:) !! upper bounds; components may be infinite or reversed
        real(dp), intent(out) :: result !! cubature estimate
        real(dp), intent(out) :: abserr !! absolute-error estimate
        character(len=*), intent(in), optional :: rule !! gk21, gk15, gauss-kronrod, or genz-malik
        real(dp), intent(in), optional :: rtol !! relative tolerance; default 1e-8
        real(dp), intent(in), optional :: atol !! absolute tolerance; default 0
        integer, intent(in), optional :: max_subdivisions !! adaptive refinements; default 10000
        real(dp), intent(in), optional :: points(:, :) !! points(dim, npoint) not to evaluate
        integer, intent(out), optional :: neval !! total callback evaluations
        integer, intent(out), optional :: subdivisions !! number of adaptive refinements
        integer, intent(out), optional :: status !! integrate_status_* or cubature_status_* result

        real(dp) :: result_vector(1)
        real(dp) :: error_vector(1)

        result = quiet_nan(0.0_dp)
        abserr = quiet_nan(0.0_dp)
        if (present(points)) then
            call cubature_core(vector_callback, a, b, result_vector, error_vector, &
                rule, rtol, atol, max_subdivisions, points, neval, subdivisions, status)
        else
            call cubature_core(vector_callback, a, b, result_vector, error_vector, &
                rule, rtol, atol, max_subdivisions, neval=neval, &
                subdivisions=subdivisions, status=status)
        end if
        result = result_vector(1)
        abserr = error_vector(1)

    contains

        subroutine vector_callback(x, values)
            real(dp), intent(in) :: x(:) !! one multidimensional integration point
            real(dp), intent(out) :: values(:) !! one-element output buffer

            values(1) = f(x)
        end subroutine vector_callback

    end subroutine cubature_scalar

    subroutine cubature_vector(f, a, b, result, abserr, rule, rtol, atol, &
            max_subdivisions, points, neval, subdivisions, status)
        procedure(cubature_vector_integrand) :: f !! vector-valued integrand callback
        real(dp), intent(in) :: a(:) !! lower bounds; components may be infinite or reversed
        real(dp), intent(in) :: b(:) !! upper bounds; components may be infinite or reversed
        real(dp), intent(out) :: result(:) !! cubature estimate for each output component
        real(dp), intent(out) :: abserr(:) !! absolute-error estimate for each output component
        character(len=*), intent(in), optional :: rule !! gk21, gk15, gauss-kronrod, or genz-malik
        real(dp), intent(in), optional :: rtol !! relative tolerance; default 1e-8
        real(dp), intent(in), optional :: atol !! absolute tolerance; default 0
        integer, intent(in), optional :: max_subdivisions !! adaptive refinements; default 10000
        real(dp), intent(in), optional :: points(:, :) !! points(dim, npoint) not to evaluate
        integer, intent(out), optional :: neval !! total callback evaluations
        integer, intent(out), optional :: subdivisions !! number of adaptive refinements
        integer, intent(out), optional :: status !! integrate_status_* or cubature_status_* result

        if (present(points)) then
            call cubature_core(f, a, b, result, abserr, rule, rtol, atol, &
                max_subdivisions, points, neval, subdivisions, status)
        else
            call cubature_core(f, a, b, result, abserr, rule, rtol, atol, &
                max_subdivisions, neval=neval, subdivisions=subdivisions, status=status)
        end if
    end subroutine cubature_vector

    subroutine cubature_core(f, a, b, result, abserr, rule, rtol, atol, &
            max_subdivisions, points, neval, subdivisions, status)
        procedure(cubature_vector_integrand) :: f !! vector-valued integrand callback
        real(dp), intent(in) :: a(:) !! original lower bounds
        real(dp), intent(in) :: b(:) !! original upper bounds
        real(dp), intent(out) :: result(:) !! cubature estimate
        real(dp), intent(out) :: abserr(:) !! componentwise error estimate
        character(len=*), intent(in), optional :: rule !! selected embedded cubature rule
        real(dp), intent(in), optional :: rtol !! requested relative tolerance
        real(dp), intent(in), optional :: atol !! requested absolute tolerance
        integer, intent(in), optional :: max_subdivisions !! maximum adaptive refinements
        real(dp), intent(in), optional :: points(:, :) !! physical-space points to avoid
        integer, intent(out), optional :: neval !! total callback evaluations
        integer, intent(out), optional :: subdivisions !! completed adaptive refinements
        integer, intent(out), optional :: status !! termination status

        character(len=20) :: selected_rule
        integer :: d
        integer :: evaluations
        integer :: i
        integer :: max_refinements
        integer :: nregions
        integer :: refine_count
        integer :: selected
        integer :: work_status
        real(dp) :: absolute_tolerance
        real(dp) :: orientation
        real(dp) :: relative_tolerance
        real(dp), allocatable :: lower(:)
        real(dp), allocatable :: upper(:)
        real(dp), allocatable :: transformed_lower(:)
        real(dp), allocatable :: transformed_upper(:)
        type(cubature_region), allocatable :: regions(:)

        result = quiet_nan(0.0_dp)
        abserr = quiet_nan(0.0_dp)
        evaluations = 0
        refine_count = 0
        if (present(neval)) neval = 0
        if (present(subdivisions)) subdivisions = 0
        if (present(status)) status = integrate_status_invalid_argument

        d = size(a)
        if (d < 1 .or. size(b) /= d) then
            if (present(status)) status = integrate_status_invalid_shape
            return
        end if
        if (size(result) < 1 .or. size(abserr) /= size(result)) then
            if (present(status)) status = integrate_status_invalid_shape
            return
        end if
        if (any(ieee_is_nan(a)) .or. any(ieee_is_nan(b))) return

        selected_rule = 'gk21'
        relative_tolerance = 1.0e-8_dp
        absolute_tolerance = 0.0_dp
        max_refinements = 10000
        if (present(rule)) selected_rule = trim(adjustl(rule))
        if (present(rtol)) relative_tolerance = rtol
        if (present(atol)) absolute_tolerance = atol
        if (present(max_subdivisions)) max_refinements = max_subdivisions
        if (relative_tolerance < 0.0_dp .or. absolute_tolerance < 0.0_dp) return
        if (max_refinements < 0) return
        if (.not. valid_rule(selected_rule, d)) return
        if (present(points)) then
            if (size(points, 1) /= d) then
                if (present(status)) status = integrate_status_invalid_shape
                return
            end if
            if (any(ieee_is_nan(points))) return
        end if

        allocate(lower(d), upper(d), transformed_lower(d), transformed_upper(d))
        orientation = 1.0_dp
        do i = 1, d
            if (a(i) > b(i)) then
                lower(i) = b(i)
                upper(i) = a(i)
                orientation = -orientation
            else
                lower(i) = a(i)
                upper(i) = b(i)
            end if
            if (lower(i) == upper(i)) then
                if (.not. ieee_is_finite(lower(i))) return
                result = 0.0_dp
                abserr = 0.0_dp
                if (present(status)) status = integrate_status_success
                return
            end if
        end do

        call transformed_bounds(lower, upper, transformed_lower, transformed_upper, work_status)
        if (work_status /= integrate_status_success) return

        allocate(regions(1))
        call set_empty_region(regions(1), transformed_lower, transformed_upper, size(result))
        nregions = 1

        call split_infinite_surfaces(regions, nregions, lower, upper, size(result), work_status)
        if (work_status /= integrate_status_success) then
            if (present(status)) status = work_status
            return
        end if
        if (present(points)) then
            call split_user_points(regions, nregions, lower, upper, points, size(result), work_status)
            if (work_status /= integrate_status_success) then
                if (present(status)) status = work_status
                return
            end if
        end if

        result = 0.0_dp
        abserr = 0.0_dp
        do i = 1, nregions
            call evaluate_region(f, lower, upper, selected_rule, regions(i), evaluations, work_status)
            if (work_status /= integrate_status_success) then
                if (present(status)) status = work_status
                if (present(neval)) neval = evaluations
                return
            end if
            result = result + regions(i)%estimate
            abserr = abserr + regions(i)%error
        end do

        do while (.not. converged(result, abserr, absolute_tolerance, relative_tolerance))
            if (refine_count >= max_refinements) then
                if (present(status)) status = cubature_status_max_subdivisions
                result = orientation * result
                if (present(neval)) neval = evaluations
                if (present(subdivisions)) subdivisions = refine_count
                return
            end if

            selected = largest_error_region(regions, nregions)
            result = result - regions(selected)%estimate
            abserr = max(0.0_dp, abserr - regions(selected)%error)
            call refine_region(f, lower, upper, selected_rule, regions, nregions, selected, &
                size(result), evaluations, result, abserr, work_status)
            if (work_status /= integrate_status_success) then
                if (present(status)) status = work_status
                if (present(neval)) neval = evaluations
                if (present(subdivisions)) subdivisions = refine_count
                return
            end if
            refine_count = refine_count + 1
        end do

        result = orientation * result
        if (present(neval)) neval = evaluations
        if (present(subdivisions)) subdivisions = refine_count
        if (present(status)) status = integrate_status_success
    end subroutine cubature_core

    logical function valid_rule(rule, d) result(valid)
        character(len=*), intent(in) :: rule !! requested cubature rule name
        integer, intent(in) :: d !! integration dimension

        select case (trim(rule))
        case ('gk21', 'gk15', 'gauss-kronrod')
            valid = .true.
        case ('genz-malik')
            valid = d >= 2
        case default
            valid = .false.
        end select
    end function valid_rule

    pure logical function converged(estimate, error, atol, rtol) result(ok)
        real(dp), intent(in) :: estimate(:) !! current global integral estimate
        real(dp), intent(in) :: error(:) !! current global error estimate
        real(dp), intent(in) :: atol !! requested absolute tolerance
        real(dp), intent(in) :: rtol !! requested relative tolerance

        ok = all(error <= atol + rtol * abs(estimate))
    end function converged

    subroutine transformed_bounds(lower, upper, ta, tb, status)
        real(dp), intent(in) :: lower(:) !! normalized physical lower bounds
        real(dp), intent(in) :: upper(:) !! normalized physical upper bounds
        real(dp), intent(out) :: ta(:) !! transformed finite lower bounds
        real(dp), intent(out) :: tb(:) !! transformed finite upper bounds
        integer, intent(out) :: status !! integrate_status_* result

        integer :: i

        status = integrate_status_success
        do i = 1, size(lower)
            if (ieee_is_finite(lower(i)) .and. ieee_is_finite(upper(i))) then
                ta(i) = lower(i)
                tb(i) = upper(i)
            else if (.not. ieee_is_finite(lower(i)) .and. &
                    .not. ieee_is_finite(upper(i))) then
                if (lower(i) >= 0.0_dp .or. upper(i) <= 0.0_dp) then
                    status = integrate_status_invalid_argument
                    return
                end if
                ta(i) = -1.0_dp
                tb(i) = 1.0_dp
            else if (ieee_is_finite(lower(i)) .and. .not. ieee_is_finite(upper(i))) then
                if (upper(i) < 0.0_dp) then
                    status = integrate_status_invalid_argument
                    return
                end if
                ta(i) = 0.0_dp
                tb(i) = 1.0_dp
            else if (.not. ieee_is_finite(lower(i)) .and. ieee_is_finite(upper(i))) then
                if (lower(i) > 0.0_dp) then
                    status = integrate_status_invalid_argument
                    return
                end if
                ta(i) = 0.0_dp
                tb(i) = 1.0_dp
            else
                status = integrate_status_invalid_argument
                return
            end if
        end do
    end subroutine transformed_bounds

    subroutine map_point(lower, upper, point, transformed, usable)
        real(dp), intent(in) :: lower(:) !! normalized physical lower bounds
        real(dp), intent(in) :: upper(:) !! normalized physical upper bounds
        real(dp), intent(in) :: point(:) !! physical-space point
        real(dp), intent(out) :: transformed(:) !! transformed point coordinates
        logical, intent(out) :: usable !! true when the point has a finite interior image

        integer :: i

        usable = .true.
        do i = 1, size(lower)
            if (.not. ieee_is_finite(point(i))) then
                usable = .false.
                return
            end if
            if (point(i) < lower(i) .or. point(i) > upper(i)) then
                usable = .false.
                return
            end if
            if (ieee_is_finite(lower(i)) .and. ieee_is_finite(upper(i))) then
                transformed(i) = point(i)
            else if (.not. ieee_is_finite(lower(i)) .and. &
                    .not. ieee_is_finite(upper(i))) then
                if (point(i) == 0.0_dp) then
                    transformed(i) = 1.0_dp
                else if (point(i) > 0.0_dp) then
                    transformed(i) = 1.0_dp / (point(i) + 1.0_dp)
                else
                    transformed(i) = 1.0_dp / (point(i) - 1.0_dp)
                end if
            else if (ieee_is_finite(lower(i))) then
                transformed(i) = 1.0_dp / (point(i) - lower(i) + 1.0_dp)
            else
                transformed(i) = 1.0_dp / (upper(i) - point(i) + 1.0_dp)
            end if
        end do
    end subroutine map_point

    subroutine evaluate_transformed(f, lower, upper, t, values, evaluations, status)
        procedure(cubature_vector_integrand) :: f !! physical-space vector callback
        real(dp), intent(in) :: lower(:) !! normalized physical lower bounds
        real(dp), intent(in) :: upper(:) !! normalized physical upper bounds
        real(dp), intent(in) :: t(:) !! transformed finite-space evaluation point
        real(dp), intent(out) :: values(:) !! transformed integrand values
        integer, intent(inout) :: evaluations !! cumulative callback evaluations
        integer, intent(out) :: status !! integrate_status_* result

        integer :: i
        real(dp) :: jacobian
        real(dp), allocatable :: x(:)

        allocate(x(size(t)))
        jacobian = 1.0_dp
        do i = 1, size(t)
            if (ieee_is_finite(lower(i)) .and. ieee_is_finite(upper(i))) then
                x(i) = t(i)
            else if (.not. ieee_is_finite(lower(i)) .and. &
                    .not. ieee_is_finite(upper(i))) then
                if (t(i) == 0.0_dp) then
                    values = quiet_nan(0.0_dp)
                    status = integrate_status_numerical_failure
                    return
                end if
                x(i) = (1.0_dp - abs(t(i))) / t(i)
                jacobian = jacobian / (t(i) * t(i))
            else if (ieee_is_finite(lower(i))) then
                if (t(i) == 0.0_dp) then
                    values = quiet_nan(0.0_dp)
                    status = integrate_status_numerical_failure
                    return
                end if
                x(i) = lower(i) + (1.0_dp - t(i)) / t(i)
                jacobian = jacobian / (t(i) * t(i))
            else
                if (t(i) == 0.0_dp) then
                    values = quiet_nan(0.0_dp)
                    status = integrate_status_numerical_failure
                    return
                end if
                x(i) = upper(i) - (1.0_dp - t(i)) / t(i)
                jacobian = jacobian / (t(i) * t(i))
            end if
        end do

        call f(x, values)
        evaluations = evaluations + 1
        if (.not. all(ieee_is_finite(values)) .or. .not. ieee_is_finite(jacobian)) then
            status = integrate_status_numerical_failure
            return
        end if
        values = jacobian * values
        if (.not. all(ieee_is_finite(values))) then
            status = integrate_status_numerical_failure
            return
        end if
        status = integrate_status_success
    end subroutine evaluate_transformed

    subroutine set_empty_region(region, lower, upper, nout)
        type(cubature_region), intent(out) :: region !! region object to initialize
        real(dp), intent(in) :: lower(:) !! transformed region lower corner
        real(dp), intent(in) :: upper(:) !! transformed region upper corner
        integer, intent(in) :: nout !! number of callback outputs

        allocate(region%lower(size(lower)), region%upper(size(upper)))
        allocate(region%estimate(nout), region%error(nout))
        region%lower = lower
        region%upper = upper
        region%estimate = 0.0_dp
        region%error = 0.0_dp
        region%score = 0.0_dp
    end subroutine set_empty_region

    subroutine append_empty_region(regions, nregions, lower, upper, nout)
        type(cubature_region), allocatable, intent(inout) :: regions(:) !! active region storage
        integer, intent(inout) :: nregions !! active region count
        real(dp), intent(in) :: lower(:) !! transformed region lower corner
        real(dp), intent(in) :: upper(:) !! transformed region upper corner
        integer, intent(in) :: nout !! number of callback outputs

        type(cubature_region), allocatable :: grown(:)
        integer :: old_size
        integer :: target_size

        if (nregions == size(regions)) then
            old_size = size(regions)
            target_size = max(old_size + 1, 2 * old_size)
            allocate(grown(target_size))
            grown(1:old_size) = regions
            call move_alloc(grown, regions)
        end if
        nregions = nregions + 1
        call set_empty_region(regions(nregions), lower, upper, nout)
    end subroutine append_empty_region

    subroutine split_infinite_surfaces(regions, nregions, lower, upper, nout, status)
        type(cubature_region), allocatable, intent(inout) :: regions(:) !! active initial regions
        integer, intent(inout) :: nregions !! active region count
        real(dp), intent(in) :: lower(:) !! normalized physical lower bounds
        real(dp), intent(in) :: upper(:) !! normalized physical upper bounds
        integer, intent(in) :: nout !! number of callback outputs
        integer, intent(out) :: status !! integrate_status_* result

        integer :: axis
        integer :: i
        integer :: original_count
        real(dp), allocatable :: left_lower(:)
        real(dp), allocatable :: left_upper(:)
        real(dp), allocatable :: right_lower(:)
        real(dp), allocatable :: right_upper(:)

        status = integrate_status_success
        allocate(left_lower(size(lower)), left_upper(size(lower)))
        allocate(right_lower(size(lower)), right_upper(size(lower)))
        do axis = 1, size(lower)
            if (ieee_is_finite(lower(axis)) .or. ieee_is_finite(upper(axis))) cycle
            original_count = nregions
            do i = original_count, 1, -1
                if (regions(i)%lower(axis) < 0.0_dp .and. regions(i)%upper(axis) > 0.0_dp) then
                    left_lower = regions(i)%lower
                    left_upper = regions(i)%upper
                    right_lower = regions(i)%lower
                    right_upper = regions(i)%upper
                    left_upper(axis) = 0.0_dp
                    right_lower(axis) = 0.0_dp
                    call replace_region_bounds(regions(i), left_lower, left_upper, nout)
                    call append_empty_region(regions, nregions, right_lower, right_upper, nout)
                end if
            end do
        end do
    end subroutine split_infinite_surfaces

    subroutine replace_region_bounds(region, lower, upper, nout)
        type(cubature_region), intent(inout) :: region !! region to reset
        real(dp), intent(in) :: lower(:) !! new transformed lower corner
        real(dp), intent(in) :: upper(:) !! new transformed upper corner
        integer, intent(in) :: nout !! number of callback outputs

        if (allocated(region%lower)) deallocate(region%lower)
        if (allocated(region%upper)) deallocate(region%upper)
        if (allocated(region%estimate)) deallocate(region%estimate)
        if (allocated(region%error)) deallocate(region%error)
        call set_empty_region(region, lower, upper, nout)
    end subroutine replace_region_bounds

    subroutine split_user_points(regions, nregions, lower, upper, points, nout, status)
        type(cubature_region), allocatable, intent(inout) :: regions(:) !! active initial regions
        integer, intent(inout) :: nregions !! active region count
        real(dp), intent(in) :: lower(:) !! normalized physical lower bounds
        real(dp), intent(in) :: upper(:) !! normalized physical upper bounds
        real(dp), intent(in) :: points(:, :) !! physical-space points, one per column
        integer, intent(in) :: nout !! number of callback outputs
        integer, intent(out) :: status !! integrate_status_* result

        integer :: j
        logical :: usable
        real(dp), allocatable :: transformed(:)

        status = integrate_status_success
        allocate(transformed(size(lower)))
        do j = 1, size(points, 2)
            call map_point(lower, upper, points(:, j), transformed, usable)
            if (.not. usable) cycle
            call split_regions_at_point(regions, nregions, transformed, nout, status)
            if (status /= integrate_status_success) return
        end do
    end subroutine split_user_points

    subroutine split_regions_at_point(regions, nregions, point, nout, status)
        type(cubature_region), allocatable, intent(inout) :: regions(:) !! active initial regions
        integer, intent(inout) :: nregions !! active region count
        real(dp), intent(in) :: point(:) !! transformed point to place on region boundaries
        integer, intent(in) :: nout !! number of callback outputs
        integer, intent(out) :: status !! integrate_status_* result

        integer :: i
        integer :: original_count

        status = integrate_status_success
        original_count = nregions
        do i = original_count, 1, -1
            if (all(point >= regions(i)%lower) .and. all(point <= regions(i)%upper)) then
                call split_one_region_at_point(regions, nregions, i, point, nout, status)
                if (status /= integrate_status_success) return
            end if
        end do
    end subroutine split_regions_at_point

    subroutine split_one_region_at_point(regions, nregions, index, point, nout, status)
        type(cubature_region), allocatable, intent(inout) :: regions(:) !! active initial regions
        integer, intent(inout) :: nregions !! active region count
        integer, intent(in) :: index !! region to split
        real(dp), intent(in) :: point(:) !! transformed split point
        integer, intent(in) :: nout !! number of callback outputs
        integer, intent(out) :: status !! integrate_status_* result

        integer :: axis
        integer :: child
        integer :: i
        integer :: interior_count
        integer, allocatable :: interior_axes(:)
        real(dp), allocatable :: child_lower(:)
        real(dp), allocatable :: child_upper(:)
        real(dp), allocatable :: parent_lower(:)
        real(dp), allocatable :: parent_upper(:)

        status = integrate_status_success
        allocate(interior_axes(size(point)), parent_lower(size(point)), parent_upper(size(point)))
        parent_lower = regions(index)%lower
        parent_upper = regions(index)%upper
        interior_count = 0
        do axis = 1, size(point)
            if (point(axis) > parent_lower(axis) .and. point(axis) < parent_upper(axis)) then
                interior_count = interior_count + 1
                interior_axes(interior_count) = axis
            end if
        end do
        if (interior_count == 0) return
        if (interior_count >= bit_size(1) - 2) then
            status = integrate_status_invalid_argument
            return
        end if

        allocate(child_lower(size(point)), child_upper(size(point)))
        do child = 0, 2**interior_count - 1
            child_lower = parent_lower
            child_upper = parent_upper
            do i = 1, interior_count
                axis = interior_axes(i)
                if (btest(child, i - 1)) then
                    child_lower(axis) = point(axis)
                else
                    child_upper(axis) = point(axis)
                end if
            end do
            if (child == 0) then
                call replace_region_bounds(regions(index), child_lower, child_upper, nout)
            else
                call append_empty_region(regions, nregions, child_lower, child_upper, nout)
            end if
        end do
    end subroutine split_one_region_at_point

    integer function largest_error_region(regions, nregions) result(index)
        type(cubature_region), intent(in) :: regions(:) !! active regions
        integer, intent(in) :: nregions !! active region count

        integer :: i
        real(dp) :: best

        index = 1
        best = regions(1)%score
        do i = 2, nregions
            if (regions(i)%score > best) then
                best = regions(i)%score
                index = i
            end if
        end do
    end function largest_error_region

    subroutine refine_region(f, lower, upper, rule, regions, nregions, selected, nout, &
            evaluations, total_estimate, total_error, status)
        procedure(cubature_vector_integrand) :: f !! physical-space vector callback
        real(dp), intent(in) :: lower(:) !! normalized physical lower bounds
        real(dp), intent(in) :: upper(:) !! normalized physical upper bounds
        character(len=*), intent(in) :: rule !! selected cubature rule
        type(cubature_region), allocatable, intent(inout) :: regions(:) !! active region storage
        integer, intent(inout) :: nregions !! active region count
        integer, intent(in) :: selected !! active region to replace by children
        integer, intent(in) :: nout !! number of callback outputs
        integer, intent(inout) :: evaluations !! cumulative callback evaluations
        real(dp), intent(inout) :: total_estimate(:) !! global integral estimate
        real(dp), intent(inout) :: total_error(:) !! global error estimate
        integer, intent(out) :: status !! integrate_status_* result

        integer :: child
        integer :: d
        integer :: axis
        integer :: child_count
        real(dp), allocatable :: midpoint(:)
        real(dp), allocatable :: parent_lower(:)
        real(dp), allocatable :: parent_upper(:)
        real(dp), allocatable :: child_lower(:)
        real(dp), allocatable :: child_upper(:)
        type(cubature_region) :: new_region

        status = integrate_status_success
        d = size(lower)
        if (d >= bit_size(1) - 2) then
            status = integrate_status_invalid_argument
            return
        end if
        child_count = 2**d
        allocate(midpoint(d), parent_lower(d), parent_upper(d), child_lower(d), child_upper(d))
        parent_lower = regions(selected)%lower
        parent_upper = regions(selected)%upper
        midpoint = 0.5_dp * parent_lower + 0.5_dp * parent_upper

        do child = 0, child_count - 1
            child_lower = parent_lower
            child_upper = parent_upper
            do axis = 1, d
                if (btest(child, axis - 1)) then
                    child_lower(axis) = midpoint(axis)
                else
                    child_upper(axis) = midpoint(axis)
                end if
            end do
            call set_empty_region(new_region, child_lower, child_upper, nout)
            call evaluate_region(f, lower, upper, rule, new_region, evaluations, status)
            if (status /= integrate_status_success) return
            total_estimate = total_estimate + new_region%estimate
            total_error = total_error + new_region%error
            if (child == 0) then
                regions(selected) = new_region
            else
                call append_region(regions, nregions, new_region)
            end if
        end do
    end subroutine refine_region

    subroutine append_region(regions, nregions, region)
        type(cubature_region), allocatable, intent(inout) :: regions(:) !! active region storage
        integer, intent(inout) :: nregions !! active region count
        type(cubature_region), intent(in) :: region !! evaluated region to append

        type(cubature_region), allocatable :: grown(:)
        integer :: old_size
        integer :: target_size

        if (nregions == size(regions)) then
            old_size = size(regions)
            target_size = max(old_size + 1, 2 * old_size)
            allocate(grown(target_size))
            grown(1:old_size) = regions
            call move_alloc(grown, regions)
        end if
        nregions = nregions + 1
        regions(nregions) = region
    end subroutine append_region

    subroutine evaluate_region(f, lower, upper, rule, region, evaluations, status)
        procedure(cubature_vector_integrand) :: f !! physical-space vector callback
        real(dp), intent(in) :: lower(:) !! normalized physical lower bounds
        real(dp), intent(in) :: upper(:) !! normalized physical upper bounds
        character(len=*), intent(in) :: rule !! selected cubature rule
        type(cubature_region), intent(inout) :: region !! region to evaluate
        integer, intent(inout) :: evaluations !! cumulative callback evaluations
        integer, intent(out) :: status !! integrate_status_* result

        select case (trim(rule))
        case ('genz-malik')
            call evaluate_region_genz_malik(f, lower, upper, region, evaluations, status)
        case ('gk15')
            call evaluate_region_gk(f, lower, upper, region, 15, evaluations, status)
        case ('gk21', 'gauss-kronrod')
            call evaluate_region_gk(f, lower, upper, region, 21, evaluations, status)
        case default
            status = integrate_status_invalid_argument
        end select
        if (status == integrate_status_success) region%score = maxval(abs(region%error))
    end subroutine evaluate_region

    subroutine evaluate_region_genz_malik(f, lower, upper, region, evaluations, status)
        procedure(cubature_vector_integrand) :: f !! physical-space vector callback
        real(dp), intent(in) :: lower(:) !! normalized physical lower bounds
        real(dp), intent(in) :: upper(:) !! normalized physical upper bounds
        type(cubature_region), intent(inout) :: region !! transformed region and outputs
        integer, intent(inout) :: evaluations !! cumulative callback evaluations
        integer, intent(out) :: status !! integrate_status_* result

        integer :: d
        integer :: i
        integer :: j
        integer :: mask
        integer :: sign_i
        integer :: sign_j
        integer :: ncorner
        real(dp) :: half_volume
        real(dp) :: l2
        real(dp) :: l3
        real(dp) :: l4
        real(dp) :: l5
        real(dp) :: w1
        real(dp) :: w2
        real(dp) :: w3
        real(dp) :: w4
        real(dp) :: w5
        real(dp) :: v1
        real(dp) :: v2
        real(dp) :: v3
        real(dp) :: v4
        real(dp), allocatable :: center(:)
        real(dp), allocatable :: half_width(:)
        real(dp), allocatable :: point(:)
        real(dp), allocatable :: values(:)
        real(dp), allocatable :: high(:)
        real(dp), allocatable :: low(:)

        d = size(region%lower)
        if (d < 2 .or. d >= bit_size(1) - 2) then
            status = integrate_status_invalid_argument
            return
        end if
        ncorner = 2**d
        allocate(center(d), half_width(d), point(d), values(size(region%estimate)))
        allocate(high(size(region%estimate)), low(size(region%estimate)))
        center = 0.5_dp * region%lower + 0.5_dp * region%upper
        half_width = 0.5_dp * (region%upper - region%lower)
        half_volume = product(half_width)
        if (.not. ieee_is_finite(half_volume)) then
            status = integrate_status_numerical_failure
            return
        end if

        l2 = sqrt(9.0_dp / 70.0_dp)
        l3 = sqrt(9.0_dp / 10.0_dp)
        l4 = l3
        l5 = sqrt(9.0_dp / 19.0_dp)
        w1 = real(2**d, dp) * real(12824 - 9120*d + 400*d*d, dp) / 19683.0_dp
        w2 = real(2**d, dp) * 980.0_dp / 6561.0_dp
        w3 = real(2**d, dp) * real(1820 - 400*d, dp) / 19683.0_dp
        w4 = real(2**d, dp) * 200.0_dp / 19683.0_dp
        w5 = 6859.0_dp / 19683.0_dp
        v1 = real(2**d, dp) * real(729 - 950*d + 50*d*d, dp) / 729.0_dp
        v2 = real(2**d, dp) * 245.0_dp / 486.0_dp
        v3 = real(2**d, dp) * real(265 - 100*d, dp) / 1458.0_dp
        v4 = real(2**d, dp) * 25.0_dp / 729.0_dp
        high = 0.0_dp
        low = 0.0_dp

        point = center
        call evaluate_transformed(f, lower, upper, point, values, evaluations, status)
        if (status /= integrate_status_success) return
        high = high + w1 * values
        low = low + v1 * values

        do i = 1, d
            point = center
            point(i) = center(i) + l2 * half_width(i)
            call evaluate_transformed(f, lower, upper, point, values, evaluations, status)
            if (status /= integrate_status_success) return
            high = high + w2 * values
            low = low + v2 * values
            point(i) = center(i) - l2 * half_width(i)
            call evaluate_transformed(f, lower, upper, point, values, evaluations, status)
            if (status /= integrate_status_success) return
            high = high + w2 * values
            low = low + v2 * values
        end do

        do i = 1, d
            point = center
            point(i) = center(i) + l3 * half_width(i)
            call evaluate_transformed(f, lower, upper, point, values, evaluations, status)
            if (status /= integrate_status_success) return
            high = high + w3 * values
            low = low + v3 * values
            point(i) = center(i) - l3 * half_width(i)
            call evaluate_transformed(f, lower, upper, point, values, evaluations, status)
            if (status /= integrate_status_success) return
            high = high + w3 * values
            low = low + v3 * values
        end do

        do i = 1, d - 1
            do j = i + 1, d
                do sign_i = -1, 1, 2
                    do sign_j = -1, 1, 2
                        point = center
                        point(i) = center(i) + real(sign_i, dp) * l4 * half_width(i)
                        point(j) = center(j) + real(sign_j, dp) * l4 * half_width(j)
                        call evaluate_transformed(f, lower, upper, point, values, evaluations, status)
                        if (status /= integrate_status_success) return
                        high = high + w4 * values
                        low = low + v4 * values
                    end do
                end do
            end do
        end do

        do mask = 0, ncorner - 1
            point = center
            do i = 1, d
                if (btest(mask, i - 1)) then
                    point(i) = center(i) + l5 * half_width(i)
                else
                    point(i) = center(i) - l5 * half_width(i)
                end if
            end do
            call evaluate_transformed(f, lower, upper, point, values, evaluations, status)
            if (status /= integrate_status_success) return
            high = high + w5 * values
        end do

        region%estimate = half_volume * high
        region%error = abs(half_volume * (high - low))
        if (.not. all(ieee_is_finite(region%estimate)) .or. &
                .not. all(ieee_is_finite(region%error))) then
            status = integrate_status_numerical_failure
            return
        end if
        status = integrate_status_success
    end subroutine evaluate_region_genz_malik

    subroutine evaluate_region_gk(f, lower, upper, region, order, evaluations, status)
        procedure(cubature_vector_integrand) :: f !! physical-space vector callback
        real(dp), intent(in) :: lower(:) !! normalized physical lower bounds
        real(dp), intent(in) :: upper(:) !! normalized physical upper bounds
        type(cubature_region), intent(inout) :: region !! transformed region and outputs
        integer, intent(in) :: order !! Kronrod order, 15 or 21
        integer, intent(inout) :: evaluations !! cumulative callback evaluations
        integer, intent(out) :: status !! integrate_status_* result

        integer :: root_status
        integer :: lower_order
        real(dp), allocatable :: gauss_nodes(:)
        real(dp), allocatable :: gauss_weights(:)
        real(dp), pointer :: kronrod_nodes(:)
        real(dp), pointer :: kronrod_weights(:)
        real(dp), target :: local_nodes21(21)
        real(dp), target :: local_weights21(21)
        real(dp), target :: local_nodes15(15)
        real(dp), target :: local_weights15(15)
        real(dp), allocatable :: high(:)
        real(dp), allocatable :: low(:)

        local_nodes21 = gk21_nodes
        local_weights21 = gk21_weights
        local_nodes15 = gk15_nodes
        local_weights15 = gk15_weights
        if (order == 21) then
            kronrod_nodes => local_nodes21
            kronrod_weights => local_weights21
            lower_order = 10
        else if (order == 15) then
            kronrod_nodes => local_nodes15
            kronrod_weights => local_weights15
            lower_order = 7
        else
            status = integrate_status_invalid_argument
            return
        end if

        call roots_legendre(lower_order, gauss_nodes, gauss_weights, status=root_status)
        if (root_status /= orthogonal_roots_status_success) then
            status = integrate_status_numerical_failure
            return
        end if
        allocate(high(size(region%estimate)), low(size(region%estimate)))
        call product_rule_sum(f, lower, upper, region%lower, region%upper, &
            kronrod_nodes, kronrod_weights, high, evaluations, status)
        if (status /= integrate_status_success) return
        call product_rule_sum(f, lower, upper, region%lower, region%upper, &
            gauss_nodes, gauss_weights, low, evaluations, status)
        if (status /= integrate_status_success) return
        region%estimate = high
        region%error = abs(high - low)
        status = integrate_status_success
    end subroutine evaluate_region_gk

    subroutine product_rule_sum(f, lower, upper, region_lower, region_upper, nodes, weights, &
            estimate, evaluations, status)
        procedure(cubature_vector_integrand) :: f !! physical-space vector callback
        real(dp), intent(in) :: lower(:) !! normalized physical lower bounds
        real(dp), intent(in) :: upper(:) !! normalized physical upper bounds
        real(dp), intent(in) :: region_lower(:) !! transformed region lower corner
        real(dp), intent(in) :: region_upper(:) !! transformed region upper corner
        real(dp), intent(in) :: nodes(:) !! one-dimensional fixed-rule nodes on [-1,1]
        real(dp), intent(in) :: weights(:) !! one-dimensional fixed-rule weights
        real(dp), intent(out) :: estimate(:) !! product-rule integral estimate
        integer, intent(inout) :: evaluations !! cumulative callback evaluations
        integer, intent(out) :: status !! integrate_status_* result

        integer :: d
        integer :: axis
        integer :: digit
        integer :: n
        integer :: total
        integer :: linear
        real(dp) :: scale_factor
        real(dp) :: weight_product
        real(dp), allocatable :: center(:)
        real(dp), allocatable :: half_width(:)
        real(dp), allocatable :: point(:)
        real(dp), allocatable :: values(:)

        d = size(region_lower)
        n = size(nodes)
        if (size(weights) /= n .or. n < 1) then
            status = integrate_status_invalid_argument
            return
        end if
        if (.not. safe_integer_power(n, d, total)) then
            status = integrate_status_invalid_argument
            return
        end if
        allocate(center(d), half_width(d), point(d), values(size(estimate)))
        center = 0.5_dp * region_lower + 0.5_dp * region_upper
        half_width = 0.5_dp * (region_upper - region_lower)
        scale_factor = product(half_width)
        estimate = 0.0_dp
        do linear = 0, total - 1
            digit = linear
            weight_product = 1.0_dp
            do axis = 1, d
                point(axis) = center(axis) + half_width(axis) * nodes(mod(digit, n) + 1)
                weight_product = weight_product * weights(mod(digit, n) + 1)
                digit = digit / n
            end do
            call evaluate_transformed(f, lower, upper, point, values, evaluations, status)
            if (status /= integrate_status_success) return
            estimate = estimate + weight_product * values
        end do
        estimate = scale_factor * estimate
        if (.not. all(ieee_is_finite(estimate))) then
            status = integrate_status_numerical_failure
            return
        end if
        status = integrate_status_success
    end subroutine product_rule_sum

    logical function safe_integer_power(base, exponent, value) result(ok)
        integer, intent(in) :: base !! positive integer base
        integer, intent(in) :: exponent !! nonnegative exponent
        integer, intent(out) :: value !! exact power when representable

        integer :: i

        value = 1
        ok = base > 0 .and. exponent >= 0
        if (.not. ok) return
        do i = 1, exponent
            if (value > huge(value) / base) then
                ok = .false.
                return
            end if
            value = value * base
        end do
    end function safe_integer_power

end module scifort_integrate_cubature
