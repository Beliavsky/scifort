! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_integrate_nested
    use scifort_integrate_quadpack, only : quad, quad_status_invalid_argument, &
        quad_status_success
    use scifort_kinds, only : dp
    use scifort_math, only : quiet_nan
    implicit none
    private

    public :: dblquad
    public :: nquad
    public :: tplquad

    abstract interface
        function integrand_2d(y, x) result(value)
            import :: dp
            implicit none
            real(dp), intent(in) :: y !! inner integration coordinate
            real(dp), intent(in) :: x !! outer integration coordinate
            real(dp) :: value
        end function integrand_2d

        function integrand_3d(z, y, x) result(value)
            import :: dp
            implicit none
            real(dp), intent(in) :: z !! innermost integration coordinate
            real(dp), intent(in) :: y !! middle integration coordinate
            real(dp), intent(in) :: x !! outer integration coordinate
            real(dp) :: value
        end function integrand_3d

        function bound_1d(x) result(value)
            import :: dp
            implicit none
            real(dp), intent(in) :: x !! outer integration coordinate
            real(dp) :: value
        end function bound_1d

        function bound_2d(x, y) result(value)
            import :: dp
            implicit none
            real(dp), intent(in) :: x !! outer integration coordinate
            real(dp), intent(in) :: y !! middle integration coordinate
            real(dp) :: value
        end function bound_2d

        function integrand_nd(x) result(value)
            import :: dp
            implicit none
            real(dp), intent(in) :: x(:) !! coordinates, with x(1) integrated innermost
            real(dp) :: value
        end function integrand_nd

        subroutine bounds_nd(axis, x, lower, upper)
            import :: dp
            implicit none
            integer, intent(in) :: axis !! integration axis, 1 innermost and size(x) outermost
            real(dp), intent(in) :: x(:) !! coordinate workspace; x(axis+1:) contains outer coordinates
            real(dp), intent(out) :: lower !! lower bound for this axis; may be infinite
            real(dp), intent(out) :: upper !! upper bound for this axis; may be infinite
        end subroutine bounds_nd
    end interface

contains

    subroutine dblquad(f, a, b, gfun, hfun, result, abserr, epsabs, epsrel, limit, status)
        procedure(integrand_2d) :: f !! integrand evaluated as f(y, x)
        real(dp), intent(in) :: a !! lower outer x limit; may be infinite
        real(dp), intent(in) :: b !! upper outer x limit; may be infinite
        procedure(bound_1d) :: gfun !! lower inner y bound as a function of x
        procedure(bound_1d) :: hfun !! upper inner y bound as a function of x
        real(dp), intent(out) :: result !! nested adaptive quadrature estimate
        real(dp), intent(out) :: abserr !! maximum QUADPACK error estimate over all levels
        real(dp), intent(in), optional :: epsabs !! tolerance passed to each integration level
        real(dp), intent(in), optional :: epsrel !! relative tolerance passed to each integration level
        integer, intent(in), optional :: limit !! maximum subintervals at each level; default 50
        integer, intent(out), optional :: status !! first nonzero nested QUADPACK status, else outer status

        integer :: inner_status
        integer :: max_intervals
        integer :: nested_status
        integer :: outer_status
        real(dp) :: absolute_tolerance
        real(dp) :: current_x
        real(dp) :: maximum_error
        real(dp) :: outer_error
        real(dp) :: relative_tolerance

        absolute_tolerance = 1.49e-8_dp
        relative_tolerance = 1.49e-8_dp
        max_intervals = 50
        if (present(epsabs)) absolute_tolerance = epsabs
        if (present(epsrel)) relative_tolerance = epsrel
        if (present(limit)) max_intervals = limit
        nested_status = quad_status_success
        maximum_error = 0.0_dp
        current_x = 0.0_dp

        call quad(outer_integrand, a, b, result, outer_error, epsabs=absolute_tolerance, &
            epsrel=relative_tolerance, limit=max_intervals, status=outer_status)
        abserr = max(outer_error, maximum_error)
        if (outer_status /= quad_status_success) then
            if (present(status)) status = outer_status
        else
            if (present(status)) status = nested_status
        end if

    contains

        function outer_integrand(x) result(value)
            real(dp), intent(in) :: x !! current outer integration coordinate
            real(dp) :: value
            real(dp) :: inner_error

            current_x = x
            call quad(inner_integrand, gfun(x), hfun(x), value, inner_error, &
                epsabs=absolute_tolerance, epsrel=relative_tolerance, &
                limit=max_intervals, status=inner_status)
            maximum_error = max(maximum_error, inner_error)
            if (nested_status == quad_status_success .and. &
                    inner_status /= quad_status_success) nested_status = inner_status
        end function outer_integrand

        function inner_integrand(y) result(value)
            real(dp), intent(in) :: y !! current inner integration coordinate
            real(dp) :: value

            value = f(y, current_x)
        end function inner_integrand

    end subroutine dblquad

    subroutine nquad(f, ndim, ranges, result, abserr, epsabs, epsrel, limit, status)
        procedure(integrand_nd) :: f !! integrand evaluated at x(1:ndim), x(1) innermost
        integer, intent(in) :: ndim !! number of integration dimensions, at least one
        procedure(bounds_nd) :: ranges !! callback returning bounds for each axis
        real(dp), intent(out) :: result !! nested adaptive quadrature estimate
        real(dp), intent(out) :: abserr !! maximum QUADPACK error estimate over all levels
        real(dp), intent(in), optional :: epsabs !! tolerance passed to each integration level
        real(dp), intent(in), optional :: epsrel !! relative tolerance passed to each integration level
        integer, intent(in), optional :: limit !! maximum subintervals at each level; default 50
        integer, intent(out), optional :: status !! first nonzero QUADPACK status across all levels

        integer :: max_intervals
        integer :: nested_status
        real(dp) :: absolute_tolerance
        real(dp) :: maximum_error
        real(dp) :: relative_tolerance
        real(dp), allocatable :: x(:)

        result = quiet_nan(0.0_dp)
        abserr = quiet_nan(0.0_dp)
        if (present(status)) status = quad_status_invalid_argument
        if (ndim < 1) return

        absolute_tolerance = 1.49e-8_dp
        relative_tolerance = 1.49e-8_dp
        max_intervals = 50
        if (present(epsabs)) absolute_tolerance = epsabs
        if (present(epsrel)) relative_tolerance = epsrel
        if (present(limit)) max_intervals = limit
        if (max_intervals < 1) return

        allocate(x(ndim))
        x = 0.0_dp
        maximum_error = 0.0_dp
        nested_status = quad_status_success
        call nquad_level(f, ranges, x, ndim, absolute_tolerance, relative_tolerance, &
            max_intervals, result, maximum_error, nested_status)
        abserr = maximum_error
        if (present(status)) status = nested_status
    end subroutine nquad

    recursive subroutine nquad_level(f, ranges, x, axis, epsabs, epsrel, limit, &
            result, maximum_error, nested_status)
        procedure(integrand_nd) :: f !! n-dimensional integrand callback
        procedure(bounds_nd) :: ranges !! axis-bound callback
        real(dp), intent(inout) :: x(:) !! coordinate workspace shared across integration levels
        integer, intent(in) :: axis !! current integration axis, decreasing toward one
        real(dp), intent(in) :: epsabs !! absolute tolerance for this and inner levels
        real(dp), intent(in) :: epsrel !! relative tolerance for this and inner levels
        integer, intent(in) :: limit !! maximum subintervals at each level
        real(dp), intent(out) :: result !! integral over this and all inner axes
        real(dp), intent(inout) :: maximum_error !! largest error estimate seen so far
        integer, intent(inout) :: nested_status !! first nonzero status seen so far

        integer :: level_status
        real(dp) :: level_error
        real(dp) :: lower
        real(dp) :: upper

        call ranges(axis, x, lower, upper)
        call quad(level_integrand, lower, upper, result, level_error, epsabs=epsabs, &
            epsrel=epsrel, limit=limit, status=level_status)
        maximum_error = max(maximum_error, level_error)
        if (nested_status == quad_status_success .and. &
                level_status /= quad_status_success) nested_status = level_status

    contains

        recursive function level_integrand(coordinate) result(value)
            real(dp), intent(in) :: coordinate !! coordinate on the current integration axis
            real(dp) :: value

            x(axis) = coordinate
            if (axis == 1) then
                value = f(x)
            else
                call nquad_level(f, ranges, x, axis - 1, epsabs, epsrel, limit, &
                    value, maximum_error, nested_status)
            end if
        end function level_integrand

    end subroutine nquad_level

    subroutine tplquad(f, a, b, gfun, hfun, qfun, rfun, result, abserr, &
            epsabs, epsrel, limit, status)
        procedure(integrand_3d) :: f !! integrand evaluated as f(z, y, x)
        real(dp), intent(in) :: a !! lower outer x limit; may be infinite
        real(dp), intent(in) :: b !! upper outer x limit; may be infinite
        procedure(bound_1d) :: gfun !! lower middle y bound as a function of x
        procedure(bound_1d) :: hfun !! upper middle y bound as a function of x
        procedure(bound_2d) :: qfun !! lower inner z bound as a function of x and y
        procedure(bound_2d) :: rfun !! upper inner z bound as a function of x and y
        real(dp), intent(out) :: result !! nested adaptive quadrature estimate
        real(dp), intent(out) :: abserr !! maximum QUADPACK error estimate over all levels
        real(dp), intent(in), optional :: epsabs !! tolerance passed to each integration level
        real(dp), intent(in), optional :: epsrel !! relative tolerance passed to each integration level
        integer, intent(in), optional :: limit !! maximum subintervals at each level; default 50
        integer, intent(out), optional :: status !! first nonzero nested QUADPACK status, else outer status

        integer :: inner_status
        integer :: max_intervals
        integer :: middle_status
        integer :: nested_status
        integer :: outer_status
        real(dp) :: absolute_tolerance
        real(dp) :: current_x
        real(dp) :: current_y
        real(dp) :: maximum_error
        real(dp) :: outer_error
        real(dp) :: relative_tolerance

        absolute_tolerance = 1.49e-8_dp
        relative_tolerance = 1.49e-8_dp
        max_intervals = 50
        if (present(epsabs)) absolute_tolerance = epsabs
        if (present(epsrel)) relative_tolerance = epsrel
        if (present(limit)) max_intervals = limit
        nested_status = quad_status_success
        maximum_error = 0.0_dp
        current_x = 0.0_dp
        current_y = 0.0_dp

        call quad(outer_integrand, a, b, result, outer_error, epsabs=absolute_tolerance, &
            epsrel=relative_tolerance, limit=max_intervals, status=outer_status)
        abserr = max(outer_error, maximum_error)
        if (outer_status /= quad_status_success) then
            if (present(status)) status = outer_status
        else
            if (present(status)) status = nested_status
        end if

    contains

        function outer_integrand(x) result(value)
            real(dp), intent(in) :: x !! current outer integration coordinate
            real(dp) :: value
            real(dp) :: middle_error

            current_x = x
            call quad(middle_integrand, gfun(x), hfun(x), value, middle_error, &
                epsabs=absolute_tolerance, epsrel=relative_tolerance, &
                limit=max_intervals, status=middle_status)
            maximum_error = max(maximum_error, middle_error)
            if (nested_status == quad_status_success .and. &
                    middle_status /= quad_status_success) nested_status = middle_status
        end function outer_integrand

        function middle_integrand(y) result(value)
            real(dp), intent(in) :: y !! current middle integration coordinate
            real(dp) :: value
            real(dp) :: inner_error

            current_y = y
            call quad(inner_integrand, qfun(current_x, y), rfun(current_x, y), &
                value, inner_error, epsabs=absolute_tolerance, epsrel=relative_tolerance, &
                limit=max_intervals, status=inner_status)
            maximum_error = max(maximum_error, inner_error)
            if (nested_status == quad_status_success .and. &
                    inner_status /= quad_status_success) nested_status = inner_status
        end function middle_integrand

        function inner_integrand(z) result(value)
            real(dp), intent(in) :: z !! current innermost integration coordinate
            real(dp) :: value

            value = f(z, current_y, current_x)
        end function inner_integrand

    end subroutine tplquad

end module scifort_integrate_nested
