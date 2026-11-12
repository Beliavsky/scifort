! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Correctness-first vector-valued adaptive quadrature. Each component is
! integrated through SciFort's scalar QUADPACK facade. This preserves the
! scipy.integrate.quad_vec mathematical contract while deliberately trading
! shared function evaluations for a small, idiomatic Fortran implementation.

module scifort_integrate_quad_vec
    use scifort_integrate_quadpack, only : quad, quad_status_invalid_argument, &
        quad_status_success
    use scifort_kinds, only : dp
    use scifort_math, only : quiet_nan
    implicit none
    private

    public :: quad_vec

    abstract interface
        subroutine vector_integrand(x, values)
            import :: dp
            implicit none
            real(dp), intent(in) :: x !! scalar integration coordinate
            real(dp), intent(out) :: values(:) !! vector-valued integrand at x
        end subroutine vector_integrand
    end interface

contains

    subroutine quad_vec(f, a, b, result, abserr, epsabs, epsrel, limit, points, norm, &
            neval, status)
        procedure(vector_integrand) :: f !! vector-valued integrand callback
        real(dp), intent(in) :: a !! lower integration limit; may be infinite
        real(dp), intent(in) :: b !! upper integration limit; may be infinite
        real(dp), intent(out) :: result(:) !! integral of each callback component
        real(dp), intent(out) :: abserr !! norm of componentwise absolute-error estimates
        real(dp), intent(in), optional :: epsabs !! vector absolute tolerance; default 1e-200
        real(dp), intent(in), optional :: epsrel !! vector relative tolerance; default 1e-8
        integer, intent(in), optional :: limit !! maximum scalar subintervals per component; default 10000
        real(dp), intent(in), optional :: points(:) !! optional break points passed to each scalar integral
        character(len=*), intent(in), optional :: norm !! error norm, '2' or 'max'; default '2'
        integer, intent(out), optional :: neval !! total callback evaluations across all components
        integer, intent(out), optional :: status !! first nonzero scalar status, or zero on success

        character(len=8) :: selected_norm
        integer :: component
        integer :: component_evaluations
        integer :: component_status
        integer :: max_intervals
        integer :: total_evaluations
        integer :: overall_status
        real(dp) :: absolute_tolerance
        real(dp) :: component_absolute_tolerance
        real(dp) :: relative_tolerance
        real(dp), allocatable :: component_errors(:)
        real(dp), allocatable :: work(:)

        result = quiet_nan(0.0_dp)
        abserr = quiet_nan(0.0_dp)
        total_evaluations = 0
        overall_status = quad_status_invalid_argument
        if (present(neval)) neval = total_evaluations
        if (present(status)) status = overall_status
        if (size(result) < 1) return

        absolute_tolerance = 1.0e-200_dp
        relative_tolerance = 1.0e-8_dp
        max_intervals = 10000
        selected_norm = '2'
        if (present(epsabs)) absolute_tolerance = epsabs
        if (present(epsrel)) relative_tolerance = epsrel
        if (present(limit)) max_intervals = limit
        if (present(norm)) selected_norm = trim(adjustl(norm))
        if (max_intervals < 1) return
        if (absolute_tolerance < 0.0_dp .or. relative_tolerance < 0.0_dp) return
        if (trim(selected_norm) /= '2' .and. trim(selected_norm) /= 'max') return

        allocate(component_errors(size(result)), work(size(result)))
        component_errors = 0.0_dp
        overall_status = quad_status_success
        if (trim(selected_norm) == '2') then
            component_absolute_tolerance = absolute_tolerance / sqrt(real(size(result), dp))
        else
            component_absolute_tolerance = absolute_tolerance
        end if

        do component = 1, size(result)
            if (present(points)) then
                call quad(component_integrand, a, b, result(component), component_errors(component), &
                    epsabs=component_absolute_tolerance, epsrel=relative_tolerance, &
                    limit=max_intervals, neval=component_evaluations, status=component_status, &
                    points=points)
            else
                call quad(component_integrand, a, b, result(component), component_errors(component), &
                    epsabs=component_absolute_tolerance, epsrel=relative_tolerance, &
                    limit=max_intervals, neval=component_evaluations, status=component_status)
            end if
            total_evaluations = total_evaluations + component_evaluations
            if (overall_status == quad_status_success .and. &
                    component_status /= quad_status_success) overall_status = component_status
        end do

        if (trim(selected_norm) == 'max') then
            abserr = maxval(component_errors)
        else
            abserr = sqrt(sum(component_errors * component_errors))
        end if
        if (present(neval)) neval = total_evaluations
        if (present(status)) status = overall_status

    contains

        function component_integrand(x) result(value)
            real(dp), intent(in) :: x !! scalar coordinate supplied by QUADPACK
            real(dp) :: value

            call f(x, work)
            value = work(component)
        end function component_integrand

    end subroutine quad_vec

end module scifort_integrate_quad_vec
