! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! SciPy-style adaptive one-dimensional integration facade over the external
! BSD-3-Clause modern QUADPACK package by Jacob Williams and contributors.
! The QUADPACK implementation is not copied into SciFort; fpm resolves it as
! the pinned dependency documented in fpm.toml and CODE_PROVENANCE.md.

module scifort_integrate_quadpack
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use quadpack, only : dqagi, dqagp, dqags, dqawc, dqawf, dqawo, dqaws
    use scifort_kinds, only : dp
    use scifort_math, only : quiet_nan
    implicit none
    private

    integer, parameter, public :: quad_status_success = 0
    integer, parameter, public :: quad_status_max_subdivisions = 1
    integer, parameter, public :: quad_status_roundoff = 2
    integer, parameter, public :: quad_status_bad_integrand = 3
    integer, parameter, public :: quad_status_nonconvergence = 4
    integer, parameter, public :: quad_status_divergent = 5
    integer, parameter, public :: quad_status_invalid_argument = 6
    integer, parameter, public :: quad_status_cycle_failure = 7

    public :: quad

    abstract interface
        function quad_integrand(x) result(value)
            import :: dp
            implicit none
            real(dp), intent(in) :: x !! integration point
            real(dp) :: value
        end function quad_integrand
    end interface

contains

    recursive subroutine quad(f, a, b, result, abserr, epsabs, epsrel, limit, neval, status, &
            points, weight, wvar, wvar2, limlst, maxp1)
        procedure(quad_integrand) :: f !! scalar integrand function
        real(dp), intent(in) :: a !! lower integration limit; may be infinite
        real(dp), intent(in) :: b !! upper integration limit; may be infinite
        real(dp), intent(out) :: result !! adaptive quadrature estimate including any requested weight
        real(dp), intent(out) :: abserr !! QUADPACK absolute-error estimate
        real(dp), intent(in), optional :: epsabs !! absolute tolerance; default 1.49e-8
        real(dp), intent(in), optional :: epsrel !! relative tolerance; default 1.49e-8
        integer, intent(in), optional :: limit !! maximum adaptive subintervals; default 50
        integer, intent(out), optional :: neval !! number of integrand evaluations
        integer, intent(out), optional :: status !! QUADPACK-style status; zero on success
        real(dp), intent(in), optional :: points(:) !! finite-interval break points for local difficulties
        character(len=*), intent(in), optional :: weight !! optional SciPy weight name
        real(dp), intent(in), optional :: wvar !! omega, alpha, or Cauchy pole associated with weight
        real(dp), intent(in), optional :: wvar2 !! beta for algebraic endpoint weights
        integer, intent(in), optional :: limlst !! maximum Fourier cycles for infinite sin/cos weights; default 50
        integer, intent(in), optional :: maxp1 !! Chebyshev-moment storage limit for sin/cos weights; default 50

        character(len=16) :: weight_name
        integer :: inf
        integer :: ier
        integer :: last
        integer :: leniw
        integer :: lenw
        integer :: max_cycles
        integer :: max_intervals
        integer :: max_moments
        integer :: npts2
        integer :: evaluations
        integer :: integr
        integer :: lst
        integer, allocatable :: iwork(:)
        real(dp) :: absolute_tolerance
        real(dp) :: relative_tolerance
        real(dp) :: bound
        real(dp) :: lower
        real(dp) :: orientation
        real(dp) :: upper
        real(dp), allocatable :: points_work(:)
        real(dp), allocatable :: work(:)

        result = quiet_nan(0.0_dp)
        abserr = quiet_nan(0.0_dp)
        evaluations = 0
        ier = quad_status_invalid_argument
        if (present(neval)) neval = evaluations
        if (present(status)) status = ier

        absolute_tolerance = 1.49e-8_dp
        relative_tolerance = 1.49e-8_dp
        max_intervals = 50
        max_cycles = 50
        max_moments = 50
        if (present(epsabs)) absolute_tolerance = epsabs
        if (present(epsrel)) relative_tolerance = epsrel
        if (present(limit)) max_intervals = limit
        if (present(limlst)) max_cycles = limlst
        if (present(maxp1)) max_moments = maxp1

        if (ieee_is_nan(a) .or. ieee_is_nan(b) .or. &
                ieee_is_nan(absolute_tolerance) .or. ieee_is_nan(relative_tolerance) .or. &
                max_intervals < 1) then
            return
        end if
        if (present(points)) then
            if (any(ieee_is_nan(points))) return
        end if
        if (a == b) then
            result = 0.0_dp
            abserr = 0.0_dp
            ier = quad_status_success
            if (present(status)) status = ier
            return
        end if

        orientation = 1.0_dp
        lower = a
        upper = b
        if (b < a) then
            lower = b
            upper = a
            orientation = -1.0_dp
        end if

        if (present(weight)) then
            if (present(points)) return
            weight_name = trim(adjustl(weight))
            select case (trim(weight_name))
            case ('cos', 'sin')
                if (.not. present(wvar)) return
                if (ieee_is_nan(wvar)) return
                if (max_moments < 1) return
                if (trim(weight_name) == 'cos') then
                    integr = 1
                else
                    integr = 2
                end if
                if (ieee_is_finite(lower) .and. ieee_is_finite(upper)) then
                    leniw = 2 * max_intervals
                    lenw = 2 * leniw + 25 * max_moments
                    allocate(iwork(leniw), work(lenw))
                    call dqawo(f, lower, upper, wvar, integr, absolute_tolerance, &
                        relative_tolerance, result, abserr, evaluations, ier, leniw, &
                        max_moments, lenw, last, iwork, work)
                elseif (.not. ieee_is_finite(lower) .and. .not. ieee_is_finite(upper)) then
                    return
                else
                    if (absolute_tolerance <= 0.0_dp .or. max_cycles < 3) return
                    leniw = 2 * max_intervals + max_cycles
                    lenw = 2 * leniw + 25 * max_moments
                    allocate(iwork(leniw), work(lenw))
                    if (ieee_is_finite(lower)) then
                        call dqawf(f, lower, wvar, integr, absolute_tolerance, result, &
                            abserr, evaluations, ier, max_cycles, lst, leniw, max_moments, &
                            lenw, iwork, work)
                    elseif (integr == 1) then
                        call dqawf(reflected_cos_integrand, -upper, wvar, integr, &
                            absolute_tolerance, result, abserr, evaluations, ier, &
                            max_cycles, lst, leniw, max_moments, lenw, iwork, work)
                    else
                        call dqawf(reflected_sin_integrand, -upper, wvar, integr, &
                            absolute_tolerance, result, abserr, evaluations, ier, &
                            max_cycles, lst, leniw, max_moments, lenw, iwork, work)
                    end if
                end if
            case ('alg', 'alg-loga', 'alg-logb', 'alg-log')
                if (.not. present(wvar) .or. .not. present(wvar2)) return
                if (ieee_is_nan(wvar) .or. ieee_is_nan(wvar2)) return
                if (.not. ieee_is_finite(lower) .or. .not. ieee_is_finite(upper)) return
                if (wvar <= -1.0_dp .or. wvar2 <= -1.0_dp .or. max_intervals < 2) return
                select case (trim(weight_name))
                case ('alg')
                    integr = 1
                case ('alg-loga')
                    integr = 2
                case ('alg-logb')
                    integr = 3
                case default
                    integr = 4
                end select
                lenw = 4 * max_intervals
                allocate(iwork(max_intervals), work(lenw))
                call dqaws(f, lower, upper, wvar, wvar2, integr, absolute_tolerance, &
                    relative_tolerance, result, abserr, evaluations, ier, max_intervals, &
                    lenw, last, iwork, work)
            case ('cauchy')
                if (.not. present(wvar)) return
                if (ieee_is_nan(wvar)) return
                if (.not. ieee_is_finite(lower) .or. .not. ieee_is_finite(upper)) return
                if (wvar == lower .or. wvar == upper) return
                lenw = 4 * max_intervals
                allocate(iwork(max_intervals), work(lenw))
                call dqawc(f, lower, upper, wvar, absolute_tolerance, relative_tolerance, &
                    result, abserr, evaluations, ier, max_intervals, lenw, last, &
                    iwork, work)
            case default
                return
            end select
        elseif (present(points)) then
            if (.not. ieee_is_finite(lower) .or. .not. ieee_is_finite(upper)) return
            if (size(points) >= max_intervals) return
            if (size(points) > 0) then
                if (minval(points) < lower .or. maxval(points) > upper) return
            end if
            call prepare_break_points(points, lower, upper, points_work)
            npts2 = size(points_work)
            if (max_intervals < npts2 - 1) return
            leniw = 2 * max_intervals + npts2
            lenw = 2 * leniw - npts2
            allocate(iwork(leniw), work(lenw))
            call dqagp(f, lower, upper, npts2, points_work, absolute_tolerance, &
                relative_tolerance, result, abserr, evaluations, ier, leniw, lenw, &
                last, iwork, work)
        elseif (ieee_is_finite(lower) .and. ieee_is_finite(upper)) then
            lenw = 4 * max_intervals
            allocate(iwork(max_intervals), work(lenw))
            call dqags(f, lower, upper, absolute_tolerance, relative_tolerance, result, &
                abserr, evaluations, ier, max_intervals, lenw, last, iwork, work)
        else
            if (.not. ieee_is_finite(lower) .and. .not. ieee_is_finite(upper)) then
                inf = 2
                bound = 0.0_dp
            elseif (ieee_is_finite(lower)) then
                inf = 1
                bound = lower
            else
                inf = -1
                bound = upper
            end if
            lenw = 4 * max_intervals
            allocate(iwork(max_intervals), work(lenw))
            call dqagi(f, bound, inf, absolute_tolerance, relative_tolerance, result, &
                abserr, evaluations, ier, max_intervals, lenw, last, iwork, work)
        end if

        result = orientation * result
        if (present(neval)) neval = evaluations
        if (present(status)) status = ier

    contains

        function reflected_cos_integrand(x) result(value)
            real(dp), intent(in) :: x !! reflected nonnegative integration coordinate
            real(dp) :: value

            value = f(-x)
        end function reflected_cos_integrand

        function reflected_sin_integrand(x) result(value)
            real(dp), intent(in) :: x !! reflected nonnegative integration coordinate
            real(dp) :: value

            value = -f(-x)
        end function reflected_sin_integrand

    end subroutine quad

    subroutine prepare_break_points(points, lower, upper, prepared)
        real(dp), intent(in) :: points(:) !! user-provided break points including optional endpoints
        real(dp), intent(in) :: lower !! normalized finite lower integration limit
        real(dp), intent(in) :: upper !! normalized finite upper integration limit
        real(dp), allocatable, intent(out) :: prepared(:) !! sorted unique interior points plus two work slots

        integer :: i
        integer :: j
        integer :: n
        integer :: unique_count
        real(dp) :: temp
        real(dp), allocatable :: interior(:)

        n = count(points > lower .and. points < upper)
        allocate(interior(n))
        j = 0
        do i = 1, size(points)
            if (points(i) > lower .and. points(i) < upper) then
                j = j + 1
                interior(j) = points(i)
            end if
        end do

        do i = 2, n
            temp = interior(i)
            j = i - 1
            do while (j >= 1)
                if (interior(j) <= temp) exit
                interior(j + 1) = interior(j)
                j = j - 1
            end do
            interior(j + 1) = temp
        end do

        unique_count = 0
        do i = 1, n
            if (unique_count == 0) then
                unique_count = 1
                interior(unique_count) = interior(i)
            elseif (interior(i) /= interior(unique_count)) then
                unique_count = unique_count + 1
                interior(unique_count) = interior(i)
            end if
        end do

        allocate(prepared(unique_count + 2))
        if (unique_count > 0) prepared(1:unique_count) = interior(1:unique_count)
        prepared(unique_count + 1:unique_count + 2) = 0.0_dp
    end subroutine prepare_break_points

end module scifort_integrate_quadpack
