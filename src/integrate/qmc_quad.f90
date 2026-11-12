! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Quasi-Monte Carlo integration over finite hyperrectangles. Independent
! scrambled Halton replicates provide the integral estimate and standard error.
! Randomness is supplied only through an explicit SciFort RNG state.

module scifort_integrate_qmc_quad
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_integrate_quadrature, only : integrate_status_invalid_argument, &
        integrate_status_invalid_shape, integrate_status_numerical_failure, &
        integrate_status_success
    use scifort_kinds, only : dp
    use scifort_math, only : quiet_nan
    use scifort_qmc, only : qmc_halton, qmc_halton_init, qmc_status_success
    use scifort_random, only : rng_state
    implicit none
    private

    public :: qmc_quad

    abstract interface
        function qmc_integrand(x) result(value)
            import :: dp
            implicit none
            real(dp), intent(in) :: x(:) !! one point in the integration hyperrectangle
            real(dp) :: value
        end function qmc_integrand
    end interface

contains

    subroutine qmc_quad(f, a, b, result, standard_error, state, n_estimates, n_points, &
            neval, status)
        procedure(qmc_integrand) :: f !! scalar integrand callback for one multidimensional point
        real(dp), intent(in) :: a(:) !! lower integration bounds by dimension
        real(dp), intent(in) :: b(:) !! upper integration bounds by dimension; reversal is allowed
        real(dp), intent(out) :: result !! mean of independently scrambled QMC integral estimates
        real(dp), intent(out) :: standard_error !! standard error across independent QMC estimates
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced to scramble each estimate
        integer, intent(in), optional :: n_estimates !! independent scrambled estimates; default 8
        integer, intent(in), optional :: n_points !! Halton points per estimate; default 1024
        integer, intent(out), optional :: neval !! total integrand evaluations
        integer, intent(out), optional :: status !! integrate_status_* result

        integer :: d
        integer :: estimate_index
        integer :: i
        integer :: j
        integer :: qmc_status
        integer :: selected_estimates
        integer :: selected_points
        real(dp) :: volume
        real(dp) :: value
        real(dp), allocatable :: estimates(:)
        real(dp), allocatable :: point(:)
        real(dp), allocatable :: sample(:, :)
        type(qmc_halton) :: engine

        result = quiet_nan(0.0_dp)
        standard_error = quiet_nan(0.0_dp)
        if (present(neval)) neval = 0
        if (present(status)) status = integrate_status_invalid_argument

        d = size(a)
        if (d < 1 .or. size(b) /= d) then
            if (present(status)) status = integrate_status_invalid_shape
            return
        end if
        if (.not. all(ieee_is_finite(a)) .or. .not. all(ieee_is_finite(b))) return

        selected_estimates = 8
        selected_points = 1024
        if (present(n_estimates)) selected_estimates = n_estimates
        if (present(n_points)) selected_points = n_points
        if (selected_estimates < 1 .or. selected_points < 1) return

        volume = product(b - a)
        if (volume == 0.0_dp) then
            result = 0.0_dp
            standard_error = 0.0_dp
            if (present(status)) status = integrate_status_success
            return
        end if

        allocate(estimates(selected_estimates), point(d), sample(selected_points,d))
        do estimate_index = 1, selected_estimates
            call qmc_halton_init(engine, d, qmc_status, scramble=.true., state=state)
            if (qmc_status /= qmc_status_success) then
                if (present(status)) status = integrate_status_numerical_failure
                return
            end if
            call engine%random(selected_points, sample, qmc_status)
            if (qmc_status /= qmc_status_success) then
                if (present(status)) status = integrate_status_numerical_failure
                return
            end if

            estimates(estimate_index) = 0.0_dp
            do i = 1, selected_points
                do j = 1, d
                    point(j) = a(j) + sample(i,j) * (b(j) - a(j))
                end do
                value = f(point)
                if (.not. ieee_is_finite(value)) then
                    if (present(status)) status = integrate_status_numerical_failure
                    return
                end if
                estimates(estimate_index) = estimates(estimate_index) + value
            end do
            estimates(estimate_index) = volume * estimates(estimate_index) / real(selected_points, dp)
        end do

        result = sum(estimates) / real(selected_estimates, dp)
        if (selected_estimates == 1) then
            standard_error = quiet_nan(0.0_dp)
        else
            standard_error = sqrt(sum((estimates - result)**2) / &
                real(selected_estimates * (selected_estimates - 1), dp))
        end if
        if (present(neval)) neval = selected_estimates * selected_points
        if (present(status)) status = integrate_status_success
    end subroutine qmc_quad

end module scifort_integrate_qmc_quad
