! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Von Mises-Fisher distribution on the unit hypersphere. Density and entropy
! use the standard Bessel-I normalization. Sampling uses the circular inverse
! transform in 2D, the stable direct S^2 construction in 3D, and Wood's 1994
! rejection algorithm in dimensions four and higher.

module scifort_vonmises_fisher
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_bessel_i, only : besseli_ratio, log_besseli_scaled
    use scifort_beta, only : beta_ppf
    use scifort_constants, only : scifort_pi
    use scifort_kinds, only : dp
    use scifort_math, only : expm1_safe, log1p_safe, quiet_nan
    use scifort_random, only : rng_state, rng_uniform
    use scifort_uniform_direction, only : uniform_direction_rvs, &
        uniform_direction_status_success
    use scifort_vonmises, only : vonmises_ppf
    implicit none
    private

    integer, parameter, public :: vonmises_fisher_status_success = 0
    integer, parameter, public :: vonmises_fisher_status_invalid_shape = 1
    integer, parameter, public :: vonmises_fisher_status_invalid_parameter = 2
    integer, parameter, public :: vonmises_fisher_status_numerical_failure = 3

    real(dp), parameter :: unit_vector_tol = 1.1e-5_dp

    public :: vonmises_fisher_entropy
    public :: vonmises_fisher_fit
    public :: vonmises_fisher_logpdf
    public :: vonmises_fisher_pdf
    public :: vonmises_fisher_rvs
    public :: vonmises_fisher_rvs_array

contains

    function vonmises_fisher_logpdf(x, mu, kappa) result(y)
        real(dp), intent(in) :: x(:) !! unit vector at which the density is evaluated
        real(dp), intent(in) :: mu(:) !! mean direction; unit vector of matching dimension >= 2
        real(dp), intent(in) :: kappa !! positive concentration parameter
        real(dp) :: y

        integer :: dim
        real(dp) :: nu

        dim = size(mu)
        if (.not. valid_parameters(mu, kappa) .or. size(x) /= dim .or. &
                .not. all(ieee_is_finite(x)) .or. .not. is_unit_vector(x)) then
            y = quiet_nan(0.0_dp)
            return
        end if
        nu = 0.5_dp * real(dim,dp) - 1.0_dp
        y = nu * log(kappa) - 0.5_dp * real(dim,dp) * log(2.0_dp * scifort_pi) - &
            log_besseli_scaled(nu, kappa) - kappa + kappa * dot_product(mu, x)
    end function vonmises_fisher_logpdf

    function vonmises_fisher_pdf(x, mu, kappa) result(y)
        real(dp), intent(in) :: x(:) !! unit vector at which the density is evaluated
        real(dp), intent(in) :: mu(:) !! mean direction; unit vector of matching dimension >= 2
        real(dp), intent(in) :: kappa !! positive concentration parameter
        real(dp) :: y

        y = exp(vonmises_fisher_logpdf(x, mu, kappa))
    end function vonmises_fisher_pdf

    function vonmises_fisher_entropy(mu, kappa) result(y)
        real(dp), intent(in) :: mu(:) !! mean direction; unit vector of dimension >= 2
        real(dp), intent(in) :: kappa !! positive concentration parameter
        real(dp) :: y

        integer :: dim
        real(dp) :: log_norm
        real(dp) :: nu

        if (.not. valid_parameters(mu, kappa)) then
            y = quiet_nan(0.0_dp)
            return
        end if
        dim = size(mu)
        nu = 0.5_dp * real(dim,dp) - 1.0_dp
        log_norm = nu * log(kappa) - 0.5_dp * real(dim,dp) * log(2.0_dp * scifort_pi) - &
            log_besseli_scaled(nu, kappa) - kappa
        y = -log_norm - kappa * besseli_ratio(nu, kappa)
    end function vonmises_fisher_entropy

    subroutine vonmises_fisher_rvs(state, mu, kappa, sample, status)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by the draw
        real(dp), intent(in) :: mu(:) !! mean direction; unit vector of dimension >= 2
        real(dp), intent(in) :: kappa !! positive concentration parameter
        real(dp), intent(out) :: sample(:) !! sampled unit vector of the same dimension as mu
        integer, intent(out), optional :: status !! zero on success; nonzero for invalid/numerical failure

        integer :: dim
        integer :: info

        dim = size(mu)
        sample = 0.0_dp
        if (size(sample) /= dim .or. dim < 2) then
            if (present(status)) status = vonmises_fisher_status_invalid_shape
            return
        end if
        if (.not. valid_parameters(mu, kappa)) then
            if (present(status)) status = vonmises_fisher_status_invalid_parameter
            return
        end if

        select case (dim)
        case (2)
            call sample_2d(state, mu, kappa, sample)
            info = vonmises_fisher_status_success
        case (3)
            call sample_3d_north(state, kappa, sample, info)
            if (info == vonmises_fisher_status_success) call rotate_from_north(sample, mu)
        case default
            call sample_wood_north(state, kappa, sample, info)
            if (info == vonmises_fisher_status_success) call rotate_from_north(sample, mu)
        end select

        if (present(status)) status = info
    end subroutine vonmises_fisher_rvs

    subroutine vonmises_fisher_rvs_array(state, mu, kappa, samples, status)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by all draws
        real(dp), intent(in) :: mu(:) !! mean direction; unit vector of dimension >= 2
        real(dp), intent(in) :: kappa !! positive concentration parameter
        real(dp), intent(out) :: samples(:, :) !! samples by column, shape (dimension, count)
        integer, intent(out), optional :: status !! zero on success; nonzero for invalid/numerical failure

        integer :: info
        integer :: j

        samples = 0.0_dp
        if (size(samples,1) /= size(mu) .or. size(mu) < 2) then
            if (present(status)) status = vonmises_fisher_status_invalid_shape
            return
        end if
        do j = 1, size(samples,2)
            call vonmises_fisher_rvs(state, mu, kappa, samples(:,j), info)
            if (info /= vonmises_fisher_status_success) then
                if (present(status)) status = info
                return
            end if
        end do
        if (present(status)) status = vonmises_fisher_status_success
    end subroutine vonmises_fisher_rvs_array

    subroutine vonmises_fisher_fit(data, mu, kappa, status)
        real(dp), intent(in) :: data(:, :) !! unit-vector observations by column, shape (dimension, count)
        real(dp), intent(out) :: mu(:) !! fitted mean direction
        real(dp), intent(out) :: kappa !! fitted positive concentration
        integer, intent(out) :: status !! zero on success; nonzero for invalid/numerical failure

        integer :: dim
        integer :: j
        integer :: nobs
        real(dp) :: avg(size(data,1))
        real(dp) :: hi
        real(dp) :: lo
        real(dp) :: mid
        real(dp) :: norm_avg
        real(dp) :: nu
        real(dp) :: r
        real(dp) :: ratio_mid

        dim = size(data,1)
        nobs = size(data,2)
        mu = quiet_nan(0.0_dp)
        kappa = quiet_nan(0.0_dp)
        if (dim < 2 .or. nobs < 1 .or. size(mu) /= dim .or. &
                .not. all(ieee_is_finite(data))) then
            status = vonmises_fisher_status_invalid_shape
            return
        end if
        do j = 1, nobs
            if (.not. is_unit_vector(data(:,j))) then
                status = vonmises_fisher_status_invalid_parameter
                return
            end if
        end do

        avg = sum(data, dim=2) / real(nobs,dp)
        norm_avg = sqrt(dot_product(avg, avg))
        if (.not. ieee_is_finite(norm_avg) .or. norm_avg <= sqrt(epsilon(1.0_dp))) then
            status = vonmises_fisher_status_numerical_failure
            return
        end if
        mu = avg / norm_avg
        r = min(norm_avg, 1.0_dp)
        nu = 0.5_dp * real(dim,dp) - 1.0_dp

        lo = 1.0e-8_dp
        hi = 1.0e9_dp
        if (r <= besseli_ratio(nu, lo) .or. r >= besseli_ratio(nu, hi)) then
            status = vonmises_fisher_status_numerical_failure
            return
        end if
        do j = 1, 120
            mid = sqrt(lo * hi)
            ratio_mid = besseli_ratio(nu, mid)
            if (ratio_mid < r) then
                lo = mid
            else
                hi = mid
            end if
            if (abs(hi - lo) <= 2.0e-12_dp * max(1.0_dp, mid)) exit
        end do
        kappa = 0.5_dp * (lo + hi)
        status = vonmises_fisher_status_success
    end subroutine vonmises_fisher_fit

    subroutine sample_2d(state, mu, kappa, sample)
        type(rng_state), intent(inout) :: state !! RNG state advanced by one circular draw
        real(dp), intent(in) :: mu(:) !! two-dimensional unit mean direction
        real(dp), intent(in) :: kappa !! positive concentration parameter
        real(dp), intent(out) :: sample(:) !! two-dimensional unit sample

        real(dp) :: angle
        real(dp) :: mean_angle

        mean_angle = atan2(mu(2), mu(1))
        angle = vonmises_ppf(rng_uniform(state), kappa, mean_angle, 1.0_dp)
        sample(1) = cos(angle)
        sample(2) = sin(angle)
    end subroutine sample_2d

    subroutine sample_3d_north(state, kappa, sample, status)
        type(rng_state), intent(inout) :: state !! RNG state advanced by the axial and circle draws
        real(dp), intent(in) :: kappa !! positive concentration parameter
        real(dp), intent(out) :: sample(:) !! three-dimensional sample around the first coordinate axis
        integer, intent(out) :: status !! zero on success; nonzero on numerical failure

        integer :: info
        real(dp) :: a
        real(dp) :: circle(2)
        real(dp) :: radial
        real(dp) :: u
        real(dp) :: x

        u = rng_uniform(state)
        a = -expm1_safe(-2.0_dp * kappa)
        x = 1.0_dp + log1p_safe(-(1.0_dp - u) * a) / kappa
        x = min(1.0_dp, max(-1.0_dp, x))
        radial = sqrt(max(0.0_dp, 1.0_dp - x * x))
        call uniform_direction_rvs(state, circle, info)
        if (info /= uniform_direction_status_success) then
            status = vonmises_fisher_status_numerical_failure
            return
        end if
        sample = [x, radial * circle(1), radial * circle(2)]
        status = vonmises_fisher_status_success
    end subroutine sample_3d_north

    subroutine sample_wood_north(state, kappa, sample, status)
        type(rng_state), intent(inout) :: state !! RNG state advanced by rejection and direction draws
        real(dp), intent(in) :: kappa !! positive concentration parameter
        real(dp), intent(out) :: sample(:) !! sample around the first coordinate axis, dimension >= 4
        integer, intent(out) :: status !! zero on success; nonzero on numerical failure

        integer, parameter :: max_attempts = 100000
        integer :: attempt
        integer :: dim
        integer :: info
        real(dp) :: b
        real(dp) :: beta_draw
        real(dp) :: correction
        real(dp) :: criterion
        real(dp) :: d1
        real(dp) :: halfdim
        real(dp) :: node
        real(dp) :: radial
        real(dp) :: root
        real(dp) :: u
        real(dp) :: x
        real(dp), allocatable :: rest(:)

        dim = size(sample)
        d1 = real(dim - 1, dp)
        halfdim = 0.5_dp * d1
        root = hypot(2.0_dp * kappa, d1)
        b = d1 / (root + 2.0_dp * kappa)
        node = (1.0_dp - b) / (1.0_dp + b)
        correction = kappa * node + d1 * (log(4.0_dp) + log(b) - 2.0_dp * log1p_safe(b))

        x = 0.0_dp
        do attempt = 1, max_attempts
            beta_draw = beta_ppf(rng_uniform(state), halfdim, halfdim)
            x = (1.0_dp - (1.0_dp + b) * beta_draw) / &
                (1.0_dp - (1.0_dp - b) * beta_draw)
            u = rng_uniform(state)
            criterion = kappa * x + d1 * log((1.0_dp + b - x + x * b) / (1.0_dp + b)) - correction
            if (criterion > log(u)) exit
        end do
        if (attempt > max_attempts) then
            status = vonmises_fisher_status_numerical_failure
            return
        end if

        allocate(rest(dim - 1))
        call uniform_direction_rvs(state, rest, info)
        if (info /= uniform_direction_status_success) then
            status = vonmises_fisher_status_numerical_failure
            return
        end if
        radial = sqrt(max(0.0_dp, 1.0_dp - x * x))
        sample(1) = x
        sample(2:dim) = radial * rest
        status = vonmises_fisher_status_success
    end subroutine sample_wood_north

    pure subroutine rotate_from_north(sample, mu)
        real(dp), intent(inout) :: sample(:) !! sample initially centered on the first coordinate axis
        real(dp), intent(in) :: mu(:) !! target unit mean direction

        real(dp) :: denom
        real(dp) :: projection
        real(dp) :: v(size(mu))

        v = -mu
        v(1) = v(1) + 1.0_dp
        denom = dot_product(v, v)
        if (denom <= 64.0_dp * epsilon(1.0_dp)) return
        projection = dot_product(v, sample)
        sample = sample - (2.0_dp * projection / denom) * v
    end subroutine rotate_from_north

    pure logical function valid_parameters(mu, kappa) result(valid)
        real(dp), intent(in) :: mu(:) !! candidate mean direction
        real(dp), intent(in) :: kappa !! candidate concentration

        valid = size(mu) >= 2 .and. ieee_is_finite(kappa) .and. kappa > 0.0_dp .and. &
            all(ieee_is_finite(mu))
        if (valid) valid = is_unit_vector(mu)
    end function valid_parameters

    pure logical function is_unit_vector(x) result(valid)
        real(dp), intent(in) :: x(:) !! candidate direction vector

        real(dp) :: normx

        if (size(x) < 1 .or. .not. all(ieee_is_finite(x))) then
            valid = .false.
            return
        end if
        normx = sqrt(dot_product(x, x))
        valid = abs(normx - 1.0_dp) <= unit_vector_tol
    end function is_unit_vector

end module scifort_vonmises_fisher
