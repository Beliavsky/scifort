! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
module scifort_gaussian_kde_bridge
    use, intrinsic :: iso_c_binding, only : c_ptr, c_loc, c_f_pointer, c_associated, c_null_ptr
    use scifort_kinds, only : dp
    use scifort_math, only : quiet_nan
    use scifort_random, only : rng_state
    use scifort_gaussian_kde, only : gaussian_kde, gaussian_kde_init
    implicit none
    private
    public :: kde_create, kde_destroy, kde_pdf, kde_logpdf, kde_info
    public :: kde_integrate_gaussian, kde_integrate_box_1d, kde_integrate_box
    public :: kde_integrate_kde, kde_resample, kde_marginal
contains
    subroutine kde_create(dataset, weights, bandwidth, factor, handle, status)
        !! Copy observations into a heap-owned estimator; failure returns a null handle.
        !! @dimension d = size(dataset,1)
        !! @dimension n = size(dataset,2)
        !! @require d >= 1
        !! @require n >= 2
        !! d <= n; weighted covariance must be positive definite.
        real(dp), intent(in) :: dataset(:,:) !! shape (d,n); finite observations by column
        real(dp), intent(in) :: weights(:) !! shape (n); finite nonnegative weights with positive sum
        integer, intent(in) :: bandwidth !! 0 Scott; 1 Silverman; 2 constant factor
        real(dp), intent(in) :: factor !! positive finite factor for mode 2; ignored otherwise
        type(c_ptr), intent(out) :: handle !! owning opaque handle; null on failure
        integer, intent(out) :: status !! 0 success; 1-4 native error; 6 object allocation failure
        type(gaussian_kde), pointer :: kde
        integer :: info

        handle = c_null_ptr
        status = 2
        if (bandwidth < 0 .or. bandwidth > 2) return
        allocate(kde, stat=info)
        if (info /= 0) then
            status = 6
            return
        end if
        select case (bandwidth)
        case (0)
            call gaussian_kde_init(kde, dataset, status, weights=weights, bw_method='scott')
        case (1)
            call gaussian_kde_init(kde, dataset, status, weights=weights, bw_method='silverman')
        case (2)
            call gaussian_kde_init(kde, dataset, status, weights=weights, bw_factor=factor)
        end select
        if (status /= 0) then
            deallocate(kde)
            return
        end if
        handle = c_loc(kde)
    end subroutine kde_create

    subroutine kde_destroy(handle, status)
        !! Release this owning handle and null its slot; a null slot is already closed.
        type(c_ptr), intent(inout) :: handle !! owning handle slot; aliases become invalid after release
        integer, intent(out) :: status !! 0 success, including a null handle
        type(gaussian_kde), pointer :: kde

        if (c_associated(handle)) then
            call c_f_pointer(handle, kde)
            deallocate(kde)
            handle = c_null_ptr
        end if
        status = 0
    end subroutine kde_destroy

    subroutine kde_pdf(handle, points, values, status)
        !! Evaluate density in a batch; native failures fill values with NaNs.
        !! @dimension d = size(points,1)
        !! @dimension m = size(points,2)
        !! @require d >= 1
        !! @require m >= 0
        !! @when m == 0: return an empty result after handle and dimension validation
        !! d must also equal the estimator dimension.
        type(c_ptr), intent(in) :: handle !! live estimator handle; null returns status 5
        real(dp), intent(in) :: points(:,:) !! shape (d,m); finite evaluation points by column
        real(dp), intent(out) :: values(:) !! shape (m); densities, or NaNs on failure
        integer, intent(out) :: status !! 0 success; 1 shape/nonfinite point; 4 solve error; 5 null handle
        type(gaussian_kde), pointer :: kde
        real(dp), allocatable :: evaluated(:)

        values = quiet_nan(0.0_dp)
        status = 5
        if (.not. c_associated(handle)) return
        call c_f_pointer(handle, kde)
        status = 1
        if (size(points,1) /= kde%d .or. size(values) /= size(points,2)) return
        evaluated = kde%pdf(points, status)
        if (status == 0) values = evaluated
    end subroutine kde_pdf

    subroutine kde_logpdf(handle, points, values, status)
        !! Evaluate log density directly through the native stable batch method.
        !! @dimension d = size(points,1)
        !! @dimension m = size(points,2)
        !! @require d >= 1
        !! @require m >= 0
        !! @when m == 0: return an empty result after handle and dimension validation
        !! d must also equal the estimator dimension.
        type(c_ptr), intent(in) :: handle !! live estimator handle; null returns status 5
        real(dp), intent(in) :: points(:,:) !! shape (d,m); finite evaluation points by column
        real(dp), intent(out) :: values(:) !! shape (m); log densities, or NaNs on failure
        integer, intent(out) :: status !! 0 success; 1 shape/nonfinite point; 4 solve error; 5 null handle
        type(gaussian_kde), pointer :: kde
        real(dp), allocatable :: evaluated(:)

        values = quiet_nan(0.0_dp)
        status = 5
        if (.not. c_associated(handle)) return
        call c_f_pointer(handle, kde)
        status = 1
        if (size(points,1) /= kde%d .or. size(values) /= size(points,2)) return
        evaluated = kde%logpdf(points, status)
        if (status == 0) values = evaluated
    end subroutine kde_logpdf

    subroutine kde_info(handle, d, n, neff, factor, scott, silverman, status)
        !! Copy dimension, sample count and bandwidth information from a live estimator.
        type(c_ptr), intent(in) :: handle !! live estimator handle; null returns status 5
        integer, intent(out) :: d !! data dimension; zero on failure
        integer, intent(out) :: n !! observation count; zero on failure
        real(dp), intent(out) :: neff !! effective sample size; NaN on failure
        real(dp), intent(out) :: factor !! current bandwidth factor; NaN on failure
        real(dp), intent(out) :: scott !! Scott bandwidth factor; NaN on failure
        real(dp), intent(out) :: silverman !! Silverman bandwidth factor; NaN on failure
        integer, intent(out) :: status !! 0 success; 5 null handle
        type(gaussian_kde), pointer :: kde

        d = 0
        n = 0
        neff = quiet_nan(0.0_dp)
        factor = neff
        scott = neff
        silverman = neff
        status = 5
        if (.not. c_associated(handle)) return
        call c_f_pointer(handle, kde)
        d = kde%d
        n = kde%n
        neff = kde%neff
        factor = kde%factor
        scott = kde%scotts_factor()
        silverman = kde%silverman_factor()
        status = 0
    end subroutine kde_info

    subroutine kde_integrate_gaussian(handle, mean, cov, value, status)
        !! Integrate the estimator times a Gaussian density with matching dimension.
        !! @dimension d = size(mean)
        !! @require d >= 1
        !! d must equal the estimator dimension; cov plus kernel covariance must be SPD.
        type(c_ptr), intent(in) :: handle !! live estimator handle; null returns status 5
        real(dp), intent(in) :: mean(:) !! shape (d); finite Gaussian mean
        real(dp), intent(in) :: cov(:,:) !! shape (d,d); finite Gaussian covariance
        real(dp), intent(out) :: value !! integral; NaN on failure
        integer, intent(out) :: status !! 0 success; 1 invalid shape/data; 4 solve error; 5 null handle
        type(gaussian_kde), pointer :: kde
        real(dp) :: evaluated

        value = quiet_nan(0.0_dp)
        status = 5
        if (.not. c_associated(handle)) return
        call c_f_pointer(handle, kde)
        status = 1
        if (size(mean) /= kde%d .or. size(cov,1) /= kde%d .or. size(cov,2) /= kde%d) return
        evaluated = kde%integrate_gaussian(mean, cov, status)
        if (status == 0) value = evaluated
    end subroutine kde_integrate_gaussian

    subroutine kde_integrate_box_1d(handle, low, high, value, status)
        !! Integrate over a one-dimensional interval, preserving native signed limits.
        type(c_ptr), intent(in) :: handle !! live univariate estimator; null returns status 5
        real(dp), intent(in) :: low !! lower limit; may be infinite or exceed high
        real(dp), intent(in) :: high !! upper limit; may be infinite
        real(dp), intent(out) :: value !! integral; native NaN retained; NaN on failure
        integer, intent(out) :: status !! 0 native success; 1 wrong dimension; 5 null handle
        type(gaussian_kde), pointer :: kde
        real(dp) :: evaluated

        value = quiet_nan(0.0_dp)
        status = 5
        if (.not. c_associated(handle)) return
        call c_f_pointer(handle, kde)
        evaluated = kde%integrate_box_1d(low, high, status)
        if (status == 0) value = evaluated
    end subroutine kde_integrate_box_1d

    subroutine kde_integrate_box(handle, low, high, maxpts, state, value, status)
        !! Integrate over a box using caller-owned randomized-QMC state.
        !! @dimension d = size(low)
        !! @require d >= 1
        !! d must equal the estimator dimension; reversed bounds retain native signs.
        type(c_ptr), intent(in) :: handle !! live estimator handle; null returns status 5
        real(dp), intent(in) :: low(:) !! shape (d); lower limits, infinite permitted
        real(dp), intent(in) :: high(:) !! shape (d); upper limits, infinite permitted
        integer, intent(in) :: maxpts !! approximate QMC points per component; 0 native default, >= 0
        type(rng_state), intent(inout) :: state !! explicit QMC state, including shifts consumed on failure
        real(dp), intent(out) :: value !! integral; NaN on failure
        integer, intent(out) :: status !! 0 success; 1 shape; 2 budget; 4 integration error; 5 null handle
        type(gaussian_kde), pointer :: kde
        real(dp) :: evaluated

        value = quiet_nan(0.0_dp)
        status = 5
        if (.not. c_associated(handle)) return
        call c_f_pointer(handle, kde)
        status = 1
        if (size(low) /= kde%d .or. size(high) /= kde%d) return
        status = 2
        if (maxpts < 0) return
        if (maxpts == 0) then
            evaluated = kde%integrate_box(low, high, state=state, status=status)
        else
            evaluated = kde%integrate_box(low, high, maxpts=maxpts, state=state, status=status)
        end if
        if (status == 0) value = evaluated
    end subroutine kde_integrate_box

    subroutine kde_integrate_kde(handle, other, value, status)
        !! Integrate the product of two estimators with equal dimension.
        type(c_ptr), intent(in) :: handle !! first live estimator; may equal other
        type(c_ptr), intent(in) :: other !! second live estimator; null returns status 5
        real(dp), intent(out) :: value !! integral of the product; NaN on failure
        integer, intent(out) :: status !! 0 success; 1 dimension mismatch; 4 solve error; 5 null handle
        type(gaussian_kde), pointer :: kde, second
        real(dp) :: evaluated

        value = quiet_nan(0.0_dp)
        status = 5
        if (.not. c_associated(handle) .or. .not. c_associated(other)) return
        call c_f_pointer(handle, kde)
        call c_f_pointer(other, second)
        evaluated = kde%integrate_kde(second, status)
        if (status == 0) value = evaluated
    end subroutine kde_integrate_kde

    subroutine kde_resample(handle, state, samples, status)
        !! Draw samples in native column order with explicit RNG state.
        !! @dimension d = size(samples,1)
        !! @dimension n_samples = size(samples,2)
        !! @require d >= 1
        !! @require n_samples >= 0
        !! @when n_samples == 0: validate handle and dimension; leave RNG state unchanged
        !! d must equal the estimator dimension.
        type(c_ptr), intent(in) :: handle !! live estimator handle; null returns status 5
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by all draws
        real(dp), intent(out) :: samples(:,:) !! shape (d,n_samples); samples by column, NaNs on failure
        integer, intent(out) :: status !! 0 success; 1 dimension mismatch; 5 null handle
        type(gaussian_kde), pointer :: kde

        samples = quiet_nan(0.0_dp)
        status = 5
        if (.not. c_associated(handle)) return
        call c_f_pointer(handle, kde)
        status = 1
        if (size(samples,1) /= kde%d) return
        call kde%resample(state, samples, status)
    end subroutine kde_resample

    subroutine kde_marginal(handle, dimensions, marginal_handle, status)
        !! Return an independently owned marginal; preserve weights and bandwidth factor.
        !! @dimension k = size(dimensions)
        !! @require k >= 1
        !! k <= estimator dimension; coordinates are distinct and one-based at this boundary.
        type(c_ptr), intent(in) :: handle !! parent estimator; parent remains unchanged
        integer, intent(in) :: dimensions(:) !! shape (k); distinct one-based coordinates
        type(c_ptr), intent(out) :: marginal_handle !! new owning handle; null on native failure
        integer, intent(out) :: status !! 0 success; 2 invalid coordinates; 5 null parent; 6 allocation failure
        type(gaussian_kde), pointer :: kde, child
        integer :: info

        marginal_handle = c_null_ptr
        status = 5
        if (.not. c_associated(handle)) return
        call c_f_pointer(handle, kde)
        status = 2
        if (size(dimensions) < 1 .or. size(dimensions) > kde%d) return
        if (any(dimensions < 1) .or. any(dimensions > kde%d)) return
        allocate(child, stat=info)
        if (info /= 0) then
            status = 6
            return
        end if
        child = kde%marginal(dimensions, status)
        if (status /= 0) then
            deallocate(child)
            return
        end if
        marginal_handle = c_loc(child)
    end subroutine kde_marginal
end module scifort_gaussian_kde_bridge
