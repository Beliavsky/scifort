! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
module scifort_gaussian_kde_bridge
    use, intrinsic :: iso_c_binding, only : c_ptr, c_loc, c_f_pointer, c_associated, c_null_ptr
    use scifort_kinds, only : dp
    use scifort_math, only : quiet_nan
    use scifort_gaussian_kde, only : gaussian_kde, gaussian_kde_init
    implicit none
    private
    public :: kde_create, kde_destroy, kde_pdf, kde_logpdf, kde_info
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
end module scifort_gaussian_kde_bridge
