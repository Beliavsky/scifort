! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
module scifort_ecdf_bridge
    use, intrinsic :: iso_c_binding, only : c_ptr, c_null_ptr, c_loc, c_f_pointer, c_associated
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_math, only : quiet_nan
    use scifort_empirical_distribution, only : ecdf, ecdf_result, &
        empirical_distribution_function, ecdf_confidence_interval
    implicit none
    private
    public :: ecdf_create, ecdf_destroy, ecdf_info, ecdf_evaluate, ecdf_values
    public :: ecdf_interval
contains
    subroutine ecdf_create(uncensored, right, cdf_handle, sf_handle, status)
        !! Copy ordinary or right-censored observations into independent CDF/SF objects.
        !! @dimension n = size(uncensored)
        !! @dimension r = size(right)
        !! @require n >= 0
        !! @require r >= 0
        !! @when n == 0: no observed events; an all-censored or empty sample is permitted
        !! @when r == 0: ordinary ECDF; the two input lengths are independent
        real(dp), intent(in) :: uncensored(:) !! shape (n); finite event times, including ties
        real(dp), intent(in) :: right(:) !! shape (r); finite right-censoring times
        type(c_ptr), intent(out) :: cdf_handle !! new owning CDF; null on failure
        type(c_ptr), intent(out) :: sf_handle !! new owning SF; null on failure
        integer, intent(out) :: status !! 0 success; 2 nonfinite input; 6 object allocation failure
        type(ecdf_result) :: result

        cdf_handle = c_null_ptr
        sf_handle = c_null_ptr
        status = 2
        if (.not. all(ieee_is_finite(uncensored)) .or. .not. all(ieee_is_finite(right))) return
        if (size(right) == 0) then
            result = ecdf(uncensored)
        else
            result = ecdf(uncensored, right)
        end if
        call own_pair(result%cdf, result%sf, cdf_handle, sf_handle, status)
    end subroutine ecdf_create

    subroutine ecdf_destroy(handle, status)
        !! Release an owned distribution function and null its slot; null is harmless.
        type(c_ptr), intent(inout) :: handle !! owning slot; aliases become invalid after release
        integer, intent(out) :: status !! 0 success, including an already-null slot
        type(empirical_distribution_function), pointer :: edf

        if (c_associated(handle)) then
            call c_f_pointer(handle, edf)
            deallocate(edf)
            handle = c_null_ptr
        end if
        status = 0
    end subroutine ecdf_destroy

    subroutine ecdf_info(handle, count, variance, survival, status)
        !! Query exact array length and flags without exposing native descriptors.
        type(c_ptr), intent(in) :: handle !! live EDF; null returns status 5
        integer, intent(out) :: count !! number of distinct sample times; zero on failure
        integer, intent(out) :: variance !! 1 when confidence intervals are supported, otherwise 0
        integer, intent(out) :: survival !! 1 for SF, 0 for CDF; zero on failure
        integer, intent(out) :: status !! 0 success; 5 null handle
        type(empirical_distribution_function), pointer :: edf

        count = 0
        variance = 0
        survival = 0
        status = 5
        if (.not. c_associated(handle)) return
        call c_f_pointer(handle, edf)
        count = size(edf%quantiles)
        if (edf%has_variance_data) variance = 1
        if (edf%is_survival) survival = 1
        status = 0
    end subroutine ecdf_info

    subroutine ecdf_evaluate(handle, x, values, status)
        !! Evaluate right-continuous steps; preserve native NaN/infinity/empty semantics.
        !! @dimension n = size(x)
        !! @require n >= 0
        !! @when n == 0: validate the handle and return an empty result
        type(c_ptr), intent(in) :: handle !! live EDF; null returns status 5
        real(dp), intent(in) :: x(:) !! shape (n); evaluation points; NaN and infinities permitted
        real(dp), intent(out) :: values(:) !! shape (n); evaluated probabilities; NaNs on failure
        integer, intent(out) :: status !! 0 success; 1 shape; 5 null handle
        type(empirical_distribution_function), pointer :: edf

        values = quiet_nan(0.0_dp)
        status = 5
        if (.not. c_associated(handle)) return
        call c_f_pointer(handle, edf)
        status = 1
        if (size(values) /= size(x)) return
        values = edf%evaluate(x)
        status = 0
    end subroutine ecdf_evaluate

    subroutine ecdf_values(handle, quantiles, probabilities, status)
        !! Copy stored knots and probabilities; m must equal the count from ecdf_info.
        !! @dimension m = size(quantiles)
        !! @require m >= 0
        !! @when m == 0: valid only for an empty EDF
        type(c_ptr), intent(in) :: handle !! live EDF; null returns status 5
        real(dp), intent(out) :: quantiles(:) !! shape (m); distinct sorted times; NaNs on failure
        real(dp), intent(out) :: probabilities(:) !! shape (m); probabilities at knots; NaNs on failure
        integer, intent(out) :: status !! 0 success; 1 incorrect length; 5 null handle
        type(empirical_distribution_function), pointer :: edf

        quantiles = quiet_nan(0.0_dp)
        probabilities = quiet_nan(0.0_dp)
        status = 5
        if (.not. c_associated(handle)) return
        call c_f_pointer(handle, edf)
        status = 1
        if (size(quantiles) /= size(edf%quantiles)) return
        if (size(probabilities) /= size(quantiles)) return
        quantiles = edf%quantiles
        probabilities = edf%probabilities
        status = 0
    end subroutine ecdf_values

    subroutine ecdf_interval(handle, level, method, low_handle, high_handle, status)
        !! Copy Greenwood confidence bounds into independent owned EDFs.
        type(c_ptr), intent(in) :: handle !! live EDF with variance data; null returns status 5
        real(dp), intent(in) :: level !! confidence level in [0,1], finite
        integer, intent(in) :: method !! 0 linear Greenwood; 1 exponential Greenwood (log-log)
        type(c_ptr), intent(out) :: low_handle !! new owning lower bound; null on failure
        type(c_ptr), intent(out) :: high_handle !! new owning upper bound; null on failure
        integer, intent(out) :: status !! 0 success; 2 invalid options/nested bounds; 5 null; 6 allocation
        type(empirical_distribution_function), pointer :: edf
        type(ecdf_confidence_interval) :: ci

        low_handle = c_null_ptr
        high_handle = c_null_ptr
        status = 5
        if (.not. c_associated(handle)) return
        call c_f_pointer(handle, edf)
        status = 2
        if (.not. edf%has_variance_data) return
        if (.not. ieee_is_finite(level)) return
        if (level < 0.0_dp .or. level > 1.0_dp .or. method < 0 .or. method > 1) return
        if (method == 0) then
            ci = edf%confidence_interval(level, 'linear')
        else
            ci = edf%confidence_interval(level, 'log-log')
        end if
        call own_pair(ci%low, ci%high, low_handle, high_handle, status)
    end subroutine ecdf_interval

    subroutine own_pair(first, second, first_handle, second_handle, status)
        !! Allocate both owners before publishing either address.
        type(empirical_distribution_function), intent(in) :: first !! first EDF to copy
        type(empirical_distribution_function), intent(in) :: second !! second EDF to copy
        type(c_ptr), intent(out) :: first_handle !! new owning copy; null on failure
        type(c_ptr), intent(out) :: second_handle !! new owning copy; null on failure
        integer, intent(out) :: status !! 0 success; 6 outer object allocation failure
        type(empirical_distribution_function), pointer :: a, b
        integer :: info

        first_handle = c_null_ptr
        second_handle = c_null_ptr
        status = 6
        allocate(a, stat=info)
        if (info /= 0) return
        allocate(b, stat=info)
        if (info /= 0) then
            deallocate(a)
            return
        end if
        a = first
        b = second
        first_handle = c_loc(a)
        second_handle = c_loc(b)
        status = 0
    end subroutine own_pair
end module scifort_ecdf_bridge
