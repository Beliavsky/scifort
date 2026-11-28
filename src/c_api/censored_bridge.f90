! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
module scifort_censored_bridge
    !! Flat, copied-data adapters for the existing native censored likelihoods.
    use, intrinsic :: iso_fortran_env, only : int64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_math, only : quiet_nan, valid_loc_scale
    use scifort_censored, only : censored_data, make_censored_data, &
        native_normal_loglikelihood => normal_censored_loglikelihood, &
        native_normal_nnlf => normal_censored_nnlf, &
        native_exponential_loglikelihood => exponential_censored_loglikelihood, &
        native_exponential_nnlf => exponential_censored_nnlf, &
        native_weibull_loglikelihood => weibull_censored_loglikelihood, &
        native_weibull_nnlf => weibull_censored_nnlf, &
        native_lognormal_loglikelihood => lognormal_censored_loglikelihood, &
        native_lognormal_nnlf => lognormal_censored_nnlf, &
        native_gamma_loglikelihood => gamma_censored_loglikelihood, &
        native_gamma_nnlf => gamma_censored_nnlf, &
        native_beta_loglikelihood => beta_censored_loglikelihood, &
        native_beta_nnlf => beta_censored_nnlf
    implicit none
    private
    public :: censored_normal_loglikelihood
    public :: censored_normal_nnlf
    public :: censored_exponential_loglikelihood
    public :: censored_exponential_nnlf
    public :: censored_weibull_loglikelihood
    public :: censored_weibull_nnlf
    public :: censored_lognormal_loglikelihood
    public :: censored_lognormal_nnlf
    public :: censored_gamma_loglikelihood
    public :: censored_gamma_nnlf
    public :: censored_beta_loglikelihood
    public :: censored_beta_nnlf
contains
    subroutine censored_normal_loglikelihood(uncensored, left, right, lower, upper, &
            loc, scale, value, status)
        !! Evaluate normal log likelihood using native endpoint/tail conventions.
        !! @dimension nu = size(uncensored)
        !! @dimension nl = size(left)
        !! @dimension nr = size(right)
        !! @dimension ni = size(lower)
        !! size(upper) = ni; corresponding endpoints must be paired.
        !! @require nu >= 0
        !! @require nl >= 0
        !! @require nr >= 0
        !! @require ni >= 0
        !! All group lengths are independent; empty groups are permitted.
        !! Their sum must fit the native default integer kind.
        real(dp), intent(in) :: uncensored(:) !! shape (nu); finite exact observations
        real(dp), intent(in) :: left(:) !! shape (nl); finite upper censoring limits
        real(dp), intent(in) :: right(:) !! shape (nr); finite lower censoring limits
        real(dp), intent(in) :: lower(:) !! shape (ni); interval lower endpoints
        real(dp), intent(in) :: upper(:) !! shape (ni); interval upper endpoints
        real(dp), intent(in) :: loc !! finite mean
        real(dp), intent(in) :: scale !! finite positive standard deviation
        real(dp), intent(out) :: value !! log likelihood; NaN for invalid input
        integer, intent(out) :: status !! 0 evaluated; 1 count/shape; 2 data/parameters
        type(censored_data) :: data

        value = quiet_nan(0.0_dp)
        status = 2
        if (.not. valid_loc_scale(loc,scale)) return
        call pack_censored_samples(uncensored,left,right,lower,upper,data,status)
        if (status /= 0) return
        value = native_normal_loglikelihood(data,loc, scale)
    end subroutine censored_normal_loglikelihood

    subroutine censored_normal_nnlf(uncensored, left, right, lower, upper, &
            loc, scale, value, status)
        !! Evaluate normal negative log likelihood using native endpoint/tail conventions.
        !! @dimension nu = size(uncensored)
        !! @dimension nl = size(left)
        !! @dimension nr = size(right)
        !! @dimension ni = size(lower)
        !! size(upper) = ni; corresponding endpoints must be paired.
        !! @require nu >= 0
        !! @require nl >= 0
        !! @require nr >= 0
        !! @require ni >= 0
        !! All group lengths are independent; empty groups are permitted.
        !! Their sum must fit the native default integer kind.
        real(dp), intent(in) :: uncensored(:) !! shape (nu); finite exact observations
        real(dp), intent(in) :: left(:) !! shape (nl); finite upper censoring limits
        real(dp), intent(in) :: right(:) !! shape (nr); finite lower censoring limits
        real(dp), intent(in) :: lower(:) !! shape (ni); interval lower endpoints
        real(dp), intent(in) :: upper(:) !! shape (ni); interval upper endpoints
        real(dp), intent(in) :: loc !! finite mean
        real(dp), intent(in) :: scale !! finite positive standard deviation
        real(dp), intent(out) :: value !! negative log likelihood; NaN for invalid input
        integer, intent(out) :: status !! 0 evaluated; 1 count/shape; 2 data/parameters
        type(censored_data) :: data

        value = quiet_nan(0.0_dp)
        status = 2
        if (.not. valid_loc_scale(loc,scale)) return
        call pack_censored_samples(uncensored,left,right,lower,upper,data,status)
        if (status /= 0) return
        value = native_normal_nnlf(data,loc, scale)
    end subroutine censored_normal_nnlf

    subroutine censored_exponential_loglikelihood(uncensored, left, right, lower, upper, &
            loc, scale, value, status)
        !! Evaluate exponential log likelihood using native endpoint/tail conventions.
        !! @dimension nu = size(uncensored)
        !! @dimension nl = size(left)
        !! @dimension nr = size(right)
        !! @dimension ni = size(lower)
        !! size(upper) = ni; corresponding endpoints must be paired.
        !! @require nu >= 0
        !! @require nl >= 0
        !! @require nr >= 0
        !! @require ni >= 0
        !! All group lengths are independent; empty groups are permitted.
        !! Their sum must fit the native default integer kind.
        real(dp), intent(in) :: uncensored(:) !! shape (nu); finite exact observations
        real(dp), intent(in) :: left(:) !! shape (nl); finite upper censoring limits
        real(dp), intent(in) :: right(:) !! shape (nr); finite lower censoring limits
        real(dp), intent(in) :: lower(:) !! shape (ni); interval lower endpoints
        real(dp), intent(in) :: upper(:) !! shape (ni); interval upper endpoints
        real(dp), intent(in) :: loc !! finite lower support endpoint
        real(dp), intent(in) :: scale !! finite positive scale
        real(dp), intent(out) :: value !! log likelihood; NaN for invalid input
        integer, intent(out) :: status !! 0 evaluated; 1 count/shape; 2 data/parameters
        type(censored_data) :: data

        value = quiet_nan(0.0_dp)
        status = 2
        if (.not. valid_loc_scale(loc,scale)) return
        call pack_censored_samples(uncensored,left,right,lower,upper,data,status)
        if (status /= 0) return
        value = native_exponential_loglikelihood(data,loc, scale)
    end subroutine censored_exponential_loglikelihood

    subroutine censored_exponential_nnlf(uncensored, left, right, lower, upper, &
            loc, scale, value, status)
        !! Evaluate exponential negative log likelihood using native endpoint/tail conventions.
        !! @dimension nu = size(uncensored)
        !! @dimension nl = size(left)
        !! @dimension nr = size(right)
        !! @dimension ni = size(lower)
        !! size(upper) = ni; corresponding endpoints must be paired.
        !! @require nu >= 0
        !! @require nl >= 0
        !! @require nr >= 0
        !! @require ni >= 0
        !! All group lengths are independent; empty groups are permitted.
        !! Their sum must fit the native default integer kind.
        real(dp), intent(in) :: uncensored(:) !! shape (nu); finite exact observations
        real(dp), intent(in) :: left(:) !! shape (nl); finite upper censoring limits
        real(dp), intent(in) :: right(:) !! shape (nr); finite lower censoring limits
        real(dp), intent(in) :: lower(:) !! shape (ni); interval lower endpoints
        real(dp), intent(in) :: upper(:) !! shape (ni); interval upper endpoints
        real(dp), intent(in) :: loc !! finite lower support endpoint
        real(dp), intent(in) :: scale !! finite positive scale
        real(dp), intent(out) :: value !! negative log likelihood; NaN for invalid input
        integer, intent(out) :: status !! 0 evaluated; 1 count/shape; 2 data/parameters
        type(censored_data) :: data

        value = quiet_nan(0.0_dp)
        status = 2
        if (.not. valid_loc_scale(loc,scale)) return
        call pack_censored_samples(uncensored,left,right,lower,upper,data,status)
        if (status /= 0) return
        value = native_exponential_nnlf(data,loc, scale)
    end subroutine censored_exponential_nnlf

    subroutine censored_weibull_loglikelihood(uncensored, left, right, lower, upper, &
            c, loc, scale, value, status)
        !! Evaluate weibull log likelihood using native endpoint/tail conventions.
        !! @dimension nu = size(uncensored)
        !! @dimension nl = size(left)
        !! @dimension nr = size(right)
        !! @dimension ni = size(lower)
        !! size(upper) = ni; corresponding endpoints must be paired.
        !! @require nu >= 0
        !! @require nl >= 0
        !! @require nr >= 0
        !! @require ni >= 0
        !! All group lengths are independent; empty groups are permitted.
        !! Their sum must fit the native default integer kind.
        real(dp), intent(in) :: uncensored(:) !! shape (nu); finite exact observations
        real(dp), intent(in) :: left(:) !! shape (nl); finite upper censoring limits
        real(dp), intent(in) :: right(:) !! shape (nr); finite lower censoring limits
        real(dp), intent(in) :: lower(:) !! shape (ni); interval lower endpoints
        real(dp), intent(in) :: upper(:) !! shape (ni); interval upper endpoints
        real(dp), intent(in) :: c !! finite positive Weibull shape
        real(dp), intent(in) :: loc !! finite lower support endpoint
        real(dp), intent(in) :: scale !! finite positive scale
        real(dp), intent(out) :: value !! log likelihood; NaN for invalid input
        integer, intent(out) :: status !! 0 evaluated; 1 count/shape; 2 data/parameters
        type(censored_data) :: data

        value = quiet_nan(0.0_dp)
        status = 2
        if (.not. valid_loc_scale(loc,scale)) return
        if (.not. ieee_is_finite(c)) return
        if (c <= 0.0_dp) return
        call pack_censored_samples(uncensored,left,right,lower,upper,data,status)
        if (status /= 0) return
        value = native_weibull_loglikelihood(data,c, loc, scale)
    end subroutine censored_weibull_loglikelihood

    subroutine censored_weibull_nnlf(uncensored, left, right, lower, upper, &
            c, loc, scale, value, status)
        !! Evaluate weibull negative log likelihood using native endpoint/tail conventions.
        !! @dimension nu = size(uncensored)
        !! @dimension nl = size(left)
        !! @dimension nr = size(right)
        !! @dimension ni = size(lower)
        !! size(upper) = ni; corresponding endpoints must be paired.
        !! @require nu >= 0
        !! @require nl >= 0
        !! @require nr >= 0
        !! @require ni >= 0
        !! All group lengths are independent; empty groups are permitted.
        !! Their sum must fit the native default integer kind.
        real(dp), intent(in) :: uncensored(:) !! shape (nu); finite exact observations
        real(dp), intent(in) :: left(:) !! shape (nl); finite upper censoring limits
        real(dp), intent(in) :: right(:) !! shape (nr); finite lower censoring limits
        real(dp), intent(in) :: lower(:) !! shape (ni); interval lower endpoints
        real(dp), intent(in) :: upper(:) !! shape (ni); interval upper endpoints
        real(dp), intent(in) :: c !! finite positive Weibull shape
        real(dp), intent(in) :: loc !! finite lower support endpoint
        real(dp), intent(in) :: scale !! finite positive scale
        real(dp), intent(out) :: value !! negative log likelihood; NaN for invalid input
        integer, intent(out) :: status !! 0 evaluated; 1 count/shape; 2 data/parameters
        type(censored_data) :: data

        value = quiet_nan(0.0_dp)
        status = 2
        if (.not. valid_loc_scale(loc,scale)) return
        if (.not. ieee_is_finite(c)) return
        if (c <= 0.0_dp) return
        call pack_censored_samples(uncensored,left,right,lower,upper,data,status)
        if (status /= 0) return
        value = native_weibull_nnlf(data,c, loc, scale)
    end subroutine censored_weibull_nnlf

    subroutine censored_lognormal_loglikelihood(uncensored, left, right, lower, upper, &
            s, loc, scale, value, status)
        !! Evaluate lognormal log likelihood using native endpoint/tail conventions.
        !! @dimension nu = size(uncensored)
        !! @dimension nl = size(left)
        !! @dimension nr = size(right)
        !! @dimension ni = size(lower)
        !! size(upper) = ni; corresponding endpoints must be paired.
        !! @require nu >= 0
        !! @require nl >= 0
        !! @require nr >= 0
        !! @require ni >= 0
        !! All group lengths are independent; empty groups are permitted.
        !! Their sum must fit the native default integer kind.
        real(dp), intent(in) :: uncensored(:) !! shape (nu); finite exact observations
        real(dp), intent(in) :: left(:) !! shape (nl); finite upper censoring limits
        real(dp), intent(in) :: right(:) !! shape (nr); finite lower censoring limits
        real(dp), intent(in) :: lower(:) !! shape (ni); interval lower endpoints
        real(dp), intent(in) :: upper(:) !! shape (ni); interval upper endpoints
        real(dp), intent(in) :: s !! finite positive standard deviation of log(x-loc)
        real(dp), intent(in) :: loc !! finite lower support endpoint
        real(dp), intent(in) :: scale !! finite positive exp(mean(log(x-loc)))
        real(dp), intent(out) :: value !! log likelihood; NaN for invalid input
        integer, intent(out) :: status !! 0 evaluated; 1 count/shape; 2 data/parameters
        type(censored_data) :: data

        value = quiet_nan(0.0_dp)
        status = 2
        if (.not. valid_loc_scale(loc,scale)) return
        if (.not. ieee_is_finite(s)) return
        if (s <= 0.0_dp) return
        call pack_censored_samples(uncensored,left,right,lower,upper,data,status)
        if (status /= 0) return
        value = native_lognormal_loglikelihood(data,s, loc, scale)
    end subroutine censored_lognormal_loglikelihood

    subroutine censored_lognormal_nnlf(uncensored, left, right, lower, upper, &
            s, loc, scale, value, status)
        !! Evaluate lognormal negative log likelihood using native endpoint/tail conventions.
        !! @dimension nu = size(uncensored)
        !! @dimension nl = size(left)
        !! @dimension nr = size(right)
        !! @dimension ni = size(lower)
        !! size(upper) = ni; corresponding endpoints must be paired.
        !! @require nu >= 0
        !! @require nl >= 0
        !! @require nr >= 0
        !! @require ni >= 0
        !! All group lengths are independent; empty groups are permitted.
        !! Their sum must fit the native default integer kind.
        real(dp), intent(in) :: uncensored(:) !! shape (nu); finite exact observations
        real(dp), intent(in) :: left(:) !! shape (nl); finite upper censoring limits
        real(dp), intent(in) :: right(:) !! shape (nr); finite lower censoring limits
        real(dp), intent(in) :: lower(:) !! shape (ni); interval lower endpoints
        real(dp), intent(in) :: upper(:) !! shape (ni); interval upper endpoints
        real(dp), intent(in) :: s !! finite positive standard deviation of log(x-loc)
        real(dp), intent(in) :: loc !! finite lower support endpoint
        real(dp), intent(in) :: scale !! finite positive exp(mean(log(x-loc)))
        real(dp), intent(out) :: value !! negative log likelihood; NaN for invalid input
        integer, intent(out) :: status !! 0 evaluated; 1 count/shape; 2 data/parameters
        type(censored_data) :: data

        value = quiet_nan(0.0_dp)
        status = 2
        if (.not. valid_loc_scale(loc,scale)) return
        if (.not. ieee_is_finite(s)) return
        if (s <= 0.0_dp) return
        call pack_censored_samples(uncensored,left,right,lower,upper,data,status)
        if (status /= 0) return
        value = native_lognormal_nnlf(data,s, loc, scale)
    end subroutine censored_lognormal_nnlf

    subroutine censored_gamma_loglikelihood(uncensored, left, right, lower, upper, &
            a, loc, scale, value, status)
        !! Evaluate gamma log likelihood using native endpoint/tail conventions.
        !! @dimension nu = size(uncensored)
        !! @dimension nl = size(left)
        !! @dimension nr = size(right)
        !! @dimension ni = size(lower)
        !! size(upper) = ni; corresponding endpoints must be paired.
        !! @require nu >= 0
        !! @require nl >= 0
        !! @require nr >= 0
        !! @require ni >= 0
        !! All group lengths are independent; empty groups are permitted.
        !! Their sum must fit the native default integer kind.
        real(dp), intent(in) :: uncensored(:) !! shape (nu); finite exact observations
        real(dp), intent(in) :: left(:) !! shape (nl); finite upper censoring limits
        real(dp), intent(in) :: right(:) !! shape (nr); finite lower censoring limits
        real(dp), intent(in) :: lower(:) !! shape (ni); interval lower endpoints
        real(dp), intent(in) :: upper(:) !! shape (ni); interval upper endpoints
        real(dp), intent(in) :: a !! finite positive gamma shape
        real(dp), intent(in) :: loc !! finite lower support endpoint
        real(dp), intent(in) :: scale !! finite positive inverse rate
        real(dp), intent(out) :: value !! log likelihood; NaN for invalid input
        integer, intent(out) :: status !! 0 evaluated; 1 count/shape; 2 data/parameters
        type(censored_data) :: data

        value = quiet_nan(0.0_dp)
        status = 2
        if (.not. valid_loc_scale(loc,scale)) return
        if (.not. ieee_is_finite(a)) return
        if (a <= 0.0_dp) return
        call pack_censored_samples(uncensored,left,right,lower,upper,data,status)
        if (status /= 0) return
        value = native_gamma_loglikelihood(data,a, loc, scale)
    end subroutine censored_gamma_loglikelihood

    subroutine censored_gamma_nnlf(uncensored, left, right, lower, upper, &
            a, loc, scale, value, status)
        !! Evaluate gamma negative log likelihood using native endpoint/tail conventions.
        !! @dimension nu = size(uncensored)
        !! @dimension nl = size(left)
        !! @dimension nr = size(right)
        !! @dimension ni = size(lower)
        !! size(upper) = ni; corresponding endpoints must be paired.
        !! @require nu >= 0
        !! @require nl >= 0
        !! @require nr >= 0
        !! @require ni >= 0
        !! All group lengths are independent; empty groups are permitted.
        !! Their sum must fit the native default integer kind.
        real(dp), intent(in) :: uncensored(:) !! shape (nu); finite exact observations
        real(dp), intent(in) :: left(:) !! shape (nl); finite upper censoring limits
        real(dp), intent(in) :: right(:) !! shape (nr); finite lower censoring limits
        real(dp), intent(in) :: lower(:) !! shape (ni); interval lower endpoints
        real(dp), intent(in) :: upper(:) !! shape (ni); interval upper endpoints
        real(dp), intent(in) :: a !! finite positive gamma shape
        real(dp), intent(in) :: loc !! finite lower support endpoint
        real(dp), intent(in) :: scale !! finite positive inverse rate
        real(dp), intent(out) :: value !! negative log likelihood; NaN for invalid input
        integer, intent(out) :: status !! 0 evaluated; 1 count/shape; 2 data/parameters
        type(censored_data) :: data

        value = quiet_nan(0.0_dp)
        status = 2
        if (.not. valid_loc_scale(loc,scale)) return
        if (.not. ieee_is_finite(a)) return
        if (a <= 0.0_dp) return
        call pack_censored_samples(uncensored,left,right,lower,upper,data,status)
        if (status /= 0) return
        value = native_gamma_nnlf(data,a, loc, scale)
    end subroutine censored_gamma_nnlf

    subroutine censored_beta_loglikelihood(uncensored, left, right, lower, upper, &
            a, b, loc, scale, value, status)
        !! Evaluate beta log likelihood using native endpoint/tail conventions.
        !! @dimension nu = size(uncensored)
        !! @dimension nl = size(left)
        !! @dimension nr = size(right)
        !! @dimension ni = size(lower)
        !! size(upper) = ni; corresponding endpoints must be paired.
        !! @require nu >= 0
        !! @require nl >= 0
        !! @require nr >= 0
        !! @require ni >= 0
        !! All group lengths are independent; empty groups are permitted.
        !! Their sum must fit the native default integer kind.
        real(dp), intent(in) :: uncensored(:) !! shape (nu); finite exact observations
        real(dp), intent(in) :: left(:) !! shape (nl); finite upper censoring limits
        real(dp), intent(in) :: right(:) !! shape (nr); finite lower censoring limits
        real(dp), intent(in) :: lower(:) !! shape (ni); interval lower endpoints
        real(dp), intent(in) :: upper(:) !! shape (ni); interval upper endpoints
        real(dp), intent(in) :: a !! finite positive first beta shape
        real(dp), intent(in) :: b !! finite positive second beta shape
        real(dp), intent(in) :: loc !! finite lower support endpoint
        real(dp), intent(in) :: scale !! finite positive support width
        real(dp), intent(out) :: value !! log likelihood; NaN for invalid input
        integer, intent(out) :: status !! 0 evaluated; 1 count/shape; 2 data/parameters
        type(censored_data) :: data

        value = quiet_nan(0.0_dp)
        status = 2
        if (.not. valid_loc_scale(loc,scale)) return
        if (.not. ieee_is_finite(a)) return
        if (a <= 0.0_dp) return
        if (.not. ieee_is_finite(b)) return
        if (b <= 0.0_dp) return
        ! Preserve the native requirement for a finite rounded upper support edge.
        if (loc > 0.0_dp) then
            if (scale > huge(1.0_dp)-loc) return
        end if
        call pack_censored_samples(uncensored,left,right,lower,upper,data,status)
        if (status /= 0) return
        value = native_beta_loglikelihood(data,a, b, loc, scale)
    end subroutine censored_beta_loglikelihood

    subroutine censored_beta_nnlf(uncensored, left, right, lower, upper, &
            a, b, loc, scale, value, status)
        !! Evaluate beta negative log likelihood using native endpoint/tail conventions.
        !! @dimension nu = size(uncensored)
        !! @dimension nl = size(left)
        !! @dimension nr = size(right)
        !! @dimension ni = size(lower)
        !! size(upper) = ni; corresponding endpoints must be paired.
        !! @require nu >= 0
        !! @require nl >= 0
        !! @require nr >= 0
        !! @require ni >= 0
        !! All group lengths are independent; empty groups are permitted.
        !! Their sum must fit the native default integer kind.
        real(dp), intent(in) :: uncensored(:) !! shape (nu); finite exact observations
        real(dp), intent(in) :: left(:) !! shape (nl); finite upper censoring limits
        real(dp), intent(in) :: right(:) !! shape (nr); finite lower censoring limits
        real(dp), intent(in) :: lower(:) !! shape (ni); interval lower endpoints
        real(dp), intent(in) :: upper(:) !! shape (ni); interval upper endpoints
        real(dp), intent(in) :: a !! finite positive first beta shape
        real(dp), intent(in) :: b !! finite positive second beta shape
        real(dp), intent(in) :: loc !! finite lower support endpoint
        real(dp), intent(in) :: scale !! finite positive support width
        real(dp), intent(out) :: value !! negative log likelihood; NaN for invalid input
        integer, intent(out) :: status !! 0 evaluated; 1 count/shape; 2 data/parameters
        type(censored_data) :: data

        value = quiet_nan(0.0_dp)
        status = 2
        if (.not. valid_loc_scale(loc,scale)) return
        if (.not. ieee_is_finite(a)) return
        if (a <= 0.0_dp) return
        if (.not. ieee_is_finite(b)) return
        if (b <= 0.0_dp) return
        ! Preserve the native requirement for a finite rounded upper support edge.
        if (loc > 0.0_dp) then
            if (scale > huge(1.0_dp)-loc) return
        end if
        call pack_censored_samples(uncensored,left,right,lower,upper,data,status)
        if (status /= 0) return
        value = native_beta_nnlf(data,a, b, loc, scale)
    end subroutine censored_beta_nnlf

    subroutine pack_censored_samples(uncensored, left, right, lower, upper, data, status)
        !! Copy validated groups into a local native object; no storage is retained.
        !! @dimension nu = size(uncensored)
        !! @dimension nl = size(left)
        !! @dimension nr = size(right)
        !! @dimension ni = size(lower)
        !! size(upper) = ni; corresponding endpoints must be paired.
        !! @require nu >= 0
        !! @require nl >= 0
        !! @require nr >= 0
        !! @require ni >= 0
        !! All group lengths are independent; empty groups are permitted.
        !! Their sum must fit the native default integer kind.
        real(dp), intent(in) :: uncensored(:) !! shape (nu); finite exact observations
        real(dp), intent(in) :: left(:) !! shape (nl); finite upper censoring limits
        real(dp), intent(in) :: right(:) !! shape (nr); finite lower censoring limits
        real(dp), intent(in) :: lower(:) !! shape (ni); interval lower endpoints
        real(dp), intent(in) :: upper(:) !! shape (ni); interval upper endpoints
        type(censored_data), intent(out) :: data !! independent native observation storage
        integer, intent(out) :: status !! 0 valid; 1 count/shape; 2 invalid data
        integer(int64) :: total
        integer :: native_status
        real(dp), allocatable :: interval(:,:)

        status = 1
        if (size(upper) /= size(lower)) return
        total = size(uncensored,kind=int64) + size(left,kind=int64) + &
            size(right,kind=int64) + size(lower,kind=int64)
        if (total > int(huge(0),int64)) return
        allocate(interval(size(lower),2))
        interval(:,1) = lower
        interval(:,2) = upper
        data = make_censored_data(uncensored,left,right,interval,native_status)
        status = 2
        if (native_status == 0) status = 0
    end subroutine pack_censored_samples
end module scifort_censored_bridge
