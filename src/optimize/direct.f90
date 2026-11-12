! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Independent implementation of the deterministic DIRECT rectangle-partitioning
! algorithm (Jones, Perttunen, and Stuckman, JOTA 79, 157-181, 1993).
! Uses potentially optimal rectangles, not a Lipschitz constant estimate.
! SciFort extension: evaluation limits are strict (SciPy's are approximate).
module scifort_optimize_direct
    use, intrinsic :: iso_fortran_env, only : int64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    implicit none
    private

    integer, parameter, public :: direct_maxfun = 1
    integer, parameter, public :: direct_maxiter = 2
    integer, parameter, public :: direct_vol_tol = 3
    integer, parameter, public :: direct_len_tol = 4
    integer, parameter, public :: direct_target = 5
    integer, parameter, public :: direct_callback_stop = 6
    integer, parameter, public :: direct_invalid = -1
    integer, parameter, public :: direct_numerical = -2

    type, public :: direct_result
        real(dp), allocatable :: x(:) !! best feasible point in user coordinates
        real(dp) :: fun = huge(1.0_dp) !! lowest finite objective encountered
        integer :: nfev = 0 !! objective evaluations (never greater than maxfun)
        integer :: nit = 0 !! completed iterations over potential optima
        integer :: nrect = 0 !! current number of hyperrectangles
        integer :: status = direct_invalid !! termination code
        logical :: success = .false. !! true only for a stopping-tolerance criterion
        character(len=160) :: message = 'invalid input' !! termination explanation
    end type direct_result

    abstract interface
        function direct_objective(x) result(f)
            import dp
            real(dp), intent(in) :: x(:) !! parameter vector of length n
            real(dp) :: f !! objective value at x
        end function direct_objective

        subroutine direct_callback(x, value, stop)
            import dp
            real(dp), intent(in) :: x(:) !! best parameters at end of iteration
            real(dp), intent(in) :: value !! smallest objective observed so far
            logical, intent(out) :: stop !! set true to terminate exploration
        end subroutine direct_callback
    end interface

    public :: direct, direct_objective, direct_callback

contains

    function direct(fun, bounds, eps, maxfun, maxiter, locally_biased, &
            f_min, f_min_rtol, vol_tol, len_tol, callback) result(res)
        procedure(direct_objective) :: fun !! deterministic scalar objective
        real(dp), intent(in) :: bounds(:,:) !! finite domain, shape (n,2), lower < upper
        real(dp), intent(in), optional :: eps !! potential-optimality threshold, default 1e-4
        integer, intent(in), optional :: maxfun !! strict evaluation budget, default 1000*n
        integer, intent(in), optional :: maxiter !! outer iteration limit, default 1000
        logical, intent(in), optional :: locally_biased !! DIRECT-L size metric, default true
        real(dp), intent(in), optional :: f_min !! known global optimum, activates target stop
        real(dp), intent(in), optional :: f_min_rtol !! target relative tolerance, default 1e-4
        real(dp), intent(in), optional :: vol_tol !! relative best-rectangle volume, default 1e-16
        real(dp), intent(in), optional :: len_tol !! relative best-rectangle size, default 1e-6
        procedure(direct_callback), optional :: callback !! optional end-of-iteration callback
        type(direct_result) :: res !! best point, evaluation counts, and status

        real(dp), allocatable :: centers(:,:), widths(:,:), values(:), measures(:)
        real(dp), allocatable :: proposed_x(:), original(:), left(:), right(:), scoring(:)
        real(dp), allocatable :: trial(:), scale(:)
        integer, allocatable :: champions(:), selected(:), axes(:), order(:)
        integer :: n, capacity, nrect, eval_limit, iter_limit, stat
        integer :: i, j, k, g, ng, ns, idx, a, b, naxes, use_axes, nselected
        integer :: newcap, jbest, tmp, iteration
        real(dp) :: f_eps, low_k, high_k, s, di, dj
        real(dp) :: maxwidth, delta, fv, loc_eps, volume_tol, length_tol
        real(dp) :: f_target, f_target_rtol
        logical :: local_bias, active_target, eligible, stopped, failed, changed

        n = size(bounds,1)
        if (n < 1 .or. size(bounds,2) /= 2) then
            res%message = 'bounds must have shape (n,2), with n >= 1'
            return
        end if
        if (.not. all(ieee_is_finite(bounds))) then
            res%message = 'bounds must be finite'
            return
        end if
        if (any(bounds(:,1) >= bounds(:,2))) then
            res%message = 'bounds require lower < upper in every dimension'
            return
        end if
        do i = 1, n
            if (bounds(i,1) < 0.0_dp .and. bounds(i,2) > 0.0_dp) then
                if (bounds(i,2) > huge(1.0_dp) + bounds(i,1)) then
                    res%message = 'interval width overflows real(dp)'
                    return
                end if
            end if
        end do
        eval_limit = int(min(1000000_int64, 1000_int64*int(n,int64)))
        if (present(maxfun)) eval_limit = maxfun
        iter_limit = 1000
        if (present(maxiter)) iter_limit = maxiter
        loc_eps = 1.0e-4_dp
        if (present(eps)) loc_eps = eps
        volume_tol = 1.0e-16_dp
        if (present(vol_tol)) volume_tol = vol_tol
        length_tol = 1.0e-6_dp
        if (present(len_tol)) length_tol = len_tol
        f_target_rtol = 1.0e-4_dp
        if (present(f_min_rtol)) f_target_rtol = f_min_rtol
        local_bias = .true.
        if (present(locally_biased)) local_bias = locally_biased
        active_target = present(f_min)
        f_target = 0.0_dp
        if (active_target) f_target = f_min
        if (eval_limit < 1 .or. iter_limit < 0) then
            res%message = 'maxfun must be positive, maxiter nonnegative'
            return
        end if
        if (.not. ieee_is_finite(loc_eps) .or. .not. ieee_is_finite(volume_tol) .or. &
                .not. ieee_is_finite(length_tol) .or. .not. ieee_is_finite(f_target_rtol)) then
            res%message = 'eps and stop tolerances must be finite'
            return
        end if
        if (loc_eps < 0.0_dp .or. volume_tol < 0.0_dp .or. volume_tol > 1.0_dp .or. &
                length_tol < 0.0_dp .or. length_tol > 1.0_dp .or. &
                f_target_rtol < 0.0_dp .or. f_target_rtol > 1.0_dp) then
            res%message = 'invalid eps, vol_tol, len_tol, or f_min_rtol'
            return
        end if
        if (active_target) then
            if (.not. ieee_is_finite(f_target)) then
                res%message = 'f_min must be finite when supplied'
                return
            end if
        end if

        allocate(scale(n), trial(n), proposed_x(n), original(n), &
            left(n), right(n), scoring(n), axes(n), order(n), stat=stat)
        if (stat /= 0) then
            res%message = 'cannot allocate optimizer work arrays'
            return
        end if
        scale = bounds(:,2) - bounds(:,1)
        capacity = min(eval_limit, 128)
        allocate(centers(n,capacity), widths(n,capacity), values(capacity), &
            measures(capacity), champions(eval_limit), selected(eval_limit), stat=stat)
        if (stat /= 0) then
            res%message = 'cannot allocate rectangle arrays'
            return
        end if
        nrect = 1
        centers(:,1) = 0.5_dp
        widths(:,1) = 1.0_dp
        proposed_x = bounds(:,1) + 0.5_dp * scale
        fv = fun(proposed_x)
        res%nfev = 1
        if (.not. ieee_is_finite(fv)) then
            res%status = direct_numerical
            res%message = 'objective returned NaN or infinity at initial point'
            return
        end if
        values(1) = fv
        measures(1) = rectangle_measure(widths(:,1), local_bias)
        res%x = proposed_x
        res%fun = fv
        res%nrect = 1

        do iteration = 1, iter_limit
            jbest = best_rectangle()
            if (active_target) then
                if (abs(res%fun - f_target) <= f_target_rtol * max(1.0_dp, abs(f_target))) then
                    call finish(direct_target, 'known optimum reached', .true.)
                    return
                end if
            end if
            if (product(widths(:,jbest)) <= volume_tol) then
                call finish(direct_vol_tol, 'best rectangle volume tolerance reached', .true.)
                return
            end if
            if (rectangle_measure(widths(:,jbest), local_bias) <= length_tol) then
                call finish(direct_len_tol, 'best rectangle length tolerance reached', .true.)
                return
            end if
            if (eval_limit - res%nfev < 2) then
                call finish(direct_maxfun, 'evaluation budget reached', .false.)
                return
            end if

            ! Retain only the lowest-value rectangle at each distinct size.
            ! These candidates form the points used to identify the lower hull.
            ng = 0
            do i = 1, nrect
                g = 0
                do j = 1, ng
                    idx = champions(j)
                    if (same_size(measures(i), measures(idx))) then
                        g = j
                        exit
                    end if
                end do
                if (g == 0) then
                    ng = ng + 1
                    champions(ng) = i
                else if (values(i) < values(champions(g))) then
                    champions(g) = i
                end if
            end do
            f_eps = res%fun - loc_eps * abs(res%fun)
            ns = 0
            do g = 1, ng
                idx = champions(g)
                dj = measures(idx)
                low_k = max(0.0_dp, (values(idx) - f_eps) / dj)
                high_k = huge(1.0_dp)
                eligible = .true.
                do j = 1, ng
                    if (j == g) cycle
                    i = champions(j)
                    di = measures(i)
                    if (same_size(di,dj)) cycle
                    s = (values(idx) - values(i)) / (dj - di)
                    if (dj > di) then
                        low_k = max(low_k, s)
                    else
                        high_k = min(high_k, s)
                    end if
                end do
                eligible = (low_k <= high_k)
                if (.not. eligible) then
                    eligible = low_k - 64.0_dp*epsilon(1.0_dp)* &
                        max(abs(low_k),abs(high_k)) <= high_k
                end if
                if (eligible) then
                    ns = ns + 1
                    selected(ns) = idx
                end if
            end do
            if (ns == 0) then
                ! With nearly equal objective values, rounding may empty the
                ! hull. Divide a largest rectangle instead of stalling.
                ns = 1
                idx = maxloc(measures(1:nrect),dim=1)
                selected(1) = idx
            end if
            ! Split larger selected rectangles first if close to maxfun.
            do i = 2, ns
                tmp = selected(i)
                j = i - 1
                do while (j >= 1)
                    if (measures(selected(j)) >= measures(tmp)) exit
                    selected(j+1) = selected(j)
                    j = j - 1
                end do
                selected(j+1) = tmp
            end do

            changed = .false.
            do nselected = 1, ns
                idx = selected(nselected)
                original = widths(:,idx)
                maxwidth = maxval(original)
                naxes = 0
                do i = 1, n
                    if (original(i) >= maxwidth * (1.0_dp - 32.0_dp * epsilon(1.0_dp))) then
                        naxes = naxes + 1
                        axes(naxes) = i
                    end if
                end do
                use_axes = min(naxes, (eval_limit-res%nfev)/2)
                if (use_axes == 0) exit
                ! Two off-center values per longest axis, before splitting.
                do k = 1, use_axes
                    a = axes(k)
                    delta = original(a) / 3.0_dp
                    trial = centers(:,idx)
                    trial(a) = trial(a) - delta
                    call evaluate(trial, left(k), failed)
                    if (failed) then
                        call finish(direct_numerical, 'objective returned nonfinite value', .false.)
                        return
                    end if
                    trial(a) = trial(a) + 2.0_dp*delta
                    call evaluate(trial, right(k), failed)
                    if (failed) then
                        call finish(direct_numerical, 'objective returned nonfinite value', .false.)
                        return
                    end if
                    scoring(k) = min(left(k), right(k))
                    order(k) = k
                end do
                ! Split dimensions with best off-center samples first.
                do i = 2, use_axes
                    tmp = order(i)
                    j = i - 1
                    do while (j >= 1)
                        if (scoring(order(j)) <= scoring(tmp)) exit
                        order(j+1) = order(j)
                        j = j-1
                    end do
                    order(j+1) = tmp
                end do
                call grow(nrect+2*use_axes, failed)
                if (failed) then
                    call finish(direct_numerical, 'cannot grow rectangle storage', .false.)
                    return
                end if
                do k = 1, use_axes
                    b = order(k)
                    a = axes(b)
                    delta = original(a) / 3.0_dp
                    widths(a,idx) = widths(a,idx) / 3.0_dp
                    centers(:,nrect+1) = centers(:,idx)
                    centers(:,nrect+2) = centers(:,idx)
                    centers(a,nrect+1) = centers(a,nrect+1) - delta
                    centers(a,nrect+2) = centers(a,nrect+2) + delta
                    widths(:,nrect+1) = widths(:,idx)
                    widths(:,nrect+2) = widths(:,idx)
                    values(nrect+1) = left(b)
                    values(nrect+2) = right(b)
                    measures(nrect+1) = rectangle_measure(widths(:,idx), local_bias)
                    measures(nrect+2) = measures(nrect+1)
                    nrect = nrect + 2
                end do
                measures(idx) = rectangle_measure(widths(:,idx), local_bias)
                changed = .true.
            end do
            if (.not. changed) then
                call finish(direct_maxfun, 'evaluation budget reached', .false.)
                return
            end if
            res%nit = iteration
            res%nrect = nrect
            if (present(callback)) then
                stopped = .false.
                call callback(res%x, res%fun, stopped)
                if (stopped) then
                    call finish(direct_callback_stop, 'stopped by callback', .false.)
                    return
                end if
            end if
        end do
        if (iter_limit == 0) then
            call finish(direct_maxiter, 'iteration budget reached', .false.)
        else
            ! A final test is useful when the last split met a tolerance.
            jbest = best_rectangle()
            if (active_target) then
                if (abs(res%fun - f_target) <= f_target_rtol * max(1.0_dp,abs(f_target))) then
                    call finish(direct_target, 'known optimum reached', .true.)
                    return
                end if
            end if
            if (product(widths(:,jbest)) <= volume_tol) then
                call finish(direct_vol_tol, 'best rectangle volume tolerance reached', .true.)
            else if (rectangle_measure(widths(:,jbest), local_bias) <= length_tol) then
                call finish(direct_len_tol, 'best rectangle length tolerance reached', .true.)
            else
                call finish(direct_maxiter, 'iteration budget reached', .false.)
            end if
        end if

    contains

        subroutine evaluate(u, value, invalid)
            real(dp), intent(in) :: u(:) !! normalized sample point
            real(dp), intent(out) :: value !! objective at sample point
            logical, intent(out) :: invalid !! true for nonfinite evaluations
            proposed_x = bounds(:,1) + u*scale
            value = fun(proposed_x)
            res%nfev = res%nfev + 1
            invalid = .not. ieee_is_finite(value)
            if (invalid) return
            if (value < res%fun) then
                res%fun = value
                res%x = proposed_x
            end if
        end subroutine evaluate

        subroutine grow(required, failed_alloc)
            integer, intent(in) :: required !! minimum rectangle capacity
            logical, intent(out) :: failed_alloc !! allocation failed
            real(dp), allocatable :: cc(:,:), ww(:,:), vv(:), dd(:)
            integer :: st
            failed_alloc = .false.
            if (required <= capacity) return
            newcap = min(eval_limit, max(required, min(eval_limit,2*capacity)))
            allocate(cc(n,newcap), ww(n,newcap), vv(newcap), dd(newcap), stat=st)
            if (st /= 0) then
                failed_alloc = .true.
                return
            end if
            cc(:,1:nrect) = centers(:,1:nrect)
            ww(:,1:nrect) = widths(:,1:nrect)
            vv(1:nrect) = values(1:nrect)
            dd(1:nrect) = measures(1:nrect)
            call move_alloc(cc,centers)
            call move_alloc(ww,widths)
            call move_alloc(vv,values)
            call move_alloc(dd,measures)
            capacity = newcap
        end subroutine grow

        integer function best_rectangle() result(ibest)
            integer :: p
            real(dp) :: smallest_value
            ibest = 1
            smallest_value = values(1)
            do p = 2, nrect
                if (values(p) < smallest_value) then
                    ibest = p
                    smallest_value = values(p)
                else if (values(p) <= smallest_value) then
                    if (measures(p) < measures(ibest)) ibest = p
                end if
            end do
        end function best_rectangle

        subroutine finish(code, text_message, converged)
            integer, intent(in) :: code !! result termination status
            character(len=*), intent(in) :: text_message !! status description
            logical, intent(in) :: converged !! whether tolerances were satisfied
            res%status = code
            res%message = text_message
            res%success = converged
            res%nrect = nrect
        end subroutine finish

    end function direct

    pure real(dp) function rectangle_measure(width, local_bias) result(d)
        real(dp), intent(in) :: width(:) !! normalized edge lengths
        logical, intent(in) :: local_bias !! use maximum edge length for DIRECT-L
        if (local_bias) then
            d = 0.5_dp * maxval(width)
        else
            d = 0.5_dp * sqrt(sum(width**2))
        end if
    end function rectangle_measure

    pure logical function same_size(a,b) result(equal)
        real(dp), intent(in) :: a !! first rectangle size
        real(dp), intent(in) :: b !! second rectangle size
        equal = abs(a-b) <= 32.0_dp*epsilon(1.0_dp)*max(a,b)
    end function same_size

end module scifort_optimize_direct
