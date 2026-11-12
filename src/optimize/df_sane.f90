! SPDX-License-Identifier: BSD-3-Clause
! Copyright (c) 2026 SciFort contributors
! Copyright (c) 2001-2002 Enthought, Inc. 2003, SciPy Developers.
!
! Derivative-free spectral residual solver for the scipy.optimize translation.
!
! The DF-SANE iteration, spectral-coefficient safeguards, convergence rules,
! and Cruz/Cheng nonmonotone line searches are adapted from SciPy 1.17.0
! scipy/optimize/_spectral.py and scipy/optimize/_linesearch.py under SciPy's
! BSD-3-Clause license.

module scifort_optimize_df_sane
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    implicit none
    private

    abstract interface
        subroutine df_sane_function(x, fvec)
            import dp
            implicit none
            real(dp), intent(in) :: x(:) !! independent variable vector
            real(dp), intent(out) :: fvec(:) !! residual vector at x
        end subroutine df_sane_function
    end interface

    public :: df_sane_function
    public :: solve_df_sane

contains

    subroutine solve_df_sane(fun, x0, x, fvec, nfev, nit, status, success, &
            message, tol, ftol, fatol, maxfev, history, sigma_eps, sigma_0, &
            line_search)
        procedure(df_sane_function) :: fun !! residual callback
        real(dp), intent(in) :: x0(:) !! starting vector
        real(dp), intent(out) :: x(:) !! final solution estimate
        real(dp), intent(out) :: fvec(:) !! residual vector at x
        integer, intent(out) :: nfev !! residual evaluations
        integer, intent(out) :: nit !! accepted spectral iterations
        integer, intent(out) :: status !! 1 success, 2 evaluation limit, 0 invalid input
        logical, intent(out) :: success !! true when status is 1
        character(len=*), intent(out) :: message !! termination description
        real(dp), intent(in), optional :: tol !! generic root relative norm tolerance
        real(dp), intent(in), optional :: ftol !! relative residual norm tolerance
        real(dp), intent(in), optional :: fatol !! absolute residual norm tolerance
        integer, intent(in), optional :: maxfev !! maximum residual evaluations
        integer, intent(in), optional :: history !! nonmonotone window length, default 10
        real(dp), intent(in), optional :: sigma_eps !! spectral safeguard, default 1e-10
        real(dp), intent(in), optional :: sigma_0 !! initial spectral coefficient
        character(len=*), intent(in), optional :: line_search !! cruz or cheng

        real(dp), allocatable :: d(:)
        real(dp), allocatable :: fnew(:)
        real(dp), allocatable :: prev_fs(:)
        real(dp), allocatable :: s(:)
        real(dp), allocatable :: xnew(:)
        real(dp), allocatable :: y(:)
        real(dp) :: c_control
        real(dp) :: denom
        real(dp) :: eta
        real(dp) :: f0_norm
        real(dp) :: fatol_use
        real(dp) :: fk
        real(dp) :: fnorm
        real(dp) :: ftol_use
        real(dp) :: q_control
        real(dp) :: sigma
        real(dp) :: sigma_eps_use
        real(dp) :: sigma_limit
        integer :: history_count
        integer :: history_use
        integer :: maxfev_use
        integer :: n
        logical :: accepted
        character(len=16) :: search

        n = size(x0)
        x = 0.0_dp
        fvec = 0.0_dp
        nfev = 0
        nit = 0
        status = 0
        success = .false.
        message = 'improper input parameters'

        if (n <= 0 .or. size(x) /= n .or. size(fvec) /= n) return
        if (.not. all(ieee_is_finite(x0))) return

        ftol_use = 1.0e-8_dp
        if (present(tol)) ftol_use = tol
        if (present(ftol)) ftol_use = ftol
        fatol_use = 1.0e-300_dp
        if (present(fatol)) fatol_use = fatol
        if (.not. ieee_is_finite(ftol_use) .or. ftol_use < 0.0_dp) return
        if (.not. ieee_is_finite(fatol_use) .or. fatol_use < 0.0_dp) return

        maxfev_use = 1000
        if (present(maxfev)) maxfev_use = maxfev
        if (maxfev_use < 1) return

        history_use = 10
        if (present(history)) history_use = history
        if (history_use < 1) return

        sigma_eps_use = 1.0e-10_dp
        if (present(sigma_eps)) sigma_eps_use = sigma_eps
        if (.not. ieee_is_finite(sigma_eps_use)) return
        if (sigma_eps_use <= 0.0_dp .or. sigma_eps_use >= 1.0_dp) return
        sigma_limit = 1.0_dp / sigma_eps_use

        sigma = 1.0_dp
        if (present(sigma_0)) sigma = sigma_0
        if (.not. ieee_is_finite(sigma)) return

        search = 'cruz'
        if (present(line_search)) search = lowercase(trim(line_search))
        select case (trim(search))
        case ('cruz', 'cheng')
        case default
            message = 'unsupported line search; use cruz or cheng'
            return
        end select

        allocate(d(n), fnew(n), prev_fs(history_use), s(n), xnew(n), y(n))
        x = x0
        call fun(x, fvec)
        nfev = 1
        if (.not. all(ieee_is_finite(fvec))) then
            message = 'nonfinite residual encountered'
            return
        end if

        fk = dot_product(fvec, fvec)
        f0_norm = sqrt(fk)
        history_count = 1
        prev_fs = 0.0_dp
        prev_fs(1) = fk
        q_control = 1.0_dp
        c_control = fk

        do
            fnorm = sqrt(fk)
            if (fnorm < ftol_use * f0_norm + fatol_use) then
                status = 1
                success = .true.
                message = 'successful convergence'
                exit
            end if

            if (abs(sigma) > sigma_limit) then
                sigma = sign(sigma_limit, sigma)
            else if (abs(sigma) < sigma_eps_use) then
                sigma = sigma_eps_use
            end if

            d = -sigma * fvec
            ! SciPy's default eta uses the initial squared residual, not the current one.
            eta = f0_norm * f0_norm / real(nit + 1, dp)**2

            select case (trim(search))
            case ('cruz')
                call cruz_search(fun, x, d, fk, prev_fs, history_count, maxfev_use, &
                    nfev, eta, xnew, fnew, accepted)
            case ('cheng')
                call cheng_search(fun, x, d, fk, c_control, q_control, maxfev_use, &
                    nfev, eta, xnew, fnew, accepted)
            end select

            if (.not. accepted) then
                status = 2
                message = 'too many function evaluations required'
                exit
            end if

            s = xnew - x
            y = fnew - fvec
            denom = dot_product(s, y)
            if (abs(denom) > tiny(1.0_dp)) then
                sigma = dot_product(s, s) / denom
            else if (denom < 0.0_dp) then
                sigma = -sigma_limit
            else
                sigma = sigma_limit
            end if

            x = xnew
            fvec = fnew
            fk = dot_product(fvec, fvec)
            if (.not. ieee_is_finite(fk)) then
                status = 0
                message = 'nonfinite residual encountered'
                exit
            end if

            if (trim(search) == 'cruz') then
                call append_history(prev_fs, history_count, fk)
            end if
            nit = nit + 1
        end do
    end subroutine solve_df_sane

    subroutine cruz_search(fun, x, d, fk, prev_fs, history_count, maxfev, nfev, &
            eta, xnew, fnew, accepted)
        procedure(df_sane_function) :: fun !! residual callback
        real(dp), intent(in) :: x(:) !! current iterate
        real(dp), intent(in) :: d(:) !! spectral search direction
        real(dp), intent(in) :: fk !! current squared residual norm
        real(dp), intent(in) :: prev_fs(:) !! recent squared residual norms
        integer, intent(in) :: history_count !! number of valid history values
        integer, intent(in) :: maxfev !! evaluation cap
        integer, intent(inout) :: nfev !! residual evaluation counter
        real(dp), intent(in) :: eta !! nonmonotone allowance
        real(dp), intent(out) :: xnew(:) !! accepted trial point
        real(dp), intent(out) :: fnew(:) !! residual at accepted point
        logical, intent(out) :: accepted !! true when a trial is accepted

        real(dp) :: alpha_m
        real(dp) :: alpha_p
        real(dp) :: alpha_tm
        real(dp) :: alpha_tp
        real(dp) :: fbar
        real(dp) :: fp
        real(dp), parameter :: gamma = 1.0e-4_dp
        real(dp), parameter :: tau_max = 0.5_dp
        real(dp), parameter :: tau_min = 0.1_dp

        accepted = .false.
        xnew = x
        fnew = 0.0_dp
        fbar = maxval(prev_fs(1:history_count))
        alpha_p = 1.0_dp
        alpha_m = 1.0_dp

        do
            if (.not. evaluate_trial(fun, x + alpha_p * d, maxfev, nfev, &
                    xnew, fnew, fp)) return
            if (fp <= fbar + eta - gamma * alpha_p**2 * fk) then
                accepted = .true.
                return
            end if
            alpha_tp = trial_scale(alpha_p, fk, fp)

            if (.not. evaluate_trial(fun, x - alpha_m * d, maxfev, nfev, &
                    xnew, fnew, fp)) return
            if (fp <= fbar + eta - gamma * alpha_m**2 * fk) then
                accepted = .true.
                return
            end if
            alpha_tm = trial_scale(alpha_m, fk, fp)

            alpha_p = clip(alpha_tp, tau_min * alpha_p, tau_max * alpha_p)
            alpha_m = clip(alpha_tm, tau_min * alpha_m, tau_max * alpha_m)
        end do
    end subroutine cruz_search

    subroutine cheng_search(fun, x, d, fk, c_control, q_control, maxfev, nfev, &
            eta, xnew, fnew, accepted)
        procedure(df_sane_function) :: fun !! residual callback
        real(dp), intent(in) :: x(:) !! current iterate
        real(dp), intent(in) :: d(:) !! spectral search direction
        real(dp), intent(in) :: fk !! current squared residual norm
        real(dp), intent(inout) :: c_control !! nonmonotone control value
        real(dp), intent(inout) :: q_control !! nonmonotone control weight
        integer, intent(in) :: maxfev !! evaluation cap
        integer, intent(inout) :: nfev !! residual evaluation counter
        real(dp), intent(in) :: eta !! nonmonotone allowance
        real(dp), intent(out) :: xnew(:) !! accepted trial point
        real(dp), intent(out) :: fnew(:) !! residual at accepted point
        logical, intent(out) :: accepted !! true when a trial is accepted

        real(dp) :: alpha_m
        real(dp) :: alpha_p
        real(dp) :: alpha_tm
        real(dp) :: alpha_tp
        real(dp) :: fp
        real(dp) :: q_next
        real(dp), parameter :: gamma = 1.0e-4_dp
        real(dp), parameter :: nu = 0.85_dp
        real(dp), parameter :: tau_max = 0.5_dp
        real(dp), parameter :: tau_min = 0.1_dp

        accepted = .false.
        xnew = x
        fnew = 0.0_dp
        alpha_p = 1.0_dp
        alpha_m = 1.0_dp

        do
            if (.not. evaluate_trial(fun, x + alpha_p * d, maxfev, nfev, &
                    xnew, fnew, fp)) return
            if (fp <= c_control + eta - gamma * alpha_p**2 * fk) exit
            alpha_tp = trial_scale(alpha_p, fk, fp)

            if (.not. evaluate_trial(fun, x - alpha_m * d, maxfev, nfev, &
                    xnew, fnew, fp)) return
            if (fp <= c_control + eta - gamma * alpha_m**2 * fk) exit
            alpha_tm = trial_scale(alpha_m, fk, fp)

            alpha_p = clip(alpha_tp, tau_min * alpha_p, tau_max * alpha_p)
            alpha_m = clip(alpha_tm, tau_min * alpha_m, tau_max * alpha_m)
        end do

        q_next = nu * q_control + 1.0_dp
        c_control = (nu * q_control * (c_control + eta) + fp) / q_next
        q_control = q_next
        accepted = .true.
    end subroutine cheng_search

    logical function evaluate_trial(fun, trial, maxfev, nfev, xnew, fnew, fp) &
            result(ok)
        procedure(df_sane_function) :: fun !! residual callback
        real(dp), intent(in) :: trial(:) !! trial point
        integer, intent(in) :: maxfev !! evaluation cap
        integer, intent(inout) :: nfev !! residual evaluation counter
        real(dp), intent(out) :: xnew(:) !! copied trial point
        real(dp), intent(out) :: fnew(:) !! trial residual
        real(dp), intent(out) :: fp !! squared trial residual norm

        if (nfev >= maxfev) then
            ok = .false.
            return
        end if
        xnew = trial
        call fun(xnew, fnew)
        nfev = nfev + 1
        fp = dot_product(fnew, fnew)
        ok = all(ieee_is_finite(fnew)) .and. ieee_is_finite(fp)
    end function evaluate_trial

    pure function trial_scale(alpha, fk, fp) result(value)
        real(dp), intent(in) :: alpha !! current positive/negative trial magnitude
        real(dp), intent(in) :: fk !! current squared residual norm
        real(dp), intent(in) :: fp !! trial squared residual norm
        real(dp) :: denom
        real(dp) :: value

        denom = fp + (2.0_dp * alpha - 1.0_dp) * fk
        if (abs(denom) <= tiny(1.0_dp)) then
            value = 0.5_dp * alpha
        else
            value = alpha**2 * fk / denom
        end if
    end function trial_scale

    subroutine append_history(history, count, value)
        real(dp), intent(inout) :: history(:) !! rolling history storage
        integer, intent(inout) :: count !! valid history entries
        real(dp), intent(in) :: value !! newest merit value

        if (count < size(history)) then
            count = count + 1
            history(count) = value
        else if (size(history) > 1) then
            history(1:size(history)-1) = history(2:size(history))
            history(size(history)) = value
        else
            history(1) = value
        end if
    end subroutine append_history

    pure function clip(value, lower, upper) result(clipped)
        real(dp), intent(in) :: value !! value to clip
        real(dp), intent(in) :: lower !! lower bound
        real(dp), intent(in) :: upper !! upper bound
        real(dp) :: clipped

        clipped = max(lower, min(upper, value))
    end function clip

    pure function lowercase(text) result(lower)
        character(len=*), intent(in) :: text !! source text
        character(len=len(text)) :: lower
        integer :: code
        integer :: i

        lower = text
        do i = 1, len(text)
            code = iachar(text(i:i))
            if (code >= iachar('A') .and. code <= iachar('Z')) then
                lower(i:i) = achar(code + iachar('a') - iachar('A'))
            end if
        end do
    end function lowercase

end module scifort_optimize_df_sane
