! SPDX-License-Identifier: BSD-3-Clause
! Copyright (c) 2026 SciFort contributors
! Copyright (c) 2001-2002 Enthought, Inc. 2003, SciPy Developers.
!
! Strong-Wolfe scalar line search adapted from SciPy 1.17.0
! scipy/optimize/_linesearch.py and scipy/optimize/_dcsrch.py.  The latter
! is SciPy's translation of the MINPACK-2 DCSRCH/DCSTEP algorithm by
! More, Thuente, Averick, and Carter.

module scifort_optimize_wolfe
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    implicit none
    private

    abstract interface
        function scalar_line_function(step) result(value)
            import dp
            implicit none
            real(dp), intent(in) :: step !! scalar trial step
            real(dp) :: value !! scalar objective or directional derivative
        end function scalar_line_function

        logical function scalar_line_condition(step, value) result(accept)
            import dp
            implicit none
            real(dp), intent(in) :: step !! strong-Wolfe candidate step
            real(dp), intent(in) :: value !! objective at the candidate step
        end function scalar_line_condition
    end interface

    public :: scalar_line_function
    public :: scalar_line_condition
    public :: scalar_search_wolfe1
    public :: scalar_search_wolfe2

contains

    subroutine scalar_search_wolfe1(phi, derphi, phi0, step, phi1, success, &
            derphi0, c1, c2, amax, amin, xtol, old_phi0)
        procedure(scalar_line_function) :: phi !! scalar objective along search line
        procedure(scalar_line_function) :: derphi !! directional derivative callback
        real(dp), intent(in) :: phi0 !! objective at step zero
        real(dp), intent(out) :: step !! accepted step, undefined when success is false
        real(dp), intent(out) :: phi1 !! objective at accepted/current step
        logical, intent(out) :: success !! true when strong-Wolfe conditions are met
        real(dp), intent(in), optional :: derphi0 !! derivative at step zero
        real(dp), intent(in), optional :: c1 !! sufficient-decrease parameter
        real(dp), intent(in), optional :: c2 !! curvature parameter
        real(dp), intent(in), optional :: amax !! upper step bound
        real(dp), intent(in), optional :: amin !! lower step bound
        real(dp), intent(in), optional :: xtol !! relative interval tolerance
        real(dp), intent(in), optional :: old_phi0 !! objective value at previous outer iterate

        real(dp) :: brack_hi
        real(dp) :: brack_lo
        real(dp) :: c1_use
        real(dp) :: c2_use
        real(dp) :: finit
        real(dp) :: ftest
        real(dp) :: g
        real(dp) :: ginit
        real(dp) :: gtest
        real(dp) :: gx
        real(dp) :: gy
        real(dp) :: phi_current
        real(dp) :: stmax
        real(dp) :: stmin
        real(dp) :: stpmax
        real(dp) :: stpmin
        real(dp) :: width
        real(dp) :: width1
        real(dp) :: xtol_use
        real(dp) :: fx
        real(dp) :: fy
        integer :: iter
        integer :: stage
        logical :: brackt
        logical :: request_fg
        logical :: terminated
        logical :: warning

        c1_use = 1.0e-4_dp
        if (present(c1)) c1_use = c1
        c2_use = 0.9_dp
        if (present(c2)) c2_use = c2
        stpmax = 50.0_dp
        if (present(amax)) stpmax = amax
        stpmin = 1.0e-8_dp
        if (present(amin)) stpmin = amin
        xtol_use = 1.0e-14_dp
        if (present(xtol)) xtol_use = xtol

        step = 1.0_dp
        if (present(old_phi0) .and. present(derphi0)) then
            if (.not. same_real(derphi0, 0.0_dp)) then
                step = min(1.0_dp, 1.01_dp * 2.0_dp * (phi0 - old_phi0) / derphi0)
                if (step < 0.0_dp) step = 1.0_dp
            end if
        end if
        phi1 = phi0
        success = .false.
        if (.not. ieee_is_finite(phi0)) return
        if (present(derphi0)) then
            ginit = derphi0
        else
            ginit = derphi(0.0_dp)
        end if
        if (.not. ieee_is_finite(ginit)) return
        if (step < stpmin .or. step > stpmax .or. ginit >= 0.0_dp) return
        if (c1_use < 0.0_dp .or. c2_use < 0.0_dp .or. c1_use >= c2_use .or. &
                c2_use >= 1.0_dp .or. xtol_use < 0.0_dp .or. stpmin < 0.0_dp .or. &
                stpmax < stpmin) return

        brackt = .false.
        stage = 1
        finit = phi0
        gtest = c1_use * ginit
        width = stpmax - stpmin
        width1 = width / 0.5_dp
        brack_lo = 0.0_dp
        fx = finit
        gx = ginit
        brack_hi = 0.0_dp
        fy = finit
        gy = ginit
        stmin = 0.0_dp
        stmax = step + 4.0_dp * step
        request_fg = .true.
        terminated = .false.
        warning = .false.
        phi_current = phi0
        g = ginit

        do iter = 1, 100
            if (request_fg) then
                phi_current = phi(step)
                g = derphi(step)
                if (.not. ieee_is_finite(step) .or. .not. ieee_is_finite(phi_current) .or. &
                        .not. ieee_is_finite(g)) then
                    warning = .true.
                    exit
                end if
                request_fg = .false.
            end if

            ftest = finit + step * gtest
            if (stage == 1 .and. phi_current <= ftest .and. g >= 0.0_dp) stage = 2

            warning = .false.
            if (brackt .and. (step <= stmin .or. step >= stmax)) warning = .true.
            if (brackt .and. stmax - stmin <= xtol_use * stmax) warning = .true.
            if (same_real(step, stpmax) .and. phi_current <= ftest .and. g <= gtest) &
                warning = .true.
            if (same_real(step, stpmin) .and. (phi_current > ftest .or. g >= gtest)) &
                warning = .true.

            if (phi_current <= ftest .and. abs(g) <= c2_use * (-ginit)) then
                success = .true.
                terminated = .true.
                exit
            end if
            if (warning) then
                terminated = .true.
                exit
            end if

            if (stage == 1 .and. phi_current <= fx .and. phi_current > ftest) then
                call dcstep(brack_lo, fx, gx, brack_hi, fy, gy, step, &
                    phi_current, g, brackt, stmin, stmax, gtest)
            else
                call dcstep(brack_lo, fx, gx, brack_hi, fy, gy, step, &
                    phi_current, g, brackt, stmin, stmax)
            end if

            if (.not. ieee_is_finite(step)) then
                warning = .true.
                exit
            end if

            if (brackt) then
                if (abs(brack_hi - brack_lo) >= 0.66_dp * width1) then
                    step = brack_lo + 0.5_dp * (brack_hi - brack_lo)
                end if
                width1 = width
                width = abs(brack_hi - brack_lo)
            end if

            if (brackt) then
                stmin = min(brack_lo, brack_hi)
                stmax = max(brack_lo, brack_hi)
            else
                stmin = step + 1.1_dp * (step - brack_lo)
                stmax = step + 4.0_dp * (step - brack_lo)
            end if

            step = min(stpmax, max(stpmin, step))
            if (brackt .and. ((step <= stmin .or. step >= stmax) .or. &
                    stmax - stmin <= xtol_use * stmax)) step = brack_lo
            request_fg = .true.
        end do

        if (.not. terminated .and. .not. success) warning = .true.
        phi1 = phi_current
        if (warning) success = .false.
    end subroutine scalar_search_wolfe1


    subroutine scalar_search_wolfe2(phi, derphi, phi0, step, phi1, success, &
            derphi0, c1, c2, amax, old_phi0, extra_condition, maxiter)
        procedure(scalar_line_function) :: phi !! scalar objective along search line
        procedure(scalar_line_function) :: derphi !! directional derivative callback
        real(dp), intent(in) :: phi0 !! objective at step zero
        real(dp), intent(out) :: step !! accepted step, undefined when success is false
        real(dp), intent(out) :: phi1 !! objective at accepted/current step
        logical, intent(out) :: success !! true when strong-Wolfe conditions are met
        real(dp), intent(in), optional :: derphi0 !! derivative at step zero
        real(dp), intent(in), optional :: c1 !! sufficient-decrease parameter
        real(dp), intent(in), optional :: c2 !! curvature parameter
        real(dp), intent(in), optional :: amax !! upper step bound
        real(dp), intent(in), optional :: old_phi0 !! objective at previous outer iterate
        procedure(scalar_line_condition), optional :: extra_condition !! extra acceptance rule
        integer, intent(in), optional :: maxiter !! outer bracketing iterations, default 10

        real(dp) :: alpha0
        real(dp) :: alpha1
        real(dp) :: alpha2
        real(dp) :: amax_use
        real(dp) :: c1_use
        real(dp) :: c2_use
        real(dp) :: derphi_a0
        real(dp) :: derphi_a1
        real(dp) :: derphi0_use
        real(dp) :: phi_a0
        real(dp) :: phi_a1
        real(dp) :: phi_star
        integer :: i
        integer :: maxiter_use
        logical :: accepted
        logical :: zoom_ok

        c1_use = 1.0e-4_dp
        if (present(c1)) c1_use = c1
        c2_use = 0.9_dp
        if (present(c2)) c2_use = c2
        amax_use = huge(1.0_dp)
        if (present(amax)) amax_use = amax
        maxiter_use = 10
        if (present(maxiter)) maxiter_use = maxiter

        step = 0.0_dp
        phi1 = phi0
        success = .false.
        if (present(derphi0)) then
            derphi0_use = derphi0
        else
            derphi0_use = derphi(0.0_dp)
        end if
        if (.not. ieee_is_finite(phi0) .or. .not. ieee_is_finite(derphi0_use)) return
        if (c1_use <= 0.0_dp .or. c1_use >= c2_use .or. c2_use >= 1.0_dp) return
        if (amax_use <= 0.0_dp .or. maxiter_use < 0) return

        alpha0 = 0.0_dp
        alpha1 = 1.0_dp
        if (present(old_phi0) .and. .not. same_real(derphi0_use, 0.0_dp)) then
            alpha1 = min(1.0_dp, 1.01_dp * 2.0_dp * (phi0 - old_phi0) / derphi0_use)
        end if
        if (alpha1 < 0.0_dp) alpha1 = 1.0_dp
        alpha1 = min(alpha1, amax_use)
        phi_a1 = phi(alpha1)
        phi_a0 = phi0
        derphi_a0 = derphi0_use

        do i = 0, maxiter_use - 1
            if (same_real(alpha1, 0.0_dp) .or. alpha0 > amax_use) return
            if (phi_a1 > phi0 + c1_use * alpha1 * derphi0_use .or. &
                    (i > 0 .and. phi_a1 >= phi_a0)) then
                call wolfe_zoom(alpha0, alpha1, phi_a0, phi_a1, derphi_a0, phi, derphi, &
                    phi0, derphi0_use, c1_use, c2_use, step, phi_star, zoom_ok, &
                    extra_condition)
                if (zoom_ok) then
                    phi1 = phi_star
                    success = .true.
                end if
                return
            end if

            derphi_a1 = derphi(alpha1)
            if (abs(derphi_a1) <= -c2_use * derphi0_use) then
                accepted = .true.
                if (present(extra_condition)) accepted = extra_condition(alpha1, phi_a1)
                if (accepted) then
                    step = alpha1
                    phi1 = phi_a1
                    success = .true.
                    return
                end if
            end if

            if (derphi_a1 >= 0.0_dp) then
                call wolfe_zoom(alpha1, alpha0, phi_a1, phi_a0, derphi_a1, phi, derphi, &
                    phi0, derphi0_use, c1_use, c2_use, step, phi_star, zoom_ok, &
                    extra_condition)
                if (zoom_ok) then
                    phi1 = phi_star
                    success = .true.
                end if
                return
            end if

            alpha2 = min(2.0_dp * alpha1, amax_use)
            alpha0 = alpha1
            alpha1 = alpha2
            phi_a0 = phi_a1
            phi_a1 = phi(alpha1)
            derphi_a0 = derphi_a1
        end do

        step = alpha1
        phi1 = phi_a1
    end subroutine scalar_search_wolfe2

    subroutine wolfe_zoom(a_lo_in, a_hi_in, phi_lo_in, phi_hi_in, derphi_lo_in, &
            phi, derphi, phi0, derphi0, c1, c2, step, phi_star, success, &
            extra_condition)
        real(dp), intent(in) :: a_lo_in !! initial low bracket endpoint
        real(dp), intent(in) :: a_hi_in !! initial high bracket endpoint
        real(dp), intent(in) :: phi_lo_in !! objective at low endpoint
        real(dp), intent(in) :: phi_hi_in !! objective at high endpoint
        real(dp), intent(in) :: derphi_lo_in !! derivative at low endpoint
        procedure(scalar_line_function) :: phi !! scalar objective along search line
        procedure(scalar_line_function) :: derphi !! directional derivative callback
        real(dp), intent(in) :: phi0 !! objective at step zero
        real(dp), intent(in) :: derphi0 !! derivative at step zero
        real(dp), intent(in) :: c1 !! sufficient-decrease parameter
        real(dp), intent(in) :: c2 !! curvature parameter
        real(dp), intent(out) :: step !! accepted step
        real(dp), intent(out) :: phi_star !! objective at accepted step
        logical, intent(out) :: success !! true when zoom found an acceptable step
        procedure(scalar_line_condition), optional :: extra_condition !! extra acceptance rule

        real(dp), parameter :: delta1 = 0.2_dp
        real(dp), parameter :: delta2 = 0.1_dp
        real(dp) :: a
        real(dp) :: a_hi
        real(dp) :: a_j
        real(dp) :: a_lo
        real(dp) :: a_rec
        real(dp) :: b
        real(dp) :: cchk
        real(dp) :: dalpha
        real(dp) :: derphi_aj
        real(dp) :: derphi_lo
        real(dp) :: phi_aj
        real(dp) :: phi_hi
        real(dp) :: phi_lo
        real(dp) :: phi_rec
        real(dp) :: qchk
        integer :: i
        logical :: accept
        logical :: interp_ok

        a_lo = a_lo_in
        a_hi = a_hi_in
        phi_lo = phi_lo_in
        phi_hi = phi_hi_in
        derphi_lo = derphi_lo_in
        phi_rec = phi0
        a_rec = 0.0_dp
        success = .false.
        step = 0.0_dp
        phi_star = phi0

        do i = 0, 10
            dalpha = a_hi - a_lo
            if (dalpha < 0.0_dp) then
                a = a_hi
                b = a_lo
            else
                a = a_lo
                b = a_hi
            end if

            interp_ok = .false.
            if (i > 0) then
                cchk = delta1 * dalpha
                call cubic_minimizer(a_lo, phi_lo, derphi_lo, a_hi, phi_hi, &
                    a_rec, phi_rec, a_j, interp_ok)
                if (interp_ok) then
                    if (a_j > b - cchk .or. a_j < a + cchk) interp_ok = .false.
                end if
            end if
            if (.not. interp_ok) then
                qchk = delta2 * dalpha
                call quadratic_minimizer(a_lo, phi_lo, derphi_lo, a_hi, phi_hi, &
                    a_j, interp_ok)
                if (interp_ok) then
                    if (a_j > b - qchk .or. a_j < a + qchk) interp_ok = .false.
                end if
                if (.not. interp_ok) a_j = a_lo + 0.5_dp * dalpha
            end if

            phi_aj = phi(a_j)
            if (phi_aj > phi0 + c1 * a_j * derphi0 .or. phi_aj >= phi_lo) then
                phi_rec = phi_hi
                a_rec = a_hi
                a_hi = a_j
                phi_hi = phi_aj
            else
                derphi_aj = derphi(a_j)
                if (abs(derphi_aj) <= -c2 * derphi0) then
                    accept = .true.
                    if (present(extra_condition)) accept = extra_condition(a_j, phi_aj)
                    if (accept) then
                        step = a_j
                        phi_star = phi_aj
                        success = .true.
                        return
                    end if
                end if
                if (derphi_aj * (a_hi - a_lo) >= 0.0_dp) then
                    phi_rec = phi_hi
                    a_rec = a_hi
                    a_hi = a_lo
                    phi_hi = phi_lo
                else
                    phi_rec = phi_lo
                    a_rec = a_lo
                end if
                a_lo = a_j
                phi_lo = phi_aj
                derphi_lo = derphi_aj
            end if
        end do
    end subroutine wolfe_zoom

    subroutine cubic_minimizer(a, fa, fpa, b, fb, c, fc, xmin, success)
        real(dp), intent(in) :: a !! first interpolation abscissa
        real(dp), intent(in) :: fa !! objective at a
        real(dp), intent(in) :: fpa !! derivative at a
        real(dp), intent(in) :: b !! second interpolation abscissa
        real(dp), intent(in) :: fb !! objective at b
        real(dp), intent(in) :: c !! third interpolation abscissa
        real(dp), intent(in) :: fc !! objective at c
        real(dp), intent(out) :: xmin !! cubic minimizer when successful
        logical, intent(out) :: success !! true when a finite minimizer exists

        real(dp) :: aa
        real(dp) :: bb
        real(dp) :: db
        real(dp) :: dc
        real(dp) :: denom
        real(dp) :: radical
        real(dp) :: rhs1
        real(dp) :: rhs2

        db = b - a
        dc = c - a
        denom = (db * dc)**2 * (db - dc)
        success = .false.
        xmin = a
        if (same_real(denom, 0.0_dp)) return
        rhs1 = fb - fa - fpa * db
        rhs2 = fc - fa - fpa * dc
        aa = (dc**2 * rhs1 - db**2 * rhs2) / denom
        bb = (-dc**3 * rhs1 + db**3 * rhs2) / denom
        if (same_real(aa, 0.0_dp)) return
        radical = bb * bb - 3.0_dp * aa * fpa
        if (radical < 0.0_dp) return
        xmin = a + (-bb + sqrt(radical)) / (3.0_dp * aa)
        success = ieee_is_finite(xmin)
    end subroutine cubic_minimizer

    subroutine quadratic_minimizer(a, fa, fpa, b, fb, xmin, success)
        real(dp), intent(in) :: a !! first interpolation abscissa
        real(dp), intent(in) :: fa !! objective at a
        real(dp), intent(in) :: fpa !! derivative at a
        real(dp), intent(in) :: b !! second interpolation abscissa
        real(dp), intent(in) :: fb !! objective at b
        real(dp), intent(out) :: xmin !! quadratic minimizer when successful
        logical, intent(out) :: success !! true when a finite minimizer exists

        real(dp) :: bb
        real(dp) :: db

        db = b - a
        success = .false.
        xmin = a
        if (same_real(db, 0.0_dp)) return
        bb = (fb - fa - fpa * db) / (db * db)
        if (same_real(bb, 0.0_dp)) return
        xmin = a - fpa / (2.0_dp * bb)
        success = ieee_is_finite(xmin)
    end subroutine quadratic_minimizer

    subroutine dcstep(stx, fx, dx, sty, fy, dy, stp, fp, dp_value, brackt, &
            stpmin, stpmax, gtest)
        real(dp), intent(inout) :: stx !! best step endpoint
        real(dp), intent(inout) :: fx !! objective at stx
        real(dp), intent(inout) :: dx !! derivative at stx
        real(dp), intent(inout) :: sty !! opposite bracket endpoint
        real(dp), intent(inout) :: fy !! objective at sty
        real(dp), intent(inout) :: dy !! derivative at sty
        real(dp), intent(inout) :: stp !! current/new trial step
        real(dp), intent(in) :: fp !! objective at current stp
        real(dp), intent(in) :: dp_value !! derivative at current stp
        logical, intent(inout) :: brackt !! whether a minimizer is bracketed
        real(dp), intent(in) :: stpmin !! lower bound on next step
        real(dp), intent(in) :: stpmax !! upper bound on next step
        real(dp), intent(in), optional :: gtest !! stage-one modified derivative offset

        real(dp) :: dp_use
        real(dp) :: dx_use
        real(dp) :: dy_use
        real(dp) :: fp_use
        real(dp) :: fx_use
        real(dp) :: fy_use
        real(dp) :: gamma
        real(dp) :: p
        real(dp) :: q
        real(dp) :: r
        real(dp) :: s
        real(dp) :: sgnd
        real(dp) :: stpc
        real(dp) :: stpf
        real(dp) :: stpq
        real(dp) :: theta

        fx_use = fx
        fy_use = fy
        fp_use = fp
        dx_use = dx
        dy_use = dy
        dp_use = dp_value
        if (present(gtest)) then
            fx_use = fx - stx * gtest
            fy_use = fy - sty * gtest
            fp_use = fp - stp * gtest
            dx_use = dx - gtest
            dy_use = dy - gtest
            dp_use = dp_value - gtest
        end if

        if (same_real(dp_use, 0.0_dp) .or. same_real(dx_use, 0.0_dp)) then
            sgnd = 0.0_dp
        else
            sgnd = sign(1.0_dp, dp_use) * sign(1.0_dp, dx_use)
        end if

        if (fp_use > fx_use) then
            theta = 3.0_dp * (fx_use - fp_use) / (stp - stx) + dx_use + dp_use
            s = max(abs(theta), abs(dx_use), abs(dp_use))
            gamma = s * sqrt(max(0.0_dp, (theta / s)**2 - &
                (dx_use / s) * (dp_use / s)))
            if (stp < stx) gamma = -gamma
            p = (gamma - dx_use) + theta
            q = ((gamma - dx_use) + gamma) + dp_use
            r = p / q
            stpc = stx + r * (stp - stx)
            stpq = stx + ((dx_use / ((fx_use - fp_use) / (stp - stx) + &
                dx_use)) / 2.0_dp) * (stp - stx)
            if (abs(stpc - stx) <= abs(stpq - stx)) then
                stpf = stpc
            else
                stpf = stpc + (stpq - stpc) / 2.0_dp
            end if
            brackt = .true.
        else if (sgnd < 0.0_dp) then
            theta = 3.0_dp * (fx_use - fp_use) / (stp - stx) + dx_use + dp_use
            s = max(abs(theta), abs(dx_use), abs(dp_use))
            gamma = s * sqrt(max(0.0_dp, (theta / s)**2 - &
                (dx_use / s) * (dp_use / s)))
            if (stp > stx) gamma = -gamma
            p = (gamma - dp_use) + theta
            q = ((gamma - dp_use) + gamma) + dx_use
            r = p / q
            stpc = stp + r * (stx - stp)
            stpq = stp + (dp_use / (dp_use - dx_use)) * (stx - stp)
            if (abs(stpc - stp) > abs(stpq - stp)) then
                stpf = stpc
            else
                stpf = stpq
            end if
            brackt = .true.
        else if (abs(dp_use) < abs(dx_use)) then
            theta = 3.0_dp * (fx_use - fp_use) / (stp - stx) + dx_use + dp_use
            s = max(abs(theta), abs(dx_use), abs(dp_use))
            gamma = s * sqrt(max(0.0_dp, (theta / s)**2 - &
                (dx_use / s) * (dp_use / s)))
            if (stp > stx) gamma = -gamma
            p = (gamma - dp_use) + theta
            q = (gamma + (dx_use - dp_use)) + gamma
            r = p / q
            if (r < 0.0_dp .and. .not. same_real(gamma, 0.0_dp)) then
                stpc = stp + r * (stx - stp)
            else if (stp > stx) then
                stpc = stpmax
            else
                stpc = stpmin
            end if
            stpq = stp + (dp_use / (dp_use - dx_use)) * (stx - stp)
            if (brackt) then
                if (abs(stpc - stp) < abs(stpq - stp)) then
                    stpf = stpc
                else
                    stpf = stpq
                end if
                if (stp > stx) then
                    stpf = min(stp + 0.66_dp * (sty - stp), stpf)
                else
                    stpf = max(stp + 0.66_dp * (sty - stp), stpf)
                end if
            else
                if (abs(stpc - stp) > abs(stpq - stp)) then
                    stpf = stpc
                else
                    stpf = stpq
                end if
                stpf = min(stpmax, max(stpmin, stpf))
            end if
        else
            if (brackt) then
                theta = 3.0_dp * (fp_use - fy_use) / (sty - stp) + dy_use + dp_use
                s = max(abs(theta), abs(dy_use), abs(dp_use))
                gamma = s * sqrt(max(0.0_dp, (theta / s)**2 - &
                    (dy_use / s) * (dp_use / s)))
                if (stp > sty) gamma = -gamma
                p = (gamma - dp_use) + theta
                q = ((gamma - dp_use) + gamma) + dy_use
                r = p / q
                stpc = stp + r * (sty - stp)
                stpf = stpc
            else if (stp > stx) then
                stpf = stpmax
            else
                stpf = stpmin
            end if
        end if

        if (fp_use > fx_use) then
            sty = stp
            fy_use = fp_use
            dy_use = dp_use
        else
            if (sgnd < 0.0_dp) then
                sty = stx
                fy_use = fx_use
                dy_use = dx_use
            end if
            stx = stp
            fx_use = fp_use
            dx_use = dp_use
        end if
        stp = stpf

        if (present(gtest)) then
            fx = fx_use + stx * gtest
            fy = fy_use + sty * gtest
            dx = dx_use + gtest
            dy = dy_use + gtest
        else
            fx = fx_use
            fy = fy_use
            dx = dx_use
            dy = dy_use
        end if
    end subroutine dcstep

    pure logical function same_real(a, b) result(equal)
        real(dp), intent(in) :: a !! first value
        real(dp), intent(in) :: b !! second value

        equal = a <= b .and. a >= b
    end function same_real

end module scifort_optimize_wolfe
