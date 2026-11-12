! SPDX-License-Identifier: BSD-3-Clause
! Copyright (c) 2026 SciFort contributors
! Copyright (c) 2020-2025, Zaikun ZHANG
! Copyright (c) 2001-2002 Enthought, Inc. 2003, SciPy Developers.
!
! Internal COBYLA numerical kernel for the scipy.optimize translation.
!
! The algorithm, trust-region LP subproblem, simplex geometry updates, and
! stopping rules are adapted primarily from SciPy 1.17.0
! scipy/_lib/pyprima/cobyla, a faithful translation of Zaikun Zhang's modern
! Fortran PRIMA reference implementation of Powell's COBYLA algorithm. Both
! sources are distributed under BSD-3-Clause. The implementation here is an
! independent modern-Fortran adaptation for SciFort. It intentionally keeps
! all constraints in the PRIMA convention c(x) <= 0; the public SciFort wrapper
! converts its SciPy-like equality, inequality, and bound conventions.

module scifort_optimize_cobyla_prima
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    implicit none
    private

    real(dp), parameter :: funcmax = 1.0e30_dp
    real(dp), parameter :: constrmax = 1.0e30_dp
    real(dp), parameter :: eps = epsilon(1.0_dp)
    real(dp), parameter :: realmax = huge(1.0_dp)
    real(dp), parameter :: realmin = tiny(1.0_dp)

    integer, parameter, public :: cobyla_small_tr_radius = 0
    integer, parameter, public :: cobyla_ftarget_achieved = 1
    integer, parameter, public :: cobyla_trsubp_failed = 2
    integer, parameter, public :: cobyla_maxfun_reached = 3
    integer, parameter, public :: cobyla_damaging_rounding = 7
    integer, parameter, public :: cobyla_maxtr_reached = 20

    abstract interface
        subroutine cobyla_evaluator(x, f, c)
            import dp
            implicit none
            real(dp), intent(in) :: x(:)
            real(dp), intent(out) :: f
            real(dp), intent(out) :: c(:)
        end subroutine cobyla_evaluator
    end interface

    public :: cobyla_evaluator
    public :: prima_cobyla

contains

    subroutine prima_cobyla(calcfc, x0, m, rhobeg, rhoend, maxfun, ctol, ftarget, &
            x, f, constr, cstrv, nf, info, nit)
        procedure(cobyla_evaluator) :: calcfc
        real(dp), intent(in) :: x0(:)
        integer, intent(in) :: m
        real(dp), intent(in) :: rhobeg
        real(dp), intent(in) :: rhoend
        integer, intent(in) :: maxfun
        real(dp), intent(in) :: ctol
        real(dp), intent(in) :: ftarget
        real(dp), intent(out) :: x(:)
        real(dp), intent(out) :: f
        real(dp), intent(out) :: constr(:)
        real(dp), intent(out) :: cstrv
        integer, intent(out) :: nf
        integer, intent(out) :: info
        integer, intent(out) :: nit

        real(dp), allocatable :: a(:, :)
        real(dp), allocatable :: chist(:)
        real(dp), allocatable :: conmat(:, :)
        real(dp), allocatable :: conhist(:, :)
        real(dp), allocatable :: cval(:)
        real(dp), allocatable :: d(:)
        real(dp), allocatable :: distsq(:)
        real(dp), allocatable :: fhist(:)
        real(dp), allocatable :: fval(:)
        real(dp), allocatable :: g(:)
        real(dp), allocatable :: modelc(:)
        real(dp), allocatable :: sim(:, :)
        real(dp), allocatable :: simi(:, :)
        real(dp), allocatable :: trialc(:)
        real(dp), allocatable :: trialx(:)
        real(dp), allocatable :: xhist(:, :)
        real(dp) :: actrem
        real(dp) :: cpen
        real(dp) :: delta
        real(dp) :: delbar
        real(dp) :: dnorm
        real(dp) :: gamma1
        real(dp) :: gamma2
        real(dp) :: gamma3
        real(dp) :: newrho
        real(dp) :: prerec
        real(dp) :: preref
        real(dp) :: prerem
        real(dp) :: ratio
        real(dp) :: rho
        real(dp) :: trialcv
        real(dp) :: trialf
        integer :: histn
        integer :: i
        integer :: j
        integer :: jdrop_geo
        integer :: jdrop_tr
        integer :: maxfun_use
        integer :: maxtr
        integer :: n
        integer :: subinfo
        integer :: tr
        logical :: adequate_geo
        logical :: bad_trstep
        logical :: improve_geo
        logical :: reduce_rho
        logical :: reused
        logical :: shortd
        logical :: trfail
        logical :: ximproved

        n = size(x0)
        nf = 0
        nit = 0
        info = cobyla_maxtr_reached
        f = funcmax
        cstrv = constrmax
        if (size(x) /= n .or. size(constr) /= m) return
        x = x0
        if (m > 0) constr = constrmax
        if (n <= 0 .or. m < 0) return
        if (.not. all(ieee_is_finite(x0))) return
        if (.not. ieee_is_finite(rhobeg) .or. .not. ieee_is_finite(rhoend)) return
        if (rhobeg <= 0.0_dp .or. rhoend <= 0.0_dp .or. rhobeg < rhoend) return

        maxfun_use = max(maxfun, n + 2)
        allocate(a(n, m), chist(maxfun_use), conmat(m, n + 1), &
            conhist(m, maxfun_use), cval(n + 1), d(n), distsq(n + 1), &
            fhist(maxfun_use), fval(n + 1), g(n), modelc(m), &
            sim(n, n + 1), simi(n, n), trialc(m), trialx(n), &
            xhist(n, maxfun_use))
        a = 0.0_dp
        chist = constrmax
        if (m > 0) conmat = constrmax
        if (m > 0) conhist = constrmax
        cval = constrmax
        d = 0.0_dp
        fhist = funcmax
        fval = funcmax
        sim = 0.0_dp
        simi = 0.0_dp
        xhist = 0.0_dp
        histn = 0

        ! PRIMA initialization: the last column is the pole and the first N
        ! columns are displacements from the pole to the other vertices.
        do i = 1, n
            sim(i, i) = rhobeg
        end do
        sim(:, n + 1) = x0

        call evaluate_point(calcfc, sim(:, n + 1), trialf, trialc)
        nf = nf + 1
        trialcv = constraint_violation(trialc)
        fval(n + 1) = trialf
        if (m > 0) conmat(:, n + 1) = trialc
        cval(n + 1) = trialcv
        call save_history(sim(:, n + 1), trialf, trialc, trialcv, xhist, fhist, &
            conhist, chist, histn)
        subinfo = check_break(maxfun_use, nf, trialcv, ctol, trialf, ftarget)
        if (subinfo /= -999) then
            info = subinfo
            call select_history(xhist, fhist, conhist, chist, histn, 1.0e8_dp, ctol, &
                x, f, constr, cstrv)
            return
        end if

        do j = 1, n
            trialx = sim(:, n + 1)
            trialx(j) = trialx(j) + rhobeg
            call evaluate_point(calcfc, trialx, trialf, trialc)
            nf = nf + 1
            trialcv = constraint_violation(trialc)
            fval(j) = trialf
            if (m > 0) conmat(:, j) = trialc
            cval(j) = trialcv
            call save_history(trialx, trialf, trialc, trialcv, xhist, fhist, &
                conhist, chist, histn)

            subinfo = check_break(maxfun_use, nf, trialcv, ctol, trialf, ftarget)
            if (subinfo /= -999) then
                info = subinfo
                call select_history(xhist, fhist, conhist, chist, histn, 1.0e8_dp, ctol, &
                    x, f, constr, cstrv)
                return
            end if

            if (fval(j) < fval(n + 1)) then
                call swap_scalar(fval(j), fval(n + 1))
                call swap_scalar(cval(j), cval(n + 1))
                if (m > 0) call swap_columns(conmat, j, n + 1)
                sim(:, n + 1) = trialx
                sim(j, 1:j) = -rhobeg
            end if
        end do

        call invert_matrix(sim(:, 1:n), simi, adequate_geo)
        if (.not. adequate_geo) then
            info = cobyla_damaging_rounding
            call select_history(xhist, fhist, conhist, chist, histn, 1.0e8_dp, ctol, &
                x, f, constr, cstrv)
            return
        end if

        rho = rhobeg
        delta = rhobeg
        cpen = max(eps, min(1.0e3_dp, fc_ratio(conmat, fval)))
        shortd = .false.
        ratio = -1.0_dp
        jdrop_tr = 0
        gamma1 = 0.5_dp
        gamma2 = 2.0_dp
        gamma3 = max(1.0_dp, min(0.75_dp * gamma2, 1.5_dp))
        maxtr = 10 * maxfun_use
        info = cobyla_maxtr_reached

        do tr = 1, maxtr
            nit = tr - 1
            call get_cpen(conmat, cval, fval, sim, simi, delta, cpen, subinfo)
            if (subinfo == cobyla_damaging_rounding) then
                info = subinfo
                exit
            end if

            call update_pole(cpen, conmat, cval, fval, sim, simi, subinfo)
            if (subinfo == cobyla_damaging_rounding) then
                info = subinfo
                exit
            end if

            adequate_geo = simplex_geometry_adequate(sim, delta)
            call linear_models(conmat, fval, simi, g, a)
            call trstlp(a, -conmat(:, n + 1), delta, g, d)
            dnorm = min(delta, vecnorm(d))
            shortd = dnorm <= 0.1_dp * rho

            preref = -dot_product(d, g)
            if (m > 0) then
                modelc = conmat(:, n + 1) + matmul(transpose(a), d)
                prerec = cval(n + 1) - constraint_violation(modelc)
            else
                prerec = 0.0_dp
            end if
            prerem = preref + cpen * prerec
            trfail = .not. (prerem > 1.0e-6_dp * min(cpen, 1.0_dp) * rho)

            if (shortd .or. trfail) then
                delta = 0.1_dp * delta
                if (delta <= gamma3 * rho) delta = rho
            else
                trialx = sim(:, n + 1) + d
                call interpolation_distances(trialx, sim, distsq)
                j = minloc(distsq, dim=1)
                reused = distsq(j) <= (1.0e-4_dp * rhoend)**2
                if (reused) then
                    trialf = fval(j)
                    if (m > 0) trialc = conmat(:, j)
                    trialcv = cval(j)
                else
                    call evaluate_point(calcfc, trialx, trialf, trialc)
                    nf = nf + 1
                    trialcv = constraint_violation(trialc)
                    call save_history(trialx, trialf, trialc, trialcv, xhist, fhist, &
                        conhist, chist, histn)
                end if

                actrem = (fval(n + 1) + cpen * cval(n + 1)) - &
                    (trialf + cpen * trialcv)
                ratio = reduction_ratio(actrem, prerem, 0.1_dp)
                if (dnorm > 0.0_dp) then
                    delta = trust_radius(delta, dnorm, 0.1_dp, 0.7_dp, gamma1, gamma2, ratio)
                    if (delta <= gamma3 * rho) delta = rho
                end if
                ximproved = actrem > 0.0_dp
                jdrop_tr = set_drop_tr(ximproved, d, delta, rho, sim, simi)
                call update_xfc(jdrop_tr, trialc, cpen, trialcv, d, trialf, conmat, cval, &
                    fval, sim, simi, subinfo)
                if (subinfo == cobyla_damaging_rounding) then
                    info = subinfo
                    exit
                end if

                if (.not. reused) then
                    subinfo = check_break(maxfun_use, nf, trialcv, ctol, trialf, ftarget)
                    if (subinfo /= -999) then
                        info = subinfo
                        exit
                    end if
                end if
            end if

            bad_trstep = shortd .or. trfail .or. ratio <= 0.0_dp .or. jdrop_tr == 0
            improve_geo = bad_trstep .and. .not. adequate_geo
            reduce_rho = bad_trstep .and. adequate_geo .and. max(delta, dnorm) <= rho

            if (improve_geo .and. .not. simplex_geometry_adequate(sim, delta)) then
                jdrop_geo = farthest_simplex_column(sim)
                delbar = 0.5_dp * delta
                call geometry_step(jdrop_geo, conmat, cpen, delbar, fval, simi, d)
                trialx = sim(:, n + 1) + d
                call interpolation_distances(trialx, sim, distsq)
                j = minloc(distsq, dim=1)
                reused = distsq(j) <= (1.0e-4_dp * rhoend)**2
                if (reused) then
                    trialf = fval(j)
                    if (m > 0) trialc = conmat(:, j)
                    trialcv = cval(j)
                else
                    call evaluate_point(calcfc, trialx, trialf, trialc)
                    nf = nf + 1
                    trialcv = constraint_violation(trialc)
                    call save_history(trialx, trialf, trialc, trialcv, xhist, fhist, &
                        conhist, chist, histn)
                end if
                call update_xfc(jdrop_geo, trialc, cpen, trialcv, d, trialf, conmat, &
                    cval, fval, sim, simi, subinfo)
                if (subinfo == cobyla_damaging_rounding) then
                    info = subinfo
                    exit
                end if
                if (.not. reused) then
                    subinfo = check_break(maxfun_use, nf, trialcv, ctol, trialf, ftarget)
                    if (subinfo /= -999) then
                        info = subinfo
                        exit
                    end if
                end if
            end if

            if (reduce_rho) then
                if (rho <= rhoend) then
                    info = cobyla_small_tr_radius
                    exit
                end if
                newrho = reduce_rho_value(rho, rhoend)
                delta = max(0.5_dp * rho, newrho)
                rho = newrho
                cpen = max(eps, min(cpen, fc_ratio(conmat, fval)))
                call update_pole(cpen, conmat, cval, fval, sim, simi, subinfo)
                if (subinfo == cobyla_damaging_rounding) then
                    info = subinfo
                    exit
                end if
            end if
        end do

        nit = min(tr, maxtr)

        ! PRIMA tries the last short trust-region step once before returning.
        trialx = sim(:, n + 1) + d
        if (info == cobyla_small_tr_radius .and. shortd .and. &
                vecnorm(trialx - sim(:, n + 1)) > 1.0e-3_dp * rhoend .and. &
                nf < maxfun_use) then
            call evaluate_point(calcfc, trialx, trialf, trialc)
            nf = nf + 1
            trialcv = constraint_violation(trialc)
            call save_history(trialx, trialf, trialc, trialcv, xhist, fhist, conhist, &
                chist, histn)
        end if

        call select_history(xhist, fhist, conhist, chist, histn, max(cpen, 1.0e8_dp), &
            ctol, x, f, constr, cstrv)
    end subroutine prima_cobyla

    subroutine evaluate_point(calcfc, x, f, c)
        procedure(cobyla_evaluator) :: calcfc
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: f
        real(dp), intent(out) :: c(:)
        integer :: i

        call calcfc(x, f, c)
        if (ieee_is_nan(f) .or. f > funcmax) f = funcmax
        if (f < -realmax) f = -realmax
        do i = 1, size(c)
            if (ieee_is_nan(c(i)) .or. c(i) > constrmax) then
                c(i) = constrmax
            else if (c(i) < -constrmax) then
                c(i) = -constrmax
            end if
        end do
    end subroutine evaluate_point

    pure real(dp) function constraint_violation(c) result(cv)
        real(dp), intent(in) :: c(:)
        if (size(c) <= 0) then
            cv = 0.0_dp
        else
            cv = max(0.0_dp, maxval(c))
        end if
    end function constraint_violation

    subroutine save_history(x, f, c, cv, xhist, fhist, conhist, chist, histn)
        real(dp), intent(in) :: x(:)
        real(dp), intent(in) :: f
        real(dp), intent(in) :: c(:)
        real(dp), intent(in) :: cv
        real(dp), intent(inout) :: xhist(:, :)
        real(dp), intent(inout) :: fhist(:)
        real(dp), intent(inout) :: conhist(:, :)
        real(dp), intent(inout) :: chist(:)
        integer, intent(inout) :: histn

        if (histn >= size(fhist)) return
        histn = histn + 1
        xhist(:, histn) = x
        fhist(histn) = f
        if (size(c) > 0) conhist(:, histn) = c
        chist(histn) = cv
    end subroutine save_history

    integer function check_break(maxfun, nf, cv, ctol, f, ftarget) result(info)
        integer, intent(in) :: maxfun
        integer, intent(in) :: nf
        real(dp), intent(in) :: cv
        real(dp), intent(in) :: ctol
        real(dp), intent(in) :: f
        real(dp), intent(in) :: ftarget

        info = -999
        if (cv <= ctol .and. f <= ftarget) then
            info = cobyla_ftarget_achieved
        else if (nf >= maxfun) then
            info = cobyla_maxfun_reached
        end if
    end function check_break

    subroutine select_history(xhist, fhist, conhist, chist, nhist, cweight, ctol, &
            x, f, c, cv)
        real(dp), intent(in) :: xhist(:, :)
        real(dp), intent(in) :: fhist(:)
        real(dp), intent(in) :: conhist(:, :)
        real(dp), intent(in) :: chist(:)
        integer, intent(in) :: nhist
        real(dp), intent(in) :: cweight
        real(dp), intent(in) :: ctol
        real(dp), intent(out) :: x(:)
        real(dp), intent(out) :: f
        real(dp), intent(out) :: c(:)
        real(dp), intent(out) :: cv

        real(dp), allocatable :: shifted(:)
        real(dp) :: cmin
        real(dp) :: cref
        real(dp) :: phi
        real(dp) :: phibest
        integer :: i
        integer :: kbest

        if (nhist <= 0) then
            x = 0.0_dp
            f = funcmax
            if (size(c) > 0) c = constrmax
            cv = constrmax
            return
        end if

        allocate(shifted(nhist))
        shifted = max(chist(1:nhist) - max(ctol, 0.0_dp), 0.0_dp)
        cmin = minval(shifted)
        cref = max(eps, 2.0_dp * cmin)
        kbest = 0
        phibest = realmax
        do i = 1, nhist
            if (fhist(i) >= funcmax .or. chist(i) >= constrmax) cycle
            if (shifted(i) > cref) cycle
            phi = max(fhist(i), -realmax) + max(cweight, 0.0_dp) * shifted(i)
            if (kbest == 0 .or. phi < phibest) then
                kbest = i
                phibest = phi
            else if (exact_equal_real(phi, phibest)) then
                if (shifted(i) < shifted(kbest) .or. &
                        (exact_equal_real(shifted(i), shifted(kbest)) .and. &
                        fhist(i) < fhist(kbest)) .or. &
                        (exact_equal_real(shifted(i), shifted(kbest)) .and. &
                        exact_equal_real(fhist(i), fhist(kbest)) .and. &
                        chist(i) < chist(kbest))) then
                    kbest = i
                    phibest = phi
                end if
            end if
        end do
        if (kbest == 0) kbest = nhist
        x = xhist(:, kbest)
        f = fhist(kbest)
        if (size(c) > 0) c = conhist(:, kbest)
        cv = chist(kbest)
    end subroutine select_history

    pure logical function simplex_geometry_adequate(sim, delta) result(ok)
        real(dp), intent(in) :: sim(:, :)
        real(dp), intent(in) :: delta
        integer :: j
        integer :: n

        n = size(sim, 1)
        ok = .true.
        do j = 1, n
            if (dot_product(sim(:, j), sim(:, j)) > 4.0_dp * delta * delta) then
                ok = .false.
                exit
            end if
        end do
    end function simplex_geometry_adequate

    integer function farthest_simplex_column(sim) result(jdrop)
        real(dp), intent(in) :: sim(:, :)
        real(dp) :: best
        real(dp) :: dist
        integer :: j
        integer :: n

        n = size(sim, 1)
        jdrop = 1
        best = -1.0_dp
        do j = 1, n
            dist = dot_product(sim(:, j), sim(:, j))
            if (dist > best) then
                best = dist
                jdrop = j
            end if
        end do
    end function farthest_simplex_column

    subroutine interpolation_distances(x, sim, distsq)
        real(dp), intent(in) :: x(:)
        real(dp), intent(in) :: sim(:, :)
        real(dp), intent(out) :: distsq(:)
        real(dp), allocatable :: dx(:)
        integer :: j
        integer :: n

        n = size(x)
        allocate(dx(n))
        dx = x - sim(:, n + 1)
        distsq(n + 1) = dot_product(dx, dx)
        do j = 1, n
            dx = x - (sim(:, n + 1) + sim(:, j))
            distsq(j) = dot_product(dx, dx)
        end do
    end subroutine interpolation_distances

    subroutine linear_models(conmat, fval, simi, g, a)
        real(dp), intent(in) :: conmat(:, :)
        real(dp), intent(in) :: fval(:)
        real(dp), intent(in) :: simi(:, :)
        real(dp), intent(out) :: g(:)
        real(dp), intent(out) :: a(:, :)
        real(dp), allocatable :: diff(:)
        integer :: i
        integer :: m
        integer :: n

        n = size(simi, 1)
        m = size(conmat, 1)
        g = matmul(fval(1:n) - fval(n + 1), simi)
        if (m <= 0) return
        allocate(diff(n))
        do i = 1, m
            diff = conmat(i, 1:n) - conmat(i, n + 1)
            a(:, i) = matmul(diff, simi)
        end do
    end subroutine linear_models

    subroutine get_cpen(conmat_in, cval_in, fval_in, sim_in, simi_in, delta, &
            cpen, info)
        real(dp), intent(in) :: conmat_in(:, :)
        real(dp), intent(in) :: cval_in(:)
        real(dp), intent(in) :: fval_in(:)
        real(dp), intent(in) :: sim_in(:, :)
        real(dp), intent(in) :: simi_in(:, :)
        real(dp), intent(in) :: delta
        real(dp), intent(inout) :: cpen
        integer, intent(out) :: info

        real(dp), allocatable :: a(:, :)
        real(dp), allocatable :: conmat(:, :)
        real(dp), allocatable :: cval(:)
        real(dp), allocatable :: d(:)
        real(dp), allocatable :: fval(:)
        real(dp), allocatable :: g(:)
        real(dp), allocatable :: modelc(:)
        real(dp), allocatable :: sim(:, :)
        real(dp), allocatable :: simi(:, :)
        real(dp) :: prerec
        real(dp) :: preref
        integer :: iter
        integer :: m
        integer :: n

        n = size(sim_in, 1)
        m = size(conmat_in, 1)
        allocate(a(n, m), conmat(m, n + 1), cval(n + 1), d(n), fval(n + 1), &
            g(n), modelc(m), sim(n, n + 1), simi(n, n))
        if (m > 0) conmat = conmat_in
        cval = cval_in
        fval = fval_in
        sim = sim_in
        simi = simi_in
        info = -999

        do iter = 1, n + 1
            call update_pole(cpen, conmat, cval, fval, sim, simi, info)
            if (info == cobyla_damaging_rounding) return
            call linear_models(conmat, fval, simi, g, a)
            call trstlp(a, -conmat(:, n + 1), delta, g, d)
            preref = -dot_product(d, g)
            if (m > 0) then
                modelc = conmat(:, n + 1) + matmul(transpose(a), d)
                prerec = cval(n + 1) - constraint_violation(modelc)
            else
                prerec = 0.0_dp
            end if
            if (.not. (prerec > 0.0_dp .and. preref < 0.0_dp)) exit
            cpen = max(cpen, min(-2.0_dp * preref / prerec, realmax))
            if (find_pole(cpen, cval, fval) == n + 1) exit
        end do
        info = -999
    end subroutine get_cpen

    pure real(dp) function fc_ratio(conmat, fval) result(r)
        real(dp), intent(in) :: conmat(:, :)
        real(dp), intent(in) :: fval(:)
        real(dp) :: cmax
        real(dp) :: cmin
        real(dp) :: denom
        real(dp) :: fmax
        real(dp) :: fmin
        integer :: i
        logical :: have

        if (size(conmat, 1) <= 0) then
            r = 0.0_dp
            return
        end if
        fmin = minval(fval)
        fmax = maxval(fval)
        denom = realmax
        have = .false.
        do i = 1, size(conmat, 1)
            cmin = minval(-conmat(i, :))
            cmax = maxval(-conmat(i, :))
            if (cmin < 0.5_dp * cmax) then
                have = .true.
                denom = min(denom, max(cmax, 0.0_dp) - cmin)
            end if
        end do
        if (have .and. fmin < fmax .and. denom > 0.0_dp .and. ieee_is_finite(denom)) then
            r = max(0.0_dp, (fmax - fmin) / denom)
        else
            r = 0.0_dp
        end if
    end function fc_ratio

    integer function find_pole(cpen, cval, fval) result(jopt)
        real(dp), intent(in) :: cpen
        real(dp), intent(in) :: cval(:)
        real(dp), intent(in) :: fval(:)
        real(dp), allocatable :: phi(:)
        real(dp) :: phimin
        integer :: j
        integer :: n

        n = size(fval) - 1
        allocate(phi(n + 1))
        phi = fval + cpen * cval
        phimin = minval(phi)
        jopt = n + 1
        if (phimin < phi(n + 1) .or. &
                any((cval < cval(n + 1)) .and. (phi <= phi(n + 1)))) then
            do j = 1, n + 1
                if (phi(j) <= phimin) then
                    if (jopt == n + 1 .or. phi(jopt) > phimin .or. cval(j) < cval(jopt)) &
                        jopt = j
                end if
            end do
        end if
    end function find_pole

    subroutine update_pole(cpen, conmat, cval, fval, sim, simi, info)
        real(dp), intent(in) :: cpen
        real(dp), intent(inout) :: conmat(:, :)
        real(dp), intent(inout) :: cval(:)
        real(dp), intent(inout) :: fval(:)
        real(dp), intent(inout) :: sim(:, :)
        real(dp), intent(inout) :: simi(:, :)
        integer, intent(out) :: info

        real(dp), allocatable :: sim_j(:)
        real(dp), allocatable :: test_inv(:, :)
        integer :: j
        integer :: jopt
        integer :: n
        logical :: ok

        n = size(sim, 1)
        info = -999
        jopt = find_pole(cpen, cval, fval)
        if (jopt <= n) then
            allocate(sim_j(n))
            sim(:, n + 1) = sim(:, n + 1) + sim(:, jopt)
            sim_j = sim(:, jopt)
            sim(:, jopt) = 0.0_dp
            do j = 1, n
                sim(:, j) = sim(:, j) - sim_j
            end do
        end if

        allocate(test_inv(n, n))
        call invert_matrix(sim(:, 1:n), test_inv, ok)
        if (.not. ok .or. inverse_error(test_inv, sim(:, 1:n)) > 1.0_dp) then
            info = cobyla_damaging_rounding
            return
        end if
        simi = test_inv

        if (jopt <= n) then
            call swap_scalar(fval(jopt), fval(n + 1))
            call swap_scalar(cval(jopt), cval(n + 1))
            if (size(conmat, 1) > 0) call swap_columns(conmat, jopt, n + 1)
        end if
    end subroutine update_pole

    subroutine update_xfc(jdrop, constr, cpen, cstrv, d, f, conmat, cval, fval, sim, &
            simi, info)
        integer, intent(in) :: jdrop
        real(dp), intent(in) :: constr(:)
        real(dp), intent(in) :: cpen
        real(dp), intent(in) :: cstrv
        real(dp), intent(in) :: d(:)
        real(dp), intent(in) :: f
        real(dp), intent(inout) :: conmat(:, :)
        real(dp), intent(inout) :: cval(:)
        real(dp), intent(inout) :: fval(:)
        real(dp), intent(inout) :: sim(:, :)
        real(dp), intent(inout) :: simi(:, :)
        integer, intent(out) :: info

        real(dp), allocatable :: test_inv(:, :)
        integer :: j
        integer :: n
        logical :: ok

        info = -999
        n = size(sim, 1)
        if (jdrop <= 0) return
        if (jdrop <= n) then
            sim(:, jdrop) = d
        else if (jdrop == n + 1) then
            sim(:, n + 1) = sim(:, n + 1) + d
            do j = 1, n
                sim(:, j) = sim(:, j) - d
            end do
        else
            info = cobyla_damaging_rounding
            return
        end if

        allocate(test_inv(n, n))
        call invert_matrix(sim(:, 1:n), test_inv, ok)
        if (.not. ok .or. inverse_error(test_inv, sim(:, 1:n)) > 1.0_dp) then
            info = cobyla_damaging_rounding
            return
        end if
        simi = test_inv
        fval(jdrop) = f
        if (size(constr) > 0) conmat(:, jdrop) = constr
        cval(jdrop) = cstrv
        call update_pole(cpen, conmat, cval, fval, sim, simi, info)
    end subroutine update_xfc

    integer function set_drop_tr(ximproved, d, delta, rho, sim, simi) result(jdrop)
        logical, intent(in) :: ximproved
        real(dp), intent(in) :: d(:)
        real(dp), intent(in) :: delta
        real(dp), intent(in) :: rho
        real(dp), intent(in) :: sim(:, :)
        real(dp), intent(in) :: simi(:, :)

        real(dp), allocatable :: distsq(:)
        real(dp), allocatable :: score(:)
        real(dp), allocatable :: simid(:)
        real(dp) :: base
        real(dp) :: weight
        integer :: j
        integer :: n

        n = size(d)
        allocate(distsq(n + 1), score(n + 1), simid(n))
        if (ximproved) then
            do j = 1, n
                distsq(j) = dot_product(sim(:, j) - d, sim(:, j) - d)
            end do
            distsq(n + 1) = dot_product(d, d)
        else
            do j = 1, n
                distsq(j) = dot_product(sim(:, j), sim(:, j))
            end do
            distsq(n + 1) = 0.0_dp
        end if
        simid = matmul(simi, d)
        base = max(rho, delta / 10.0_dp)**2
        do j = 1, n
            weight = max(1.0_dp, distsq(j) / base)
            score(j) = weight * abs(simid(j))
            if (ieee_is_nan(score(j))) score(j) = -1.0_dp
        end do
        weight = max(1.0_dp, distsq(n + 1) / base)
        score(n + 1) = weight * abs(1.0_dp - sum(simid))
        if (.not. ximproved) score(n + 1) = -1.0_dp
        if (ieee_is_nan(score(n + 1))) score(n + 1) = -1.0_dp

        jdrop = 0
        if (any(score > 0.0_dp)) jdrop = maxloc(score, dim=1)
        if (ximproved .and. jdrop == 0) jdrop = maxloc(distsq, dim=1)
    end function set_drop_tr

    subroutine geometry_step(jdrop, conmat, cpen, delbar, fval, simi, d)
        integer, intent(in) :: jdrop
        real(dp), intent(in) :: conmat(:, :)
        real(dp), intent(in) :: cpen
        real(dp), intent(in) :: delbar
        real(dp), intent(in) :: fval(:)
        real(dp), intent(in) :: simi(:, :)
        real(dp), intent(out) :: d(:)

        real(dp), allocatable :: a(:, :)
        real(dp), allocatable :: g(:)
        real(dp), allocatable :: model(:)
        real(dp) :: cvnd
        real(dp) :: cvpd
        real(dp) :: dn
        integer :: m
        integer :: n

        n = size(d)
        m = size(conmat, 1)
        allocate(a(n, m), g(n), model(m))
        d = simi(jdrop, :)
        dn = vecnorm(d)
        if (dn > 0.0_dp) d = delbar * d / dn
        call linear_models(conmat, fval, simi, g, a)
        if (m > 0) then
            model = conmat(:, n + 1) + matmul(transpose(a), d)
            cvpd = constraint_violation(model)
            model = conmat(:, n + 1) - matmul(transpose(a), d)
            cvnd = constraint_violation(model)
        else
            cvpd = 0.0_dp
            cvnd = 0.0_dp
        end if
        if (-dot_product(d, g) + cpen * cvnd < &
                dot_product(d, g) + cpen * cvpd) d = -d
    end subroutine geometry_step

    subroutine trstlp(a, b, delta, g, d)
        real(dp), intent(in) :: a(:, :)
        real(dp), intent(in) :: b(:)
        real(dp), intent(in) :: delta
        real(dp), intent(in) :: g(:)
        real(dp), intent(out) :: d(:)

        real(dp), allocatable :: aaug(:, :)
        real(dp), allocatable :: baug(:)
        real(dp), allocatable :: vmultc(:)
        real(dp), allocatable :: z(:, :)
        real(dp) :: maxval_col
        real(dp) :: modscal
        integer, allocatable :: iact(:)
        integer :: i
        integer :: m
        integer :: n
        integer :: nact

        n = size(g)
        m = size(b)
        allocate(aaug(n, m + 1), baug(m + 1), vmultc(m + 1), iact(m + 1), z(n, n))
        if (m > 0) then
            aaug(:, 1:m) = a
            baug(1:m) = b
        end if
        aaug(:, m + 1) = g
        baug(m + 1) = 0.0_dp

        do i = 1, m + 1
            maxval_col = maxval(abs(aaug(:, i)))
            if (maxval_col > 1.0e12_dp) then
                modscal = max(2.0_dp * realmin, 1.0_dp / maxval_col)
                aaug(:, i) = modscal * aaug(:, i)
                baug(i) = modscal * baug(i)
            end if
        end do

        vmultc = 0.0_dp
        iact = 0
        nact = 0
        d = 0.0_dp
        z = 0.0_dp
        call trstlp_sub(iact(1:m), nact, 1, aaug(:, 1:m), baug(1:m), delta, d, &
            vmultc(1:m), z)
        call trstlp_sub(iact, nact, 2, aaug, baug, delta, d, vmultc, z)
    end subroutine trstlp

    subroutine trstlp_sub(iact, nact, stage, a, b, delta, d, vmultc, z)
        integer, intent(inout) :: iact(:)
        integer, intent(inout) :: nact
        integer, intent(in) :: stage
        real(dp), intent(in) :: a(:, :)
        real(dp), intent(in) :: b(:)
        real(dp), intent(in) :: delta
        real(dp), intent(inout) :: d(:)
        real(dp), intent(inout) :: vmultc(:)
        real(dp), intent(inout) :: z(:, :)

        real(dp), allocatable :: active_a(:, :)
        real(dp), allocatable :: cvshift(:)
        real(dp), allocatable :: cvsabs(:)
        real(dp), allocatable :: dnew(:)
        real(dp), allocatable :: dold(:)
        real(dp), allocatable :: sdirn(:)
        real(dp), allocatable :: vmultd(:)
        real(dp), allocatable :: zdasav(:)
        real(dp), allocatable :: zdota(:)
        real(dp) :: cviol
        real(dp) :: dd
        real(dp) :: fm
        real(dp) :: frac
        real(dp) :: optnew
        real(dp) :: optold
        real(dp) :: sd
        real(dp) :: sqrtd
        real(dp) :: ss
        real(dp) :: step
        integer :: icon
        integer :: i
        integer :: iter
        integer :: j
        integer :: maxiter
        integer :: mcon
        integer :: n
        integer :: nactold
        integer :: nactsav
        integer :: nfail
        integer :: num_constraints
        integer :: oldiact
        real(dp) :: oldvmult
        logical :: have_positive

        mcon = size(b)
        n = size(d)
        allocate(cvshift(mcon), cvsabs(mcon), dnew(n), dold(n), sdirn(n), &
            vmultd(mcon), zdasav(n), zdota(n))
        cvshift = 0.0_dp
        cvsabs = 0.0_dp
        dnew = 0.0_dp
        dold = 0.0_dp
        sdirn = 0.0_dp
        vmultd = 0.0_dp
        zdasav = 0.0_dp
        zdota = 0.0_dp

        if (stage == 1) then
            do i = 1, mcon
                iact(i) = i
            end do
            nact = 0
            d = 0.0_dp
            if (mcon > 0) then
                cviol = max(0.0_dp, maxval(-b))
                vmultc = cviol + b
            else
                cviol = 0.0_dp
            end if
            z = 0.0_dp
            do i = 1, n
                z(i, i) = 1.0_dp
            end do
            if (mcon == 0 .or. cviol <= 0.0_dp) return
            icon = maxloc(-b, dim=1)
            num_constraints = mcon
            sdirn = 0.0_dp
        else
            if (dot_product(d, d) >= delta * delta) return
            iact(mcon) = mcon
            vmultc(mcon) = 0.0_dp
            num_constraints = mcon - 1
            icon = mcon
            if (num_constraints > 0) then
                cviol = 0.0_dp
                do i = 1, num_constraints
                    cviol = max(cviol, dot_product(d, a(:, i)) - b(i))
                end do
            else
                cviol = 0.0_dp
            end if
        end if

        do i = 1, nact
            zdota(i) = dot_product(z(:, i), a(:, iact(i)))
        end do
        optold = realmax
        nactold = nact
        nfail = 0
        maxiter = min(10000, 100 * max(num_constraints, n))

        do iter = 1, maxiter
            if (stage == 1) then
                optnew = cviol
            else
                optnew = dot_product(d, a(:, mcon))
            end if
            if (optnew < optold .or. nact > nactold) then
                nactold = nact
                nfail = 0
            else
                nfail = nfail + 1
            end if
            optold = min(optold, optnew)
            if (nfail == 3) exit

            if (icon > nact) then
                if (nact > 0) zdasav(1:nact) = zdota(1:nact)
                nactsav = nact
                call qradd_rdiag(a(:, iact(icon)), z, zdota, nact)
                if (nact == nactsav + 1) then
                    if (nact /= icon) then
                        oldvmult = vmultc(nact)
                        vmultc(nact) = 0.0_dp
                        vmultc(icon) = oldvmult
                        oldiact = iact(icon)
                        iact(icon) = iact(nact)
                        iact(nact) = oldiact
                    else
                        vmultc(nact) = 0.0_dp
                    end if
                else
                    if (nact <= 0) exit
                    allocate(active_a(n, nact))
                    do j = 1, nact
                        active_a(:, j) = a(:, iact(j))
                    end do
                    call lsqr_active(active_a, a(:, iact(icon)), z(:, 1:nact), &
                        zdasav(1:nact), vmultd(1:nact))
                    deallocate(active_a)
                    have_positive = .false.
                    do i = 1, nact
                        if (vmultd(i) > 0.0_dp .and. iact(i) <= num_constraints) then
                            have_positive = .true.
                            exit
                        end if
                    end do
                    if (.not. have_positive) exit
                    frac = realmax
                    do i = 1, nact
                        if (vmultd(i) > 0.0_dp .and. iact(i) <= num_constraints) &
                            frac = min(frac, vmultc(i) / vmultd(i))
                    end do
                    vmultc(1:nact) = max(0.0_dp, vmultc(1:nact) - frac * vmultd(1:nact))
                    if (ieee_is_nan(zdota(nact)) .or. abs(zdota(nact)) <= eps**2) exit
                    vmultc(icon) = 0.0_dp
                    vmultc(nact) = frac
                    oldiact = iact(icon)
                    iact(icon) = iact(nact)
                    iact(nact) = oldiact
                end if

                if (stage == 2 .and. iact(nact) /= mcon) then
                    if (nact <= 1) exit
                    allocate(active_a(n, nact))
                    do j = 1, nact
                        active_a(:, j) = a(:, iact(j))
                    end do
                    call qrexc_rdiag(active_a, z, zdota(1:nact), nact - 1)
                    deallocate(active_a)
                    oldiact = iact(nact - 1)
                    iact(nact - 1) = iact(nact)
                    iact(nact) = oldiact
                    oldvmult = vmultc(nact - 1)
                    vmultc(nact - 1) = vmultc(nact)
                    vmultc(nact) = oldvmult
                end if
                if (nact <= 0) exit
                if (ieee_is_nan(zdota(nact)) .or. abs(zdota(nact)) <= eps**2) exit
                if (stage == 1) then
                    sdirn = sdirn - &
                        ((dot_product(sdirn, a(:, iact(nact))) + 1.0_dp) / zdota(nact)) * &
                        z(:, nact)
                else
                    sdirn = -z(:, nact) / zdota(nact)
                end if
            else
                if (nact <= 0) exit
                allocate(active_a(n, nact))
                do j = 1, nact
                    active_a(:, j) = a(:, iact(j))
                end do
                call qrexc_rdiag(active_a, z, zdota(1:nact), icon)
                deallocate(active_a)
                oldiact = iact(icon)
                oldvmult = vmultc(icon)
                do i = icon, nact - 1
                    iact(i) = iact(i + 1)
                    vmultc(i) = vmultc(i + 1)
                end do
                iact(nact) = oldiact
                vmultc(nact) = oldvmult
                nact = nact - 1
                if (stage == 2 .and. nact <= 0) exit
                if (nact > 0) then
                    if (ieee_is_nan(zdota(nact)) .or. abs(zdota(nact)) <= eps**2) exit
                end if
                if (stage == 1) then
                    sdirn = sdirn - dot_product(sdirn, z(:, nact + 1)) * z(:, nact + 1)
                else
                    sdirn = -z(:, nact) / zdota(nact)
                end if
            end if

            dd = delta * delta - dot_product(d, d)
            ss = dot_product(sdirn, sdirn)
            sd = dot_product(sdirn, d)
            if (dd <= 0.0_dp .or. ss <= eps * delta * delta .or. ieee_is_nan(sd)) exit
            sqrtd = max(sqrt(max(0.0_dp, ss * dd + sd * sd)), abs(sd), &
                sqrt(max(0.0_dp, ss * dd)))
            if (sd > 0.0_dp) then
                step = dd / (sqrtd + sd)
            else
                step = (sqrtd - sd) / ss
            end if
            if (step <= 0.0_dp .or. .not. ieee_is_finite(step)) exit
            if (stage == 1) then
                if (is_minor(cviol, step)) exit
                step = min(step, cviol)
            end if

            dnew = d + step * sdirn
            if (stage == 1) then
                cviol = 0.0_dp
                do i = 1, nact
                    cviol = max(cviol, dot_product(dnew, a(:, iact(i))) - b(iact(i)))
                end do
            end if

            if (nact > 0) then
                allocate(active_a(n, nact))
                do j = 1, nact
                    active_a(:, j) = a(:, iact(j))
                end do
                call lsqr_active(active_a, dnew, z(:, 1:nact), zdota(1:nact), vmultd(1:nact))
                vmultd(1:nact) = -vmultd(1:nact)
                deallocate(active_a)
                if (stage == 2) vmultd(nact) = max(0.0_dp, vmultd(nact))
            end if

            do i = 1, mcon
                cvshift(i) = cviol - (dot_product(dnew, a(:, iact(i))) - b(iact(i)))
                cvsabs(i) = dot_product(abs(dnew), abs(a(:, iact(i)))) + &
                    abs(b(iact(i))) + cviol
                if (is_minor(cvshift(i), cvsabs(i))) cvshift(i) = 0.0_dp
            end do
            if (nact < mcon) vmultd(nact + 1:mcon) = cvshift(nact + 1:mcon)

            frac = 1.0_dp
            icon = 0
            do i = 1, mcon
                if (vmultd(i) < 0.0_dp) then
                    fm = vmultc(i) / (vmultc(i) - vmultd(i))
                    if (fm < frac) then
                        frac = fm
                        icon = i
                    end if
                end if
            end do
            dold = d
            d = (1.0_dp - frac) * d + frac * dnew
            vmultc = max(0.0_dp, (1.0_dp - frac) * vmultc + frac * vmultd)
            if (.not. all(ieee_is_finite(d)) .or. .not. all(ieee_is_finite(vmultc))) then
                d = dold
                exit
            end if
            if (stage == 1) then
                cviol = 0.0_dp
                do i = 1, mcon
                    cviol = max(cviol, dot_product(d, a(:, i)) - b(i))
                end do
            end if
            if (icon <= 0 .or. icon > mcon) exit
        end do
    end subroutine trstlp_sub

    subroutine qradd_rdiag(c, q, rdiag, nact)
        real(dp), intent(in) :: c(:)
        real(dp), intent(inout) :: q(:, :)
        real(dp), intent(inout) :: rdiag(:)
        integer, intent(inout) :: nact

        real(dp), allocatable :: cq(:)
        real(dp), allocatable :: cqa(:)
        real(dp), allocatable :: qk(:)
        real(dp), allocatable :: qkp1(:)
        real(dp) :: cc
        real(dp) :: ss
        integer :: k
        integer :: m

        m = size(q, 2)
        allocate(cq(m), cqa(m), qk(size(c)), qkp1(size(c)))
        cq = matmul(transpose(q), c)
        cqa = matmul(transpose(abs(q)), abs(c))
        do k = 1, m
            if (is_minor(cq(k), cqa(k))) cq(k) = 0.0_dp
        end do
        do k = m - 1, nact + 1, -1
            if (abs(cq(k + 1)) > 0.0_dp) then
                call plane_rotation(cq(k), cq(k + 1), cc, ss)
                qk = q(:, k)
                qkp1 = q(:, k + 1)
                q(:, k) = cc * qk + ss * qkp1
                q(:, k + 1) = -ss * qk + cc * qkp1
                cq(k) = stable_hypot(cq(k), cq(k + 1))
                cq(k + 1) = 0.0_dp
            end if
        end do
        if (nact < m) then
            if (abs(cq(nact + 1)) > eps**2 .and. &
                    .not. is_minor(cq(nact + 1), cqa(nact + 1))) nact = nact + 1
        end if
        if (nact >= 1 .and. nact <= size(rdiag)) rdiag(nact) = cq(nact)
    end subroutine qradd_rdiag

    subroutine qrexc_rdiag(a, q, rdiag, istart)
        real(dp), intent(in) :: a(:, :)
        real(dp), intent(inout) :: q(:, :)
        real(dp), intent(inout) :: rdiag(:)
        integer, intent(in) :: istart

        real(dp), allocatable :: qk(:)
        real(dp), allocatable :: qkp1(:)
        real(dp) :: cc
        real(dp) :: ss
        integer :: k
        integer :: n

        n = size(a, 2)
        if (istart < 1 .or. istart > n) return
        allocate(qk(size(a, 1)), qkp1(size(a, 1)))
        do k = istart, n - 1
            call plane_rotation(rdiag(k + 1), dot_product(q(:, k), a(:, k + 1)), cc, ss)
            qk = q(:, k)
            qkp1 = q(:, k + 1)
            q(:, k) = cc * qkp1 + ss * qk
            q(:, k + 1) = -ss * qkp1 + cc * qk
        end do
        do k = istart, n - 1
            rdiag(k) = dot_product(q(:, k), a(:, k + 1))
        end do
        rdiag(n) = dot_product(q(:, n), a(:, istart))
    end subroutine qrexc_rdiag

    subroutine lsqr_active(a, b, q, rdiag, x)
        real(dp), intent(in) :: a(:, :)
        real(dp), intent(in) :: b(:)
        real(dp), intent(in) :: q(:, :)
        real(dp), intent(in) :: rdiag(:)
        real(dp), intent(out) :: x(:)

        real(dp), allocatable :: y(:)
        real(dp) :: yq
        real(dp) :: yqa
        integer :: i
        integer :: rank

        x = 0.0_dp
        if (size(x) <= 0) return
        allocate(y(size(b)))
        y = b
        rank = min(size(a, 1), size(a, 2))
        do i = rank, 1, -1
            yq = dot_product(y, q(:, i))
            yqa = dot_product(abs(y), abs(q(:, i)))
            if (is_minor(yq, yqa) .or. abs(rdiag(i)) <= realmin) then
                x(i) = 0.0_dp
            else
                x(i) = yq / rdiag(i)
                y = y - x(i) * a(:, i)
            end if
        end do
    end subroutine lsqr_active

    pure logical function is_minor(x, ref) result(minor)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: ref
        real(dp) :: refa
        real(dp) :: refb

        refa = abs(ref) + 0.1_dp * abs(x)
        refb = abs(ref) + 0.2_dp * abs(x)
        minor = abs(ref) >= refa .or. refa >= refb
    end function is_minor

    subroutine plane_rotation(x1, x2, c, s)
        real(dp), intent(in) :: x1
        real(dp), intent(in) :: x2
        real(dp), intent(out) :: c
        real(dp), intent(out) :: s
        real(dp) :: r

        if (ieee_is_nan(x1) .or. ieee_is_nan(x2)) then
            c = 1.0_dp
            s = 0.0_dp
        else
            r = stable_hypot(x1, x2)
            if (r <= 0.0_dp .or. .not. ieee_is_finite(r)) then
                if (abs(x1) >= abs(x2) .and. abs(x1) > 0.0_dp) then
                    c = sign(1.0_dp, x1)
                    s = 0.0_dp
                else if (abs(x2) > 0.0_dp) then
                    c = 0.0_dp
                    s = sign(1.0_dp, x2)
                else
                    c = 1.0_dp
                    s = 0.0_dp
                end if
            else
                c = x1 / r
                s = x2 / r
            end if
        end if
    end subroutine plane_rotation

    pure real(dp) function stable_hypot(x1, x2) result(r)
        real(dp), intent(in) :: x1
        real(dp), intent(in) :: x2
        real(dp) :: hi
        real(dp) :: lo

        hi = max(abs(x1), abs(x2))
        lo = min(abs(x1), abs(x2))
        if (hi > 0.0_dp) then
            r = hi * sqrt(1.0_dp + (lo / hi)**2)
        else
            r = 0.0_dp
        end if
    end function stable_hypot

    pure real(dp) function vecnorm(x) result(r)
        real(dp), intent(in) :: x(:)
        r = sqrt(max(0.0_dp, dot_product(x, x)))
    end function vecnorm

    pure real(dp) function reduction_ratio(ared, pred, rshrink) result(ratio)
        real(dp), intent(in) :: ared
        real(dp), intent(in) :: pred
        real(dp), intent(in) :: rshrink

        if (ieee_is_nan(ared)) then
            ratio = -realmax
        else if (ieee_is_nan(pred) .or. pred <= 0.0_dp) then
            if (ared > 0.0_dp) then
                ratio = 0.5_dp * rshrink
            else
                ratio = -realmax
            end if
        else if (.not. ieee_is_finite(pred)) then
            if (ared > 0.0_dp .and. .not. ieee_is_finite(ared)) then
                ratio = 1.0_dp
            else if (ared < 0.0_dp .and. .not. ieee_is_finite(ared)) then
                ratio = -realmax
            else
                ratio = ared / pred
            end if
        else
            ratio = ared / pred
        end if
    end function reduction_ratio

    pure real(dp) function trust_radius(delta_in, dnorm, eta1, eta2, gamma1, gamma2, &
            ratio) result(delta)
        real(dp), intent(in) :: delta_in
        real(dp), intent(in) :: dnorm
        real(dp), intent(in) :: eta1
        real(dp), intent(in) :: eta2
        real(dp), intent(in) :: gamma1
        real(dp), intent(in) :: gamma2
        real(dp), intent(in) :: ratio

        if (ratio <= eta1) then
            delta = gamma1 * dnorm
        else if (ratio <= eta2) then
            delta = max(gamma1 * delta_in, dnorm)
        else
            delta = max(gamma1 * delta_in, gamma2 * dnorm)
        end if
    end function trust_radius

    pure real(dp) function reduce_rho_value(rho, rhoend) result(newrho)
        real(dp), intent(in) :: rho
        real(dp), intent(in) :: rhoend
        real(dp) :: ratio

        ratio = rho / rhoend
        if (ratio > 250.0_dp) then
            newrho = 0.1_dp * rho
        else if (ratio <= 16.0_dp) then
            newrho = rhoend
        else
            newrho = sqrt(ratio) * rhoend
        end if
    end function reduce_rho_value

    subroutine invert_matrix(a, ainv, ok)
        real(dp), intent(in) :: a(:, :)
        real(dp), intent(out) :: ainv(:, :)
        logical, intent(out) :: ok

        real(dp), allocatable :: aug(:, :)
        real(dp), allocatable :: temp(:)
        real(dp) :: factor
        real(dp) :: pivot
        real(dp) :: scale
        integer :: i
        integer :: k
        integer :: n
        integer :: p

        n = size(a, 1)
        ok = .false.
        ainv = 0.0_dp
        if (size(a, 2) /= n .or. size(ainv, 1) /= n .or. size(ainv, 2) /= n) return
        allocate(aug(n, 2 * n), temp(2 * n))
        aug = 0.0_dp
        aug(:, 1:n) = a
        do i = 1, n
            aug(i, n + i) = 1.0_dp
        end do

        do k = 1, n
            p = k
            pivot = abs(aug(k, k))
            do i = k + 1, n
                if (abs(aug(i, k)) > pivot) then
                    p = i
                    pivot = abs(aug(i, k))
                end if
            end do
            scale = max(1.0_dp, maxval(abs(aug(k:n, k:n))))
            if (.not. ieee_is_finite(pivot) .or. pivot <= realmin * scale) return
            if (p /= k) then
                temp = aug(k, :)
                aug(k, :) = aug(p, :)
                aug(p, :) = temp
            end if
            factor = aug(k, k)
            aug(k, :) = aug(k, :) / factor
            do i = 1, n
                if (i == k) cycle
                factor = aug(i, k)
                aug(i, :) = aug(i, :) - factor * aug(k, :)
            end do
        end do
        ainv = aug(:, n + 1:2 * n)
        ok = all(ieee_is_finite(ainv))
    end subroutine invert_matrix

    pure real(dp) function inverse_error(ainv, a) result(err)
        real(dp), intent(in) :: ainv(:, :)
        real(dp), intent(in) :: a(:, :)
        real(dp), allocatable :: prod(:, :)
        integer :: i
        integer :: n

        n = size(a, 1)
        allocate(prod(n, n))
        prod = matmul(ainv, a)
        do i = 1, n
            prod(i, i) = prod(i, i) - 1.0_dp
        end do
        if (all(ieee_is_finite(prod))) then
            err = maxval(abs(prod))
        else
            err = realmax
        end if
    end function inverse_error

    pure elemental logical function exact_equal_real(a, b) result(equal)
        real(dp), intent(in) :: a
        real(dp), intent(in) :: b

        equal = (a <= b .and. a >= b)
    end function exact_equal_real

    subroutine swap_scalar(a, b)
        real(dp), intent(inout) :: a
        real(dp), intent(inout) :: b
        real(dp) :: t
        t = a
        a = b
        b = t
    end subroutine swap_scalar

    subroutine swap_columns(a, i, j)
        real(dp), intent(inout) :: a(:, :)
        integer, intent(in) :: i
        integer, intent(in) :: j
        real(dp), allocatable :: t(:)
        if (i == j) return
        allocate(t(size(a, 1)))
        t = a(:, i)
        a(:, i) = a(:, j)
        a(:, j) = t
    end subroutine swap_columns

end module scifort_optimize_cobyla_prima
