! SPDX-License-Identifier: MIT
! Progressive verified Rademacher diagonal probing, exact fallback and API checks.
program test_optimize_root_krylov_progressive_diagonal
    use scifort_kinds, only : dp
    use scifort_optimize_root, only : root, newton_krylov, vector_root_result, &
        vector_root_linear_operator
    implicit none

    integer, parameter :: n = 64
    real(dp) :: diagonal(n), x0(n), observed(n), saved(n)
    real(dp) :: distant
    logical :: single_row_coupling, two_row_coupling
    integer :: failures, f_calls, jv_calls, solver_calls, fail_after, i
    integer :: coupled_row, coupled_col
    type(vector_root_result) :: res

    failures = 0
    single_row_coupling = .false.
    two_row_coupling = .false.
    coupled_row = 1
    coupled_col = 43
    x0 = 0.25_dp
    do i = 1, n
        diagonal(i) = 1.0_dp + 0.01_dp * real(i, dp)
    end do

    ! At eight probes the two independent groups disagree. With a
    ! sixteen-probe limit the estimates pass both consistency and a
    ! coordinate check, avoiding 64 exact coordinate products.
    distant = 0.3_dp
    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_probes=8, inner_m_diag_probe_scheme='rademacher', &
        inner_m_diag_probe_seed=17, inner_m_diag_verify_samples=1, &
        inner_m_diag_verify_rtol=0.25_dp, inner_m_diag_probe_max=16, &
        jacobian_product=jvp, inner_solver=inspect)
    call check(jv_calls == 17, 'progressive estimate needs sixteen signs and one spot check')
    call check(solver_calls == 1 .and. res%nfev == 1, 'analytic estimate evaluation accounting')
    call check(maxval(abs(observed * diagonal - 1.0_dp)) < 0.18_dp, &
        'expanded probes keep moderate local diagonal error')
    saved = observed

    call reset()
    res = newton_krylov(residual, x0, inner_m_diag_auto=.true., &
        inner_m_diag_probes=8, inner_m_diag_probe_scheme='rademacher', &
        inner_m_diag_probe_seed=17, inner_m_diag_verify_samples=1, &
        inner_m_diag_verify_rtol=0.25_dp, inner_m_diag_probe_max=16, &
        jacobian_product=jvp, inner_solver=inspect)
    call check(jv_calls == 17 .and. all(observed == saved), &
        'same adaptive result through newton_krylov')

    ! Without progressive growth, the pre-existing path falls back to
    ! exact n-product estimation. This intentionally tests default behavior.
    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_probes=8, inner_m_diag_probe_scheme='rademacher', &
        inner_m_diag_probe_seed=17, inner_m_diag_verify_samples=1, &
        inner_m_diag_verify_rtol=0.25_dp, &
        jacobian_product=jvp, inner_solver=inspect)
    call check(jv_calls == n + 8, 'unchanged fixed-probe fallback uses exact estimation')
    call check(maxval(abs(observed * diagonal - 1.0_dp)) < 2.0e-14_dp, &
        'unchanged exact fallback gives true diagonal')

    ! Persistent disagreement triggers exact fallback after the specified
    ! maximum randomized budget, never silently accepting a bad estimate.
    distant = 0.5_dp
    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_probes=8, inner_m_diag_probe_scheme='rademacher', &
        inner_m_diag_probe_seed=17, inner_m_diag_verify_samples=1, &
        inner_m_diag_verify_rtol=0.1_dp, inner_m_diag_probe_max=16, &
        jacobian_product=jvp, inner_solver=inspect)
    call check(jv_calls == 16 + n, 'maximum budget triggers exact fallback')
    call check(maxval(abs(observed * diagonal - 1.0_dp)) < 2.0e-14_dp, &
        'progressive fallback recovers exact inverse diagonal')

    ! Equal initial and maximum budgets must use the original implementation.
    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_probes=8, inner_m_diag_probe_scheme='rademacher', &
        inner_m_diag_probe_seed=17, inner_m_diag_verify_samples=1, &
        inner_m_diag_verify_rtol=0.1_dp, inner_m_diag_probe_max=8, &
        jacobian_product=jvp, inner_solver=inspect)
    call check(jv_calls == 8 + n, 'equal budget retains exact previous fixed path')

    distant = 0.0_dp
    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_probes=8, inner_m_diag_probe_scheme='rademacher', &
        inner_m_diag_probe_seed=17, inner_m_diag_verify_samples=2, &
        inner_m_diag_probe_max=16, inner_solver=inspect)
    call check(f_calls == 11 .and. res%nfev == 11, 'two-point difference accounting')
    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_probes=8, inner_m_diag_probe_scheme='rademacher', &
        inner_m_diag_probe_seed=17, inner_m_diag_verify_samples=2, &
        inner_m_diag_probe_max=16, jacobian_diff_scheme='3-point', &
        inner_solver=inspect)
    call check(f_calls == 21 .and. res%nfev == 21, 'centered difference accounting')

    distant = 0.3_dp
    call reset()
    fail_after = 12
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_probes=8, inner_m_diag_probe_scheme='rademacher', &
        inner_m_diag_probe_seed=17, inner_m_diag_verify_samples=1, &
        inner_m_diag_verify_rtol=0.25_dp, inner_m_diag_probe_max=16, &
        jacobian_product=jvp, inner_solver=inspect)
    call check(.not. res%success .and. solver_calls == 0 .and. jv_calls == 12, &
        'failure during expansion prevents inner-solver dispatch')
    fail_after = 0

    ! This fixed seeded sign pattern gives the first eight probes a
    ! biased but batch-consistent estimate for coordinate one. A direct
    ! coordinate check rejects it; sixteen probes pass both checks.
    ! Reusing the exact first-coordinate value saves a second Jv call.
    single_row_coupling = .true.
    distant = 0.6_dp
    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_probes=8, inner_m_diag_probe_scheme='rademacher', &
        inner_m_diag_probe_seed=17, inner_m_diag_verify_samples=1, &
        inner_m_diag_verify_rtol=0.2_dp, inner_m_diag_probe_max=16, &
        jacobian_product=jvp, inner_solver=inspect)
    call check(jv_calls == 18 .and. solver_calls == 1, &
        'reuse old spot while inspecting a new coordinate after growth')
    call check(abs(observed(1) * diagonal(1) - 1.0_dp / 0.85_dp) < 1.0e-12_dp, &
        'expanded probe estimate retains expected diagonal')
    call check(maxval(abs(observed(2:n) * diagonal(2:n) - 1.0_dp)) < 2.0e-14_dp, &
        'unaffected rows remain exact')
    saved = observed

    call reset()
    res = newton_krylov(residual, x0, inner_m_diag_auto=.true., &
        inner_m_diag_probes=8, inner_m_diag_probe_scheme='rademacher', &
        inner_m_diag_probe_seed=17, inner_m_diag_verify_samples=1, &
        inner_m_diag_verify_rtol=0.2_dp, inner_m_diag_probe_max=16, &
        jacobian_product=jvp, inner_solver=inspect)
    call check(jv_calls == 18 .and. all(observed == saved), &
        'rotating spot checks through newton_krylov match root')

    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_probes=8, inner_m_diag_probe_scheme='rademacher', &
        inner_m_diag_probe_seed=17, inner_m_diag_verify_samples=1, &
        inner_m_diag_verify_rtol=0.1_dp, inner_m_diag_probe_max=16, &
        jacobian_product=jvp, inner_solver=inspect)
    call check(jv_calls == 16 + n, &
        'reuse rejected spot during exact fallback after exhausted budget')
    call check(maxval(abs(observed * diagonal - 1.0_dp)) < 2.0e-14_dp, &
        'cached exact coordinate produces the full exact diagonal')

    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_probes=8, inner_m_diag_probe_scheme='rademacher', &
        inner_m_diag_probe_seed=17, inner_m_diag_verify_samples=1, &
        inner_m_diag_verify_rtol=0.2_dp, inner_m_diag_probe_max=16, &
        inner_solver=inspect)
    call check(f_calls == 19 .and. res%nfev == 19, &
        'two-point differences reuse old and inspect new spots')
    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_probes=8, inner_m_diag_probe_scheme='rademacher', &
        inner_m_diag_probe_seed=17, inner_m_diag_verify_samples=1, &
        inner_m_diag_verify_rtol=0.2_dp, inner_m_diag_probe_max=16, &
        jacobian_diff_scheme='3-point', inner_solver=inspect)
    call check(f_calls == 37 .and. res%nfev == 37, &
        'three-point differences reuse old and inspect new spots')
    ! Three selected coordinates are 1, 22, and 43. An off-diagonal
    ! coupling in row 22 causes the second check to fail at eight probes.
    ! Cached exact coordinates 1 and 22 must be reused at sixteen probes
    ! and on fallback; the third spot is evaluated only if needed.
    coupled_row = 22
    coupled_col = 4
    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_probes=8, inner_m_diag_probe_scheme='rademacher', &
        inner_m_diag_probe_seed=17, inner_m_diag_verify_samples=3, &
        inner_m_diag_verify_rtol=0.2_dp, inner_m_diag_probe_max=16, &
        jacobian_product=jvp, inner_solver=inspect)
    call check(jv_calls == 21, 'previous and rotated spot coordinates each sampled once')
    call check(abs(observed(22) * diagonal(22) - 1.0_dp / 0.85_dp) < 1.0e-12_dp, &
        'multiple cached spot checks preserve expanded diagonal')

    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_probes=8, inner_m_diag_probe_scheme='rademacher', &
        inner_m_diag_probe_seed=17, inner_m_diag_verify_samples=3, &
        inner_m_diag_verify_rtol=0.1_dp, inner_m_diag_probe_max=16, &
        jacobian_product=jvp, inner_solver=inspect)
    call check(jv_calls == n + 16, &
        'both cached spot coordinates reused in exact fallback')
    call check(maxval(abs(observed * diagonal - 1.0_dp)) < 2.0e-14_dp, &
        'multiple cached coordinates remain exact during fallback')

    ! With seed 17, row 33 coupled to column 43 has the same biased
    ! eight- and sixteen-probe means, but its first four probes differ.
    ! Stage 1 fails the batch check before sampling row 1. At stage 2 the
    ! previous fixed spot would accept, while the rotated spot at row 33
    ! discovers the 30% bias and forces exact fallback. Its exact product
    ! must be reused during that fallback, not evaluated twice.
    coupled_row = 33
    coupled_col = 43
    distant = 0.6_dp
    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_probes=8, inner_m_diag_probe_scheme='rademacher', &
        inner_m_diag_probe_seed=17, inner_m_diag_verify_samples=1, &
        inner_m_diag_verify_rtol=0.2_dp, inner_m_diag_probe_max=16, &
        jacobian_product=jvp, inner_solver=inspect)
    call check(jv_calls == n + 16 .and. solver_calls == 1, &
        'rotated second-stage check detects previously hidden biased row')
    call check(maxval(abs(observed * diagonal - 1.0_dp)) < 2.0e-14_dp, &
        'rotated check triggers exact, correctly cached diagonal')

    ! The 2-point residual callback must count both stages and the n exact
    ! products, but must not count the already cached coordinate twice.
    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_probes=8, inner_m_diag_probe_scheme='rademacher', &
        inner_m_diag_probe_seed=17, inner_m_diag_verify_samples=1, &
        inner_m_diag_verify_rtol=0.2_dp, inner_m_diag_probe_max=16, &
        inner_solver=inspect)
    call check(f_calls == n + 17 .and. res%nfev == n + 17, &
        'rotated fallback has no duplicate two-point residual evaluation')

    ! At a 32-probe cap, the same row is still biased at sixteen but
    ! improves enough after another expansion to be safely approximated.
    ! Stage 3 inspects a third location (row 17) using a quarter-stratum
    ! offset; row 33's cached exact value is rechecked without another Jv.
    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_probes=8, inner_m_diag_probe_scheme='rademacher', &
        inner_m_diag_probe_seed=17, inner_m_diag_verify_samples=1, &
        inner_m_diag_verify_rtol=0.2_dp, inner_m_diag_probe_max=32, &
        jacobian_product=jvp, inner_solver=inspect)
    call check(jv_calls == 34 .and. solver_calls == 1, &
        'third stage checks a new coordinate without repeating cached checks')
    call check(abs(observed(33) * diagonal(33) - 1.0_dp / 0.85_dp) < 1.0e-12_dp, &
        'third-stage approximate diagonal preserves seeded-probe result')

    call reset()
    fail_after = 17
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_probes=8, inner_m_diag_probe_scheme='rademacher', &
        inner_m_diag_probe_seed=17, inner_m_diag_verify_samples=1, &
        inner_m_diag_verify_rtol=0.2_dp, inner_m_diag_probe_max=16, &
        jacobian_product=jvp, inner_solver=inspect)
    call check(.not. res%success .and. jv_calls == 17 .and. solver_calls == 0, &
        'failure in a rotated spot check aborts before inner solver')
    fail_after = 0

    single_row_coupling = .false.
    distant = 0.0_dp
    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_probes=8, inner_m_diag_probe_scheme='rademacher', &
        inner_m_diag_probe_seed=17, inner_m_diag_verify_samples=n, &
        inner_m_diag_probe_max=16, jacobian_product=jvp, inner_solver=inspect)
    call check(jv_calls == n + 8, 'full-coordinate spot sampling visits each row once')
    call check(maxval(abs(observed * diagonal - 1.0_dp)) < 2.0e-14_dp, &
        'full-coordinate progressive verification preserves exact diagonal')
    distant = 0.3_dp

    ! With more than half the coordinates requested as spot checks, the
    ! old integer-rounded stratum offset was zero at the second stage.
    ! It revisited its first set of sites and missed row five. This seed
    ! makes row five 30% too large in all three (4, 8, 16) probe averages.
    ! Row one fails verification at eight probes but passes at sixteen;
    ! the improved coverage must uncover row five and fall back to exact.
    two_row_coupling = .true.
    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_probes=8, inner_m_diag_probe_scheme='rademacher', &
        inner_m_diag_probe_seed=17, inner_m_diag_verify_samples=50, &
        inner_m_diag_verify_rtol=0.2_dp, inner_m_diag_probe_max=16, &
        jacobian_product=jvp, inner_solver=inspect)
    call check(jv_calls == n + 16 .and. solver_calls == 1, &
        'high-coverage stages inspect unvisited rows and reuse exact fallback spots')
    call check(maxval(abs(observed * diagonal - 1.0_dp)) < 2.0e-14_dp, &
        'previously hidden coupled row recovers its exact diagonal')
    ! The same coverage/fallback behavior holds when Jv is approximated
    ! through the residual, with exact function-evaluation accounting.
    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_probes=8, inner_m_diag_probe_scheme='rademacher', &
        inner_m_diag_probe_seed=17, inner_m_diag_verify_samples=50, &
        inner_m_diag_verify_rtol=0.2_dp, inner_m_diag_probe_max=16, &
        inner_solver=inspect)
    call check(f_calls == n + 17 .and. res%nfev == f_calls, &
        'high-coverage fallback two-point Jv count')
    call check(maxval(abs(observed * diagonal - 1.0_dp)) < 1.0e-7_dp, &
        'high-coverage two-point fallback recovers diagonal')
    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_probes=8, inner_m_diag_probe_scheme='rademacher', &
        inner_m_diag_probe_seed=17, inner_m_diag_verify_samples=50, &
        inner_m_diag_verify_rtol=0.2_dp, inner_m_diag_probe_max=16, &
        jacobian_diff_scheme='3-point', inner_solver=inspect)
    call check(f_calls == 2 * (n + 16) + 1 .and. res%nfev == f_calls, &
        'high-coverage fallback centered Jv count')
    call check(maxval(abs(observed * diagonal - 1.0_dp)) < 1.0e-7_dp, &
        'high-coverage centered fallback recovers diagonal')
    two_row_coupling = .false.

    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_probe_max=16)
    call rejected(res, 'maximum requires automatic preconditioner')
    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_probe_scheme='rademacher', inner_m_diag_probes=8, &
        inner_m_diag_probe_max=16)
    call rejected(res, 'maximum requires verification samples')
    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_probes=8, inner_m_diag_verify_samples=1, &
        inner_m_diag_probe_max=16)
    call rejected(res, 'maximum requires randomized scheme')
    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_probes=8, inner_m_diag_probe_scheme='rademacher', &
        inner_m_diag_verify_samples=1, inner_m_diag_probe_max=7)
    call rejected(res, 'maximum cannot be smaller than initial budget')
    call reset()
    res = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
        inner_m_diag_probes=8, inner_m_diag_probe_scheme='rademacher', &
        inner_m_diag_verify_samples=1, inner_m_diag_probe_max=n+1)
    call rejected(res, 'maximum cannot exceed problem dimension')
    call reset()
    res = root(residual, x0, method='hybr', inner_m_diag_probe_max=16)
    call rejected(res, 'non-Krylov method rejects adaptive probe control')

    if (failures /= 0) then
        print '(a,i0)', 'FAIL progressive diagonal probing: ', failures
        error stop 1
    end if
    print '(a)', 'PASS progressive verified diagonal probing'
contains
    subroutine reset()
        f_calls = 0
        jv_calls = 0
        solver_calls = 0
        fail_after = 0
    end subroutine reset

    subroutine check(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description
        if (.not. condition) then
            failures = failures + 1
            print '(a,a)', 'FAIL: ', description
        end if
    end subroutine check

    subroutine rejected(result, description)
        type(vector_root_result), intent(in) :: result
        character(len=*), intent(in) :: description
        call check(.not. result%success .and. result%nfev == 0 .and. f_calls == 0, description)
    end subroutine rejected

    subroutine apply_jacobian(v, product)
        real(dp), intent(in) :: v(:)
        real(dp), intent(out) :: product(:)
        integer :: k
        product = diagonal * v
        if (two_row_coupling) then
            product(1) = product(1) + 0.6_dp * diagonal(1) * v(43)
            product(5) = product(5) + 0.6_dp * diagonal(5) * v(28)
            return
        end if
        if (single_row_coupling) then
            product(coupled_row) = product(coupled_row) + &
                distant * diagonal(coupled_row) * v(coupled_col)
            return
        end if
        do k = 1, n - 8
            product(k) = product(k) + distant * diagonal(k) * v(k + 8)
        end do
    end subroutine apply_jacobian

    subroutine residual(x, f)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: f(:)
        f_calls = f_calls + 1
        call apply_jacobian(x, f)
    end subroutine residual

    subroutine jvp(x, vec, product, ok)
        real(dp), intent(in) :: x(:), vec(:)
        real(dp), intent(out) :: product(:)
        logical, intent(out) :: ok
        jv_calls = jv_calls + 1
        ok = fail_after == 0 .or. jv_calls < fail_after
        product = 0.0_dp
        if (ok) call apply_jacobian(vec, product)
    end subroutine jvp

    subroutine inspect(matvec, preconditioner, rhs, rtol, maxiter, solution, ok)
        procedure(vector_root_linear_operator) :: matvec, preconditioner
        real(dp), intent(in) :: rhs(:), rtol
        integer, intent(in) :: maxiter
        real(dp), intent(out) :: solution(:)
        logical, intent(out) :: ok
        logical :: valid
        solver_calls = solver_calls + 1
        call preconditioner(spread(1.0_dp, 1, size(rhs)), observed, valid)
        call check(valid, 'preconditioner output finite')
        solution = 0.0_dp
        ok = .false.
    end subroutine inspect
end program test_optimize_root_krylov_progressive_diagonal
