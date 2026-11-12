! SPDX-License-Identifier: MIT
! Compare fixed, operation-count, and CPU-timed Krylov refresh policies.
! The analytic nonlinear system is unchanged across modes. Callback CPU
! workloads are intentionally synthetic and are not physical timing estimates.
program benchmark_krylov_timed_refresh
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_optimize_root, only : root, vector_root_result
    implicit none

    integer, parameter :: n = 128, nmodes = 8, nprofiles = 2
    real(dp), parameter :: rates(2) = [0.9_dp, 2.0_dp]
    character(len=20), parameter :: modes(nmodes) = [character(len=20) :: &
        'fixed_initial', 'fixed_every', 'fixed_four', 'work_manual', &
        'work_measured_1x', 'work_measured_16x', 'cpu_timed', 'cpu_timed_off']
    character(len=20), parameter :: profiles(nprofiles) = [character(len=20) :: &
        'forward_expensive', 'transpose_expensive']
    real(dp) :: a(n), u(n,2), v(n,2), target(n), x0(n), target_value(n)
    real(dp) :: forward_seconds, transpose_seconds, base_seconds, rate
    real(dp) :: start_time, end_time, residual_norm, error, theta
    integer :: i, irate, iprofile, imode, interval, failures
    integer :: jv_calls, jtv_calls, base_calls, nfev_calls
    integer :: initial_jv, initial_jtv, initial_base, initial_nfev, initial_nit
    type(vector_root_result) :: fit

    print '(a)', 'profile,nonlinearity,mode,success,nit,nfev,jv_calls,jtv_calls,'// &
        'base_calls,forward_budget_seconds,transpose_budget_seconds,solve_cpu_seconds,'// &
        'scaled_residual,solution_error'
    do i = 1, n
        theta = real(i,dp)
        a(i) = 10.0_dp**(-5.0_dp*real(i-1,dp)/real(n-1,dp))
        u(i,1) = 2.0_dp*a(i)*(0.8_dp+0.22_dp*sin(0.53_dp*theta))
        u(i,2) = 2.0_dp*a(i)*(-0.6_dp+0.17_dp*cos(0.61_dp*theta))
        v(i,1) = (0.9_dp+0.30_dp*cos(0.29_dp*theta))/real(n,dp)
        v(i,2) = (0.7_dp+0.42_dp*sin(0.37_dp*theta))/real(n,dp)
        target(i) = 0.08_dp*cos(0.31_dp*theta)
        x0(i) = target(i)+0.95_dp+0.25_dp*sin(0.23_dp*theta)
    end do
    failures = 0
    do iprofile = 1, nprofiles
        ! The actual cost depends on the machine. These synthetic CPU-time
        ! budgets make the relative callback expense explicit and repeatable.
        forward_seconds = 0.0008_dp
        transpose_seconds = 0.0001_dp
        if (iprofile == 2) then
            forward_seconds = 0.00015_dp
            transpose_seconds = 0.0020_dp
        end if
        base_seconds = 0.0001_dp
        do irate = 1, size(rates)
            rate = rates(irate)
            call nonlinear_action(target, target_value)
            do imode = 1, nmodes
                interval = 12
                if (imode == 1) interval = 0
                if (imode == 2) interval = 1
                if (imode == 3) interval = 4
                nfev_calls = 0
                jv_calls = 0
                jtv_calls = 0
                base_calls = 0
                call cpu_time(start_time)
                select case (imode)
                case (1:3)
                    fit = solve_fixed()
                case (4)
                    fit = solve_work_manual()
                case (5:6)
                    fit = solve_work_measured(imode == 6)
                case (7:8)
                    fit = solve_timed(imode == 8)
                end select
                call cpu_time(end_time)
                residual_norm = scaled_norm(fit%fun)
                error = maxval(abs(fit%x-target))
                if (fit%nfev /= nfev_calls) error stop 'function evaluation count mismatch'
                if (.not. all(ieee_is_finite(fit%x))) error stop 'nonfinite solution'
                if (.not. ieee_is_finite(residual_norm)) error stop 'nonfinite residual'
                if (.not. fit%success .or. residual_norm > 1.0e-9_dp .or. &
                        error > 1.0e-7_dp) failures = failures+1
                if (imode == 1) then
                    initial_jv = jv_calls
                    initial_jtv = jtv_calls
                    initial_base = base_calls
                    initial_nfev = fit%nfev
                    initial_nit = fit%nit
                else if (imode == 8) then
                    ! A deliberately excessive minimum timing interval must
                    ! suppress timing-based refresh on this fast-progress case.
                    if (jv_calls /= initial_jv .or. jtv_calls /= initial_jtv .or. &
                            base_calls /= initial_base .or. fit%nfev /= initial_nfev .or. &
                            fit%nit /= initial_nit) then
                        error stop 'timed minimum-duration guard altered refresh scheduling'
                    end if
                end if
                print '(*(g0,:,","))', trim(profiles(iprofile)), rate, trim(modes(imode)), &
                    fit%success, fit%nit, fit%nfev, jv_calls, jtv_calls, base_calls, &
                    forward_seconds, transpose_seconds, end_time-start_time, residual_norm, error
            end do
        end do
    end do
    if (failures /= 0) then
        write (*,'(a,i0)') 'nonconverged configurations: ', failures
        error stop 1
    end if
contains
    function solve_fixed() result(res)
        type(vector_root_result) :: res
        res = root(residual,x0,method='krylov',inner_method='gmres', &
            inner_maxiter=60,maxiter=50,fatol=1.0e-10_dp,tol_norm=scaled_norm, &
            inner_m_diag_auto=.true.,inner_m_diag_update_interval=interval, &
            inner_m_diag_rank_two=.true.,inner_m_diag_rank_two_discover=.true., &
            inner_m_diag_rank_two_nonsymmetric_base_fun=base_diagonal, &
            inner_m_diag_rank_two_nonsymmetric_adaptive_transpose_cost=16.0_dp, &
            jacobian_product=jvp,jacobian_transpose_product=jtv)
    end function solve_fixed

    function solve_work_manual() result(res)
        type(vector_root_result) :: res
        res = root(residual,x0,method='krylov',inner_method='gmres', &
            inner_maxiter=60,maxiter=50,fatol=1.0e-10_dp,tol_norm=scaled_norm, &
            inner_m_diag_auto=.true.,inner_m_diag_update_interval=interval, &
            inner_m_diag_update_adaptive=.true.,inner_m_diag_update_stall_ratio=0.8_dp, &
            inner_m_diag_update_work_growth=1.2_dp,inner_m_diag_update_rebuild_cost=4, &
            inner_m_diag_rank_two=.true.,inner_m_diag_rank_two_discover=.true., &
            inner_m_diag_rank_two_nonsymmetric_base_fun=base_diagonal, &
            inner_m_diag_rank_two_nonsymmetric_adaptive_transpose_cost=16.0_dp, &
            jacobian_product=jvp,jacobian_transpose_product=jtv)
    end function solve_work_manual

    function solve_work_measured(expensive) result(res)
        logical, intent(in) :: expensive !! set relative transpose weight to sixteen
        real(dp) :: weight
        type(vector_root_result) :: res
        weight = 1.0_dp
        if (expensive) weight = 16.0_dp
        res = root(residual,x0,method='krylov',inner_method='gmres', &
            inner_maxiter=60,maxiter=50,fatol=1.0e-10_dp,tol_norm=scaled_norm, &
            inner_m_diag_auto=.true.,inner_m_diag_update_interval=interval, &
            inner_m_diag_update_adaptive=.true.,inner_m_diag_update_stall_ratio=0.8_dp, &
            inner_m_diag_update_work_growth=1.2_dp,inner_m_diag_update_measured_cost=.true., &
            inner_m_diag_update_transpose_cost=weight, &
            inner_m_diag_rank_two=.true.,inner_m_diag_rank_two_discover=.true., &
            inner_m_diag_rank_two_nonsymmetric_base_fun=base_diagonal, &
            inner_m_diag_rank_two_nonsymmetric_adaptive_transpose_cost=16.0_dp, &
            jacobian_product=jvp,jacobian_transpose_product=jtv)
    end function solve_work_measured

    function solve_timed(disable_early) result(res)
        logical, intent(in) :: disable_early !! force unreliable-clock cutoff
        type(vector_root_result) :: res
        real(dp) :: min_seconds
        min_seconds = 1.0e-4_dp
        if (disable_early) min_seconds = 100.0_dp
        res = root(residual,x0,method='krylov',inner_method='gmres', &
            inner_maxiter=60,maxiter=50,fatol=1.0e-10_dp,tol_norm=scaled_norm, &
            inner_m_diag_auto=.true.,inner_m_diag_update_interval=interval, &
            inner_m_diag_update_adaptive=.true.,inner_m_diag_update_stall_ratio=0.8_dp, &
            inner_m_diag_update_work_growth=1.2_dp,inner_m_diag_update_timed_cost=.true., &
            inner_m_diag_update_timed_min_seconds=min_seconds, &
            inner_m_diag_rank_two=.true.,inner_m_diag_rank_two_discover=.true., &
            inner_m_diag_rank_two_nonsymmetric_base_fun=base_diagonal, &
            inner_m_diag_rank_two_nonsymmetric_adaptive_transpose_cost=16.0_dp, &
            jacobian_product=jvp,jacobian_transpose_product=jtv)
    end function solve_timed

    subroutine nonlinear_action(x, value)
        real(dp), intent(in) :: x(:) !! current nonlinear state
        real(dp), intent(out) :: value(:) !! unshifted nonlinear residual
        real(dp) :: s1, s2
        s1 = dot_product(v(:,1),x)
        s2 = dot_product(v(:,2),x)
        value = a*(x+0.5_dp*rate*x**3) + &
            u(:,1)*sin(rate*s1)/rate + u(:,2)*(s2+0.25_dp*rate*s2**3)
    end subroutine nonlinear_action

    subroutine residual(x, fvec)
        real(dp), intent(in) :: x(:) !! nonlinear iterate
        real(dp), intent(out) :: fvec(:) !! residual at iterate
        nfev_calls = nfev_calls+1
        call nonlinear_action(x,fvec)
        fvec = fvec-target_value
    end subroutine residual

    subroutine base_diagonal(x, diagonal, ok)
        real(dp), intent(in) :: x(:) !! current iterate
        real(dp), intent(out) :: diagonal(:) !! uncoupled diagonal
        logical, intent(out) :: ok !! validity flag
        base_calls = base_calls+1
        call synthetic_work(base_seconds)
        diagonal = a*(1.0_dp+1.5_dp*rate*x*x)
        ok = all(ieee_is_finite(diagonal))
    end subroutine base_diagonal

    subroutine jvp(x, direction, product, ok)
        real(dp), intent(in) :: x(:) !! current nonlinear state
        real(dp), intent(in) :: direction(:) !! Jacobian product direction
        real(dp), intent(out) :: product(:) !! Jacobian times direction
        logical, intent(out) :: ok !! validity flag
        real(dp) :: s1, s2
        jv_calls = jv_calls+1
        call synthetic_work(forward_seconds)
        s1 = dot_product(v(:,1),x)
        s2 = dot_product(v(:,2),x)
        product = a*(1.0_dp+1.5_dp*rate*x*x)*direction + &
            u(:,1)*cos(rate*s1)*dot_product(v(:,1),direction) + &
            u(:,2)*(1.0_dp+0.75_dp*rate*s2*s2)*dot_product(v(:,2),direction)
        ok = all(ieee_is_finite(product))
    end subroutine jvp

    subroutine jtv(x, direction, product, ok)
        real(dp), intent(in) :: x(:) !! current nonlinear state
        real(dp), intent(in) :: direction(:) !! transpose product direction
        real(dp), intent(out) :: product(:) !! transpose Jacobian times direction
        logical, intent(out) :: ok !! validity flag
        real(dp) :: s1, s2
        jtv_calls = jtv_calls+1
        call synthetic_work(transpose_seconds)
        s1 = dot_product(v(:,1),x)
        s2 = dot_product(v(:,2),x)
        product = a*(1.0_dp+1.5_dp*rate*x*x)*direction + &
            v(:,1)*cos(rate*s1)*dot_product(u(:,1),direction) + &
            v(:,2)*(1.0_dp+0.75_dp*rate*s2*s2)*dot_product(u(:,2),direction)
        ok = all(ieee_is_finite(product))
    end subroutine jtv

    subroutine synthetic_work(seconds)
        real(dp), intent(in) :: seconds !! desired artificial callback CPU expense
        real(dp) :: begin, now
        integer :: k
        call cpu_time(begin)
        if (begin < 0.0_dp) return
        do k = 1, 2000000
            call cpu_time(now)
            if (now < begin .or. now-begin >= seconds) exit
        end do
    end subroutine synthetic_work

    function scaled_norm(fvec) result(value)
        real(dp), intent(in) :: fvec(:) !! residual to scale by nominal diagonal
        real(dp) :: value !! max row-scaled absolute residual
        value = maxval(abs(fvec)/a)
    end function scaled_norm
end program benchmark_krylov_timed_refresh
