! SPDX-License-Identifier: MIT
! Compare fixed, manual-cost and measured-cost Krylov refresh policies.
! Nonlinear rank-two systems have a changing diagonal and coupled Jacobian.
! No single schedule is presumed to minimize computational work.
program benchmark_krylov_work_refresh
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_optimize_root, only : root, vector_root_result
    implicit none

    integer, parameter :: nrates = 3, nintervals = 5, npolicies = 11
    integer, parameter :: intervals(nintervals) = [0, 1, 2, 4, 8]
    real(dp), parameter :: rates(nrates) = [0.3_dp, 0.9_dp, 2.0_dp]
    character(len=16), parameter :: modes(npolicies) = [character(len=16) :: &
        'fixed0', 'fixed1', 'fixed2', 'fixed4', 'fixed8', &
        'stall', 'work_cheap', 'work_cost4', 'work_growth2', &
        'measured_1x', 'measured_16x']
    real(dp), allocatable :: a(:), u(:,:), v(:,:), target(:), x0(:), target_value(:)
    integer :: n, i, isize, iscale, irate, ipolicy, interval
    logical :: adaptive, work_aware, measured
    real(dp) :: stall_ratio, work_growth, transpose_weight
    integer :: rebuild_cost
    integer :: nfev_calls, jv_calls, jtv_calls, base_calls, iterations, failures
    real(dp) :: log_scale, rate, theta, start_time, end_time, error, residual_norm
    type(vector_root_result) :: fit

    print '(a)', 'n,log10_scale,nonlinearity,mode,refresh_interval,stall_ratio,work_growth,rebuild_cost,'// &
        'measured_cost,transpose_weight,'// &
        'success,nit,nfev,'// &
        'base_calls,jv_calls,jtv_calls,work_1x,work_4x,work_16x,'// &
        'cpu_seconds,scaled_residual,solution_error'
    failures = 0
    do isize = 1, 2
        n = 32
        if (isize == 2) n = 128
        allocate(a(n), u(n,2), v(n,2), target(n), x0(n), target_value(n))
        do iscale = 1, 2
            log_scale = 2.0_dp
            if (iscale == 2) log_scale = 5.0_dp
            do i = 1, n
                theta = real(i,dp)
                a(i) = 10.0_dp**(-log_scale*real(i-1,dp)/real(n-1,dp))
                u(i,1) = 2.0_dp*a(i)*(0.8_dp+0.22_dp*sin(0.53_dp*theta))
                u(i,2) = 2.0_dp*a(i)*(-0.6_dp+0.17_dp*cos(0.61_dp*theta))
                v(i,1) = (0.9_dp+0.30_dp*cos(0.29_dp*theta))/real(n,dp)
                v(i,2) = (0.7_dp+0.42_dp*sin(0.37_dp*theta))/real(n,dp)
                target(i) = 0.08_dp*cos(0.31_dp*theta)
                x0(i) = target(i)+0.95_dp+0.25_dp*sin(0.23_dp*theta)
            end do
            do irate = 1, nrates
                rate = rates(irate)
                call nonlinear_action(target, target_value)
                call verify_products()
                do ipolicy = 1, npolicies
                    adaptive = ipolicy > nintervals
                    work_aware = ipolicy >= 7
                    measured = ipolicy >= 10
                    transpose_weight = 1.0_dp
                    if (ipolicy == 11) transpose_weight = 16.0_dp
                    stall_ratio = 0.8_dp
                    work_growth = 0.0_dp
                    rebuild_cost = 0
                    if (work_aware) then
                        work_growth = 1.2_dp
                        if (ipolicy == 8) rebuild_cost = 4
                        if (ipolicy == 9) work_growth = 2.0_dp
                    end if
                    if (adaptive) then
                        interval = 12
                    else
                        interval = intervals(ipolicy)
                    end if
                    nfev_calls = 0
                    jv_calls = 0
                    jtv_calls = 0
                    base_calls = 0
                    call cpu_time(start_time)
                    if (measured) then
                        fit = root(residual, x0, method='krylov', &
                        inner_method='gmres', inner_maxiter=60, &
                        maxiter=50, fatol=1.0e-10_dp, tol_norm=scaled_norm, &
                        inner_m_diag_auto=.true., &
                        inner_m_diag_update_interval=interval, &
                        inner_m_diag_update_adaptive=.true., &
                        inner_m_diag_update_stall_ratio=0.8_dp, &
                        inner_m_diag_update_work_growth=work_growth, &
                        inner_m_diag_update_measured_cost=.true., &
                        inner_m_diag_update_transpose_cost=transpose_weight, &
                        inner_m_diag_rank_two=.true., &
                        inner_m_diag_rank_two_discover=.true., &
                        inner_m_diag_rank_two_nonsymmetric_base_fun=base_diagonal, &
                        inner_m_diag_rank_two_nonsymmetric_adaptive_transpose_cost=16.0_dp, &
                        jacobian_product=jvp, jacobian_transpose_product=jtv)
                    else if (work_aware) then
                        fit = root(residual, x0, method='krylov', &
                        inner_method='gmres', inner_maxiter=60, &
                        maxiter=50, fatol=1.0e-10_dp, tol_norm=scaled_norm, &
                        inner_m_diag_auto=.true., &
                        inner_m_diag_update_interval=interval, &
                        inner_m_diag_update_adaptive=.true., &
                        inner_m_diag_update_stall_ratio=0.8_dp, &
                        inner_m_diag_update_work_growth=work_growth, &
                        inner_m_diag_update_rebuild_cost=rebuild_cost, &
                        inner_m_diag_rank_two=.true., &
                        inner_m_diag_rank_two_discover=.true., &
                        inner_m_diag_rank_two_nonsymmetric_base_fun=base_diagonal, &
                        inner_m_diag_rank_two_nonsymmetric_adaptive_transpose_cost=16.0_dp, &
                        jacobian_product=jvp, jacobian_transpose_product=jtv)
                    else if (adaptive) then
                        fit = root(residual, x0, method='krylov', &
                        inner_method='gmres', inner_maxiter=60, &
                        maxiter=50, fatol=1.0e-10_dp, tol_norm=scaled_norm, &
                        inner_m_diag_auto=.true., &
                        inner_m_diag_update_interval=interval, &
                        inner_m_diag_update_adaptive=.true., &
                        inner_m_diag_update_stall_ratio=stall_ratio, &
                        inner_m_diag_rank_two=.true., &
                        inner_m_diag_rank_two_discover=.true., &
                        inner_m_diag_rank_two_nonsymmetric_base_fun=base_diagonal, &
                        inner_m_diag_rank_two_nonsymmetric_adaptive_transpose_cost=16.0_dp, &
                        jacobian_product=jvp, jacobian_transpose_product=jtv)
                    else
                        fit = root(residual, x0, method='krylov', &
                        inner_method='gmres', inner_maxiter=60, &
                        maxiter=50, fatol=1.0e-10_dp, tol_norm=scaled_norm, &
                        inner_m_diag_auto=.true., &
                        inner_m_diag_update_interval=interval, &
                        inner_m_diag_rank_two=.true., &
                        inner_m_diag_rank_two_discover=.true., &
                        inner_m_diag_rank_two_nonsymmetric_base_fun=base_diagonal, &
                        inner_m_diag_rank_two_nonsymmetric_adaptive_transpose_cost=16.0_dp, &
                        jacobian_product=jvp, jacobian_transpose_product=jtv)
                    end if
                    call cpu_time(end_time)
                    iterations = fit%nit
                    residual_norm = scaled_norm(fit%fun)
                    error = maxval(abs(fit%x-target))
                    if (fit%nfev /= nfev_calls) error stop 'incorrect function evaluation count'
                    if (.not. all(ieee_is_finite(fit%x))) error stop 'nonfinite solution'
                    if (.not. ieee_is_finite(residual_norm)) error stop 'nonfinite residual'
                    if (.not. fit%success .or. residual_norm > 1.0e-9_dp .or. &
                            error > 1.0e-7_dp) failures = failures+1
                    print '(*(g0,:,","))', n, log_scale, rate, trim(modes(ipolicy)), interval, &
                        stall_ratio, work_growth, rebuild_cost, measured, transpose_weight, &
                        fit%success, iterations, fit%nfev, base_calls, jv_calls, &
                        jtv_calls, jv_calls+jtv_calls, jv_calls+4*jtv_calls, &
                        jv_calls+16*jtv_calls, end_time-start_time, residual_norm, error
                end do
            end do
        end do
        deallocate(a, u, v, target, x0, target_value)
    end do
    if (failures > 0) then
        write (*, '(a,i0)') 'failed benchmark configurations: ', failures
        error stop 1
    end if
contains
    subroutine verify_products()
        real(dp) :: direction(n), adjoint_direction(n), forward(n), transpose(n)
        real(dp) :: plus_value(n), minus_value(n), h, discrepancy
        logical :: ok
        integer :: j
        do j = 1, n
            direction(j) = cos(0.47_dp*real(j,dp))
            adjoint_direction(j) = sin(0.37_dp*real(j,dp))
        end do
        h = 2.0e-6_dp
        call jvp(x0, direction, forward, ok)
        if (.not. ok) error stop 'invalid benchmark forward derivative'
        call jtv(x0, adjoint_direction, transpose, ok)
        if (.not. ok) error stop 'invalid benchmark transpose derivative'
        call nonlinear_action(x0+h*direction, plus_value)
        call nonlinear_action(x0-h*direction, minus_value)
        discrepancy = maxval(abs((plus_value-minus_value)/(2.0_dp*h)-forward)/a)
        if (discrepancy > 1.0e-7_dp) error stop 'incorrect benchmark forward derivative'
        discrepancy = abs(dot_product(adjoint_direction,forward) - &
            dot_product(transpose,direction))
        if (discrepancy > 1.0e-11_dp*max(1.0_dp, &
                abs(dot_product(adjoint_direction,forward)))) then
            error stop 'incorrect benchmark transpose derivative'
        end if
    end subroutine verify_products

    subroutine nonlinear_action(x, value)
        real(dp), intent(in) :: x(:) !! current nonlinear state
        real(dp), intent(out) :: value(:) !! unshifted residual
        real(dp) :: s1, s2
        s1 = dot_product(v(:,1),x)
        s2 = dot_product(v(:,2),x)
        value = a*(x+0.5_dp*rate*x**3) + &
            u(:,1)*sin(rate*s1)/rate + &
            u(:,2)*(s2+0.25_dp*rate*s2**3)
    end subroutine nonlinear_action

    subroutine residual(x, fvec)
        real(dp), intent(in) :: x(:) !! current nonlinear iterate
        real(dp), intent(out) :: fvec(:) !! root residual vector
        nfev_calls = nfev_calls+1
        call nonlinear_action(x,fvec)
        fvec = fvec-target_value
    end subroutine residual

    subroutine base_diagonal(x, diagonal, ok)
        real(dp), intent(in) :: x(:) !! current nonlinear iterate
        real(dp), intent(out) :: diagonal(:) !! known uncoupled diagonal D(x)
        logical, intent(out) :: ok !! valid callback output
        base_calls = base_calls+1
        diagonal = a*(1.0_dp+1.5_dp*rate*x*x)
        ok = all(ieee_is_finite(diagonal))
    end subroutine base_diagonal

    subroutine jvp(x, direction, product, ok)
        real(dp), intent(in) :: x(:) !! current nonlinear iterate
        real(dp), intent(in) :: direction(:) !! forward direction
        real(dp), intent(out) :: product(:) !! J(x) times direction
        logical, intent(out) :: ok !! valid product
        real(dp) :: s1, s2
        jv_calls = jv_calls+1
        s1 = dot_product(v(:,1),x)
        s2 = dot_product(v(:,2),x)
        product = a*(1.0_dp+1.5_dp*rate*x*x)*direction + &
            u(:,1)*cos(rate*s1)*dot_product(v(:,1),direction) + &
            u(:,2)*(1.0_dp+0.75_dp*rate*s2*s2)*dot_product(v(:,2),direction)
        ok = all(ieee_is_finite(product))
    end subroutine jvp

    subroutine jtv(x, direction, product, ok)
        real(dp), intent(in) :: x(:) !! current nonlinear iterate
        real(dp), intent(in) :: direction(:) !! transpose direction
        real(dp), intent(out) :: product(:) !! J(x) transpose times direction
        logical, intent(out) :: ok !! valid product
        real(dp) :: s1, s2
        jtv_calls = jtv_calls+1
        s1 = dot_product(v(:,1),x)
        s2 = dot_product(v(:,2),x)
        product = a*(1.0_dp+1.5_dp*rate*x*x)*direction + &
            v(:,1)*cos(rate*s1)*dot_product(u(:,1),direction) + &
            v(:,2)*(1.0_dp+0.75_dp*rate*s2*s2)*dot_product(u(:,2),direction)
        ok = all(ieee_is_finite(product))
    end subroutine jtv

    function scaled_norm(fvec) result(value)
        real(dp), intent(in) :: fvec(:) !! residual to scale by the nominal diagonal
        real(dp) :: value !! maximum row-scaled absolute residual
        value = maxval(abs(fvec)/a)
    end function scaled_norm
end program benchmark_krylov_work_refresh
