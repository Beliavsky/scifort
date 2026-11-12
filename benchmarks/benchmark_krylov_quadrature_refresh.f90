! SPDX-License-Identifier: MIT
! Benchmark adaptive Krylov refresh with actual midpoint-quadrature callbacks.
! The low-rank nonlinear system is identical across refresh policies at a fixed
! quadrature order. OpenMP parallelizes real quadrature work, not delay loops.
program benchmark_krylov_quadrature_refresh
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_kinds, only : dp
    use scifort_optimize_root, only : root, vector_root_result
!$  use omp_lib, only : omp_get_max_threads
    implicit none

    integer, parameter :: n = 128
    integer, parameter :: nmodes = 7
    integer, parameter :: quad_orders(3) = [96, 640, 2048]
    real(dp), parameter :: rates(2) = [0.9_dp, 2.0_dp]
    character(len=24), parameter :: modes(nmodes) = [character(len=24) :: &
        'initial', 'every', 'fixed_four', 'measured_one', 'measured_sixteen', &
        'cpu_timed', 'wall_timed']
    real(dp), parameter :: beta = 0.7_dp
    real(dp) :: a(n), u(n,2), v(n,2), target(n), x0(n), target_value(n)
    real(dp) :: cpu_begin, cpu_end, residual_norm, solution_error, angle, rate
    real(dp) :: p(n), q(n), jp(n), jt(n), fp(n), fm(n), adjoint_error, fd_error
    real(dp) :: adjoint_scale, step
    integer(int64) :: tick_begin, tick_end, tick_rate
    integer :: nquad, threads, imode, irate, iq, i, interval
    integer :: nfev_calls, jv_calls, jtv_calls, base_calls, failures
    real(dp) :: wall_seconds
    logical :: ok
    type(vector_root_result) :: fit

    threads = 1
!$  threads = omp_get_max_threads()
    do i = 1, n
        angle = real(i,dp)
        a(i) = 10.0_dp**(-4.0_dp*real(i-1,dp)/real(n-1,dp))
        u(i,1) = 1.8_dp*a(i)*(0.8_dp+0.22_dp*sin(0.53_dp*angle))
        u(i,2) = 1.8_dp*a(i)*(-0.6_dp+0.17_dp*cos(0.61_dp*angle))
        v(i,1) = (0.9_dp+0.30_dp*cos(0.29_dp*angle))/real(n,dp)
        v(i,2) = (0.7_dp+0.42_dp*sin(0.37_dp*angle))/real(n,dp)
        target(i) = 0.08_dp*cos(0.31_dp*angle)
        x0(i) = target(i)+0.9_dp+0.2_dp*sin(0.23_dp*angle)
        p(i) = cos(0.17_dp*angle)
        q(i) = sin(0.31_dp*angle)
    end do

    failures = 0
    print '(a)', 'threads,quadrature,nonlinearity,mode,success,nit,nfev,jv_calls,'// &
        'jtv_calls,base_calls,cpu_seconds,wall_seconds,scaled_residual,solution_error'
    do iq = 1, size(quad_orders)
        nquad = quad_orders(iq)
        do irate = 1, size(rates)
            rate = rates(irate)
            call check_jacobian()
            call nonlinear_action(target,target_value)
            do imode = 1, nmodes
                interval = 10
                if (imode == 1) interval = 0
                if (imode == 2) interval = 1
                if (imode == 3) interval = 4
                nfev_calls = 0
                jv_calls = 0
                jtv_calls = 0
                base_calls = 0
                call system_clock(tick_begin, tick_rate)
                call cpu_time(cpu_begin)
                select case (imode)
                case (1:3)
                    fit = solve_fixed()
                case (4:5)
                    fit = solve_measured(imode == 5)
                case (6:7)
                    fit = solve_timed(imode == 7)
                end select
                call cpu_time(cpu_end)
                call system_clock(tick_end)
                wall_seconds = -1.0_dp
                if (tick_rate > 0_int64) &
                    wall_seconds = real(tick_end-tick_begin,dp)/real(tick_rate,dp)
                residual_norm = scaled_norm(fit%fun)
                solution_error = maxval(abs(fit%x-target))
                if (fit%nfev /= nfev_calls) error stop 'function evaluation mismatch'
                if (.not. all(ieee_is_finite(fit%x))) error stop 'nonfinite solution'
                if (.not. ieee_is_finite(residual_norm)) error stop 'nonfinite residual'
                if (.not. fit%success .or. residual_norm > 1.0e-9_dp .or. &
                    solution_error > 1.0e-7_dp) failures = failures+1
                print '(*(g0,:,","))', threads, nquad, rate, trim(modes(imode)), &
                    fit%success, fit%nit, fit%nfev, jv_calls, jtv_calls, base_calls, &
                    cpu_end-cpu_begin, wall_seconds, residual_norm, solution_error
            end do
        end do
    end do
    if (failures /= 0) then
        write (*,'(a,i0)') 'failed nonlinear solves: ', failures
        error stop 1
    end if

contains

    function solve_fixed() result(res)
        type(vector_root_result) :: res
        res = root(residual, x0, method='krylov',inner_method='gmres', &
            inner_maxiter=60,maxiter=50,fatol=1.0e-10_dp,tol_norm=scaled_norm, &
            inner_m_diag_auto=.true.,inner_m_diag_update_interval=interval, &
            inner_m_diag_rank_two=.true.,inner_m_diag_rank_two_discover=.true., &
            inner_m_diag_rank_two_nonsymmetric_base_fun=base_diagonal, &
            inner_m_diag_rank_two_nonsymmetric_adaptive_transpose_cost=16.0_dp, &
            jacobian_product=jvp,jacobian_transpose_product=jtv)
    end function solve_fixed

    function solve_measured(expensive) result(res)
        logical, intent(in) :: expensive !! use sixteen-fold relative transpose cost
        real(dp) :: weight
        type(vector_root_result) :: res
        weight = 1.0_dp
        if (expensive) weight = 16.0_dp
        res = root(residual, x0, method='krylov',inner_method='gmres', &
            inner_maxiter=60,maxiter=50,fatol=1.0e-10_dp,tol_norm=scaled_norm, &
            inner_m_diag_auto=.true.,inner_m_diag_update_interval=interval, &
            inner_m_diag_update_adaptive=.true.,inner_m_diag_update_stall_ratio=0.8_dp, &
            inner_m_diag_update_work_growth=1.2_dp, &
            inner_m_diag_update_measured_cost=.true., &
            inner_m_diag_update_transpose_cost=weight, &
            inner_m_diag_rank_two=.true.,inner_m_diag_rank_two_discover=.true., &
            inner_m_diag_rank_two_nonsymmetric_base_fun=base_diagonal, &
            inner_m_diag_rank_two_nonsymmetric_adaptive_transpose_cost=16.0_dp, &
            jacobian_product=jvp,jacobian_transpose_product=jtv)
    end function solve_measured

    function solve_timed(wall_clock) result(res)
        logical, intent(in) :: wall_clock !! choose wall time rather than process CPU time
        type(vector_root_result) :: res
        res = root(residual, x0, method='krylov',inner_method='gmres', &
            inner_maxiter=60,maxiter=50,fatol=1.0e-10_dp,tol_norm=scaled_norm, &
            inner_m_diag_auto=.true.,inner_m_diag_update_interval=interval, &
            inner_m_diag_update_adaptive=.true.,inner_m_diag_update_stall_ratio=0.8_dp, &
            inner_m_diag_update_work_growth=1.2_dp, &
            inner_m_diag_update_timed_cost=.true., &
            inner_m_diag_update_wall_clock=wall_clock, &
            inner_m_diag_update_timed_min_seconds=1.0e-5_dp, &
            inner_m_diag_rank_two=.true.,inner_m_diag_rank_two_discover=.true., &
            inner_m_diag_rank_two_nonsymmetric_base_fun=base_diagonal, &
            inner_m_diag_rank_two_nonsymmetric_adaptive_transpose_cost=16.0_dp, &
            jacobian_product=jvp,jacobian_transpose_product=jtv)
    end function solve_timed

    subroutine local_response(x, local)
        real(dp), intent(in) :: x(:) !! state vector
        real(dp), intent(out) :: local(:) !! midpoint approximation of exp integral minus one
        integer :: j, k
        real(dp) :: t, accum
!$omp parallel do default(shared) private(j,k,t,accum) schedule(static)
        do j = 1, size(x)
            accum = 0.0_dp
            do k = 1, nquad
                t = (real(k,dp)-0.5_dp)/real(nquad,dp)
                accum = accum+exp(t*x(j))-1.0_dp
            end do
            local(j) = accum/real(nquad,dp)
        end do
!$omp end parallel do
    end subroutine local_response

    subroutine local_derivative(x, deriv)
        real(dp), intent(in) :: x(:) !! state vector
        real(dp), intent(out) :: deriv(:) !! midpoint integral of t times exp(t*x)
        integer :: j, k
        real(dp) :: t, accum
!$omp parallel do default(shared) private(j,k,t,accum) schedule(static)
        do j = 1, size(x)
            accum = 0.0_dp
            do k = 1, nquad
                t = (real(k,dp)-0.5_dp)/real(nquad,dp)
                accum = accum+t*exp(t*x(j))
            end do
            deriv(j) = accum/real(nquad,dp)
        end do
!$omp end parallel do
    end subroutine local_derivative

    subroutine nonlinear_action(x, value)
        real(dp), intent(in) :: x(:) !! current state
        real(dp), intent(out) :: value(:) !! nonlinear model response
        real(dp) :: local(n), s1, s2
        call local_response(x,local)
        s1 = dot_product(v(:,1),x)
        s2 = dot_product(v(:,2),x)
        value = a*(x+beta*local+0.5_dp*rate*x**3) + &
            u(:,1)*sin(rate*s1)/rate+u(:,2)*(s2+0.25_dp*rate*s2**3)
    end subroutine nonlinear_action

    subroutine residual(x, value)
        real(dp), intent(in) :: x(:) !! nonlinear iterate
        real(dp), intent(out) :: value(:) !! residual relative to exact target
        nfev_calls = nfev_calls+1
        call nonlinear_action(x,value)
        value = value-target_value
    end subroutine residual

    subroutine base_diagonal(x, diagonal, ok)
        real(dp), intent(in) :: x(:) !! nonlinear iterate
        real(dp), intent(out) :: diagonal(:) !! uncoupled diagonal D(x)
        logical, intent(out) :: ok !! validity indicator
        real(dp) :: deriv(n)
        base_calls = base_calls+1
        call local_derivative(x,deriv)
        diagonal = a*(1.0_dp+beta*deriv+1.5_dp*rate*x*x)
        ok = all(ieee_is_finite(diagonal))
    end subroutine base_diagonal

    subroutine jvp(x, direction, product, ok)
        real(dp), intent(in) :: x(:) !! evaluation point
        real(dp), intent(in) :: direction(:) !! forward direction
        real(dp), intent(out) :: product(:) !! J(x) times direction
        logical, intent(out) :: ok !! validity indicator
        real(dp) :: deriv(n), s1, s2
        jv_calls = jv_calls+1
        call local_derivative(x,deriv)
        s1 = dot_product(v(:,1),x)
        s2 = dot_product(v(:,2),x)
        product = a*(1.0_dp+beta*deriv+1.5_dp*rate*x*x)*direction + &
            u(:,1)*cos(rate*s1)*dot_product(v(:,1),direction) + &
            u(:,2)*(1.0_dp+0.75_dp*rate*s2*s2)*dot_product(v(:,2),direction)
        ok = all(ieee_is_finite(product))
    end subroutine jvp

    subroutine jtv(x, direction, product, ok)
        real(dp), intent(in) :: x(:) !! evaluation point
        real(dp), intent(in) :: direction(:) !! transpose direction
        real(dp), intent(out) :: product(:) !! J(x) transpose times direction
        logical, intent(out) :: ok !! validity indicator
        real(dp) :: deriv(n), s1, s2
        jtv_calls = jtv_calls+1
        call local_derivative(x,deriv)
        s1 = dot_product(v(:,1),x)
        s2 = dot_product(v(:,2),x)
        product = a*(1.0_dp+beta*deriv+1.5_dp*rate*x*x)*direction + &
            v(:,1)*cos(rate*s1)*dot_product(u(:,1),direction) + &
            v(:,2)*(1.0_dp+0.75_dp*rate*s2*s2)*dot_product(u(:,2),direction)
        ok = all(ieee_is_finite(product))
    end subroutine jtv

    subroutine check_jacobian()
        real(dp) :: scale
        step = 2.0e-6_dp
        call jvp(x0,p,jp,ok)
        if (.not. ok) error stop 'Jv returned invalid result'
        call jtv(x0,q,jt,ok)
        if (.not. ok) error stop 'transpose Jv returned invalid result'
        adjoint_scale = max(1.0_dp,abs(dot_product(q,jp)),abs(dot_product(p,jt)))
        adjoint_error = abs(dot_product(q,jp)-dot_product(p,jt))/adjoint_scale
        call nonlinear_action(x0+step*p,fp)
        call nonlinear_action(x0-step*p,fm)
        scale = max(1.0_dp,maxval(abs(jp)))
        fd_error = maxval(abs((fp-fm)/(2.0_dp*step)-jp))/scale
        if (adjoint_error > 1.0e-10_dp .or. fd_error > 1.0e-7_dp) then
            write (*,'(a,2es14.5)') 'Jacobian consistency errors: ', adjoint_error, fd_error
            error stop 'quadrature Jacobian consistency check failed'
        end if
    end subroutine check_jacobian

    function scaled_norm(value) result(norm)
        real(dp), intent(in) :: value(:) !! nonlinear residual
        real(dp) :: norm !! maximum nominally row-scaled residual
        norm = maxval(abs(value)/a)
    end function scaled_norm
end program benchmark_krylov_quadrature_refresh
