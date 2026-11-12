! SPDX-License-Identifier: MIT
! Deterministic, matrix-free benchmark of LGMRES restart policies.
!
! The linear part is D*Q, with Q an orthogonal rotation coupling opposite
! halves of the vector and D diagonal. Thus cond_2(D*Q)=cond_2(D) exactly.
! The additional nonlinear term is 0.04*D*(x-target)**2, giving a Jacobian
! that changes at each nonlinear iterate. There is no matrix allocation.
!
! Usage: benchmark_krylov_restart_scaling N CONDITION POLICY [PRECONDITION] [NORM]
! POLICY = full | fixed | adaptive; PRECONDITION = yes | no (default yes).
! NORM = raw | scaled (default raw); scaled uses norm(max abs(F_i / D_i)).
! The program writes one CSV row. Measure process peak RSS externally with
! /usr/bin/time -f '%M' on Linux, or an equivalent process monitor.
program benchmark_krylov_restart_scaling
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_optimize_root, only : root, vector_root_result
    implicit none

    real(dp), parameter :: gamma = 0.04_dp
    real(dp), parameter :: theta = 0.65_dp
    integer, parameter :: total_width = 48
    integer, parameter :: first_width = 6
    integer, parameter :: max_width = 24
    integer :: n, nhalf, i, jvp_count, residual_count, precond_count
    real(dp), allocatable :: diagonal(:), target(:), x0(:)
    real(dp) :: condition, t0, t1, residual_norm, scaled_residual, solution_error, diagonal_exp
    character(len=40) :: arg, policy, precondition, norm_mode
    type(vector_root_result) :: result

    if (command_argument_count() < 3) then
        print '(a)', 'usage: benchmark_krylov_restart_scaling N CONDITION full|fixed|adaptive [yes|no] [raw|scaled]'
        stop 2
    end if
    call get_command_argument(1, arg)
    read(arg, *) n
    call get_command_argument(2, arg)
    read(arg, *) condition
    call get_command_argument(3, policy)
    precondition = 'yes'
    if (command_argument_count() >= 4) call get_command_argument(4, precondition)
    norm_mode = 'raw'
    if (command_argument_count() >= 5) call get_command_argument(5, norm_mode)
    if (n < 4 .or. mod(n, 2) /= 0 .or. condition < 1.0_dp .or. &
            (trim(policy) /= 'full' .and. trim(policy) /= 'fixed' .and. &
            trim(policy) /= 'adaptive') .or. &
            (trim(precondition) /= 'yes' .and. trim(precondition) /= 'no') .or. &
            (trim(norm_mode) /= 'raw' .and. trim(norm_mode) /= 'scaled')) &
        error stop 'invalid benchmark arguments'

    nhalf = n / 2
    allocate(diagonal(n), target(n), x0(n))
    do i = 1, n
        diagonal_exp = real(i - 1, dp) / real(n - 1, dp)
        diagonal(i) = condition**(-diagonal_exp)
        target(i) = 0.4_dp * sin(0.025_dp * real(i, dp))
        x0(i) = target(i) + 0.12_dp * cos(0.031_dp * real(i, dp))
    end do

    jvp_count = 0
    residual_count = 0
    precond_count = 0
    call cpu_time(t0)
    if (trim(precondition) == 'yes') then
        select case (trim(policy))
        case ('full')
            result = root(residual, x0, method='krylov', inner_method='lgmres', &
                inner_maxiter=total_width, outer_k=2, jacobian_product=jvp, &
                inner_m=preconditioner, line_search='none', fatol=1.0e-8_dp, maxiter=12, tol_norm=benchmark_norm)
        case ('fixed')
            result = root(residual, x0, method='krylov', inner_method='lgmres', &
                inner_maxiter=total_width, inner_restart=first_width, outer_k=2, &
                jacobian_product=jvp, inner_m=preconditioner, &
                line_search='none', fatol=1.0e-8_dp, maxiter=12, tol_norm=benchmark_norm)
        case ('adaptive')
            result = root(residual, x0, method='krylov', inner_method='lgmres', &
                inner_maxiter=total_width, inner_restart=first_width, &
                inner_restart_adaptive=.true., inner_restart_max=max_width, &
                outer_k=2, jacobian_product=jvp, inner_m=preconditioner, &
                line_search='none', fatol=1.0e-8_dp, maxiter=12, tol_norm=benchmark_norm)
        end select
    else
        select case (trim(policy))
        case ('full')
            result = root(residual, x0, method='krylov', inner_method='lgmres', &
                inner_maxiter=total_width, outer_k=2, jacobian_product=jvp, &
                line_search='none', fatol=1.0e-8_dp, maxiter=12, tol_norm=benchmark_norm)
        case ('fixed')
            result = root(residual, x0, method='krylov', inner_method='lgmres', &
                inner_maxiter=total_width, inner_restart=first_width, outer_k=2, &
                jacobian_product=jvp, &
                line_search='none', fatol=1.0e-8_dp, maxiter=12, tol_norm=benchmark_norm)
        case ('adaptive')
            result = root(residual, x0, method='krylov', inner_method='lgmres', &
                inner_maxiter=total_width, inner_restart=first_width, &
                inner_restart_adaptive=.true., inner_restart_max=max_width, &
                outer_k=2, jacobian_product=jvp, &
                line_search='none', fatol=1.0e-8_dp, maxiter=12, tol_norm=benchmark_norm)
        end select
    end if
    call cpu_time(t1)
    if (result%nfev /= residual_count) error stop 'inaccurate residual evaluation count'
    if (.not. allocated(result%fun)) error stop 'missing final residual'
    if (.not. allocated(result%x)) error stop 'missing final solution'
    if (.not. all(ieee_is_finite(result%fun))) error stop 'nonfinite final residual'
    residual_norm = maxval(abs(result%fun))
    scaled_residual = maxval(abs(result%fun / diagonal))
    solution_error = maxval(abs(result%x - target))
    ! The independent final residual is not charged to solver nfev.
    block
        real(dp), allocatable :: check(:)
        allocate(check(n))
        call model_residual(result%x, check)
        if (maxval(abs((check - result%fun) / diagonal)) > &
                1.0e-12_dp * max(1.0_dp, scaled_residual)) &
            error stop 'incorrect final residual'
    end block
    print '(a)', 'n,condition,precondition,norm,policy,success,nit,nfev,njvp,nprec,maxres,scaledres,maxerr,cpu_s'
    print '(i0,",",es12.4,",",a,",",a,",",a,",",l1,",",i0,",",i0,",",i0,",",i0,",",es13.5,",",es13.5,",",es13.5,",",f12.5)', &
        n, condition, trim(precondition), trim(norm_mode), trim(policy), &
        result%success, result%nit, result%nfev, jvp_count, precond_count, &
        residual_norm, scaled_residual, solution_error, t1 - t0
contains
    function benchmark_norm(values) result(value)
        real(dp), intent(in) :: values(:)
        real(dp) :: value
        if (trim(norm_mode) == 'scaled') then
            value = maxval(abs(values / diagonal))
        else
            value = maxval(abs(values))
        end if
    end function benchmark_norm

    subroutine model_residual(x, fvec)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: fvec(:)
        integer :: k
        real(dp) :: v1, v2
        do k = 1, nhalf
            v1 = x(k) - target(k)
            v2 = x(k + nhalf) - target(k + nhalf)
            fvec(k) = diagonal(k) * (cos(theta)*v1 - sin(theta)*v2 + gamma*v1*v1)
            fvec(k + nhalf) = diagonal(k + nhalf) * &
                (sin(theta)*v1 + cos(theta)*v2 + gamma*v2*v2)
        end do
    end subroutine

    subroutine residual(x, fvec)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: fvec(:)
        residual_count = residual_count + 1
        call model_residual(x, fvec)
    end subroutine

    subroutine jvp(x, vec, product, ok)
        real(dp), intent(in) :: x(:)
        real(dp), intent(in) :: vec(:)
        real(dp), intent(out) :: product(:)
        logical, intent(out) :: ok
        integer :: k
        real(dp) :: v1, v2
        jvp_count = jvp_count + 1
        do k = 1, nhalf
            v1 = vec(k)
            v2 = vec(k + nhalf)
            product(k) = diagonal(k) * &
                (cos(theta)*v1 - sin(theta)*v2 + 2.0_dp*gamma*(x(k)-target(k))*v1)
            product(k + nhalf) = diagonal(k + nhalf) * &
                (sin(theta)*v1 + cos(theta)*v2 + &
                 2.0_dp*gamma*(x(k+nhalf)-target(k+nhalf))*v2)
        end do
        ok = .true.
    end subroutine

    subroutine preconditioner(vec, product, ok)
        real(dp), intent(in) :: vec(:)
        real(dp), intent(out) :: product(:)
        logical, intent(out) :: ok
        precond_count = precond_count + 1
        product = vec / sqrt(diagonal)
        ok = .true.
    end subroutine
end program benchmark_krylov_restart_scaling
