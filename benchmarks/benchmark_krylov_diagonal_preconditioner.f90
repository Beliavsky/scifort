! SPDX-License-Identifier: MIT
! Compare no preconditioning, a diagonal array, and a callback on an
! ill-scaled diagonal nonlinear problem.  CSV output is reproducible but
! CPU time and work counts are compiler/environment dependent.
program benchmark_krylov_diagonal_preconditioner
    use scifort_kinds, only : dp
    use scifort_optimize_root, only : root, vector_root_result
    implicit none

    integer :: n, i, j, f_calls, jv_calls, precond_calls
    character(len=32) :: arg
    character(len=16), parameter :: modes(3) = [character(len=16) :: &
        'none', 'diag_array', 'callback']
    real(dp), allocatable :: a(:), invdiag(:), target(:), x0(:)
    real(dp) :: t0, t1, residual_error, solution_error
    type(vector_root_result) :: result

    n = 256
    if (command_argument_count() > 0) then
        call get_command_argument(1, arg)
        read(arg, *, iostat=i) n
        if (i /= 0 .or. n < 4) error stop 'expected dimension >= 4'
    end if
    allocate(a(n), invdiag(n), target(n), x0(n))
    do i = 1, n
        a(i) = 10.0_dp ** (-8.0_dp * real(i - 1, dp) / real(n - 1, dp))
        target(i) = 0.8_dp + 0.1_dp * sin(real(i, dp))
    end do
    invdiag = 1.0_dp / a
    x0 = 0.0_dp

    print '(a)', 'n,mode,success,nit,nfev,njvp,nprecond,cpu_seconds,scaled_residual,solution_error'
    do j = 1, size(modes)
        f_calls = 0
        jv_calls = 0
        precond_calls = 0
        call cpu_time(t0)
        select case (j)
        case (1)
            result = root(residual, x0, method='krylov', inner_method='lgmres', &
                inner_maxiter=40, inner_restart=4, inner_restart_adaptive=.true., &
                inner_restart_max=12, outer_k=0, jacobian_product=jvp, &
                tol_norm=scaled_norm, fatol=1.0e-10_dp, &
                maxiter=12, line_search='none')
        case (2)
            result = root(residual, x0, method='krylov', inner_method='lgmres', &
                inner_maxiter=40, inner_restart=4, inner_restart_adaptive=.true., &
                inner_restart_max=12, outer_k=0, jacobian_product=jvp, &
                inner_m_diag=invdiag, tol_norm=scaled_norm, fatol=1.0e-10_dp, &
                maxiter=12, line_search='none')
        case (3)
            result = root(residual, x0, method='krylov', inner_method='lgmres', &
                inner_maxiter=40, inner_restart=4, inner_restart_adaptive=.true., &
                inner_restart_max=12, outer_k=0, jacobian_product=jvp, &
                inner_m=manual_preconditioner, tol_norm=scaled_norm, fatol=1.0e-10_dp, &
                maxiter=12, line_search='none')
        end select
        call cpu_time(t1)
        residual_error = maxval(abs(result%fun) / a)
        solution_error = maxval(abs(result%x - target))
        print '(i0,a,a,a,l1,a,i0,a,i0,a,i0,a,i0,a,es12.4,a,es12.4,a,es12.4)', &
            n, ',', trim(modes(j)), ',', result%success, ',', result%nit, ',', result%nfev, &
            ',', jv_calls, ',', precond_calls, ',', t1-t0, ',', residual_error, ',', solution_error
        if (result%nfev /= f_calls) error stop 'residual evaluation count mismatch'
    end do
contains
    subroutine residual(x, fvec)
        real(dp), intent(in) :: x(:) !! current iterate
        real(dp), intent(out) :: fvec(:) !! residual vector
        f_calls = f_calls + 1
        fvec = a * (x - target)
    end subroutine residual

    subroutine jvp(x, vec, product, ok)
        real(dp), intent(in) :: x(:) !! current iterate
        real(dp), intent(in) :: vec(:) !! direction
        real(dp), intent(out) :: product(:) !! exact Jacobian-vector product
        logical, intent(out) :: ok !! success flag
        jv_calls = jv_calls + 1
        product = a * vec
        ok = .true.
    end subroutine jvp

    subroutine manual_preconditioner(vec, product, ok)
        real(dp), intent(in) :: vec(:) !! vector to precondition
        real(dp), intent(out) :: product(:) !! M * vector
        logical, intent(out) :: ok !! success flag
        precond_calls = precond_calls + 1
        product = invdiag * vec
        ok = .true.
    end subroutine manual_preconditioner

    function scaled_norm(vec) result(value)
        real(dp), intent(in) :: vec(:) !! residual or step
        real(dp) :: value !! maximum component after known row scaling
        value = maxval(abs(vec) / a)
    end function scaled_norm
end program benchmark_krylov_diagonal_preconditioner
