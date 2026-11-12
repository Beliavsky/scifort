! SPDX-License-Identifier: MIT
! Compare full- versus short-subspace matrix-free GMRES in an entire root solve.
!
! Compile using the SciFort FPM dependencies, or use the offline-focused
! build described in docs/validation.md. Example:
!   /usr/bin/time -f 'peak_kib=%M' ./benchmark_gmres_restart full 100000
!   /usr/bin/time -f 'peak_kib=%M' ./benchmark_gmres_restart restart 100000
program benchmark_krylov_restarted_gmres
    use scifort_kinds, only : dp
    use scifort_optimize_root, only : root, vector_root_result
    implicit none

    real(dp), allocatable :: x0(:), diag(:)
    type(vector_root_result) :: r
    character(len=16) :: mode, argument
    real(dp) :: elapsed_start, elapsed_end
    integer :: n, i, jvp_calls, nfev_calls, status

    call get_command_argument(1, mode)
    call get_command_argument(2, argument)
    n = 100000
    if (len_trim(argument) /= 0) then
        read(argument, *, iostat=status) n
        if (status /= 0) error stop 'invalid system size'
    end if
    if (n < 1) error stop 'system size must be positive'
    if (trim(mode) /= 'full' .and. trim(mode) /= 'restart') then
        error stop 'usage: benchmark_krylov_restarted_gmres full|restart [n]'
    end if

    allocate(x0(n), diag(n))
    do i = 1, n
        diag(i) = 1.0_dp + 3.0_dp * real(mod(i, 997), dp) / 997.0_dp
    end do
    x0 = 0.0_dp
    jvp_calls = 0
    nfev_calls = 0

    call cpu_time(elapsed_start)
    if (trim(mode) == 'full') then
        r = root(residual, x0, method='krylov', inner_method='gmres', &
            inner_maxiter=64, jacobian_product=jvp, &
            line_search='none', fatol=1.0e-10_dp, maxiter=15)
    else
        r = root(residual, x0, method='krylov', inner_method='gmres', &
            inner_maxiter=64, inner_restart=8, jacobian_product=jvp, &
            line_search='none', fatol=1.0e-10_dp, maxiter=15)
    end if
    call cpu_time(elapsed_end)
    print '(a,1x,a)', 'mode', trim(mode)
    print '(a,1x,i0)', 'n', n
    print '(a,1x,l1)', 'success', r%success
    print '(a,1x,i0)', 'nonlinear_iterations', r%nit
    print '(a,1x,i0)', 'nfev', r%nfev
    print '(a,1x,i0)', 'analytic_jvp_calls', jvp_calls
    print '(a,1x,es16.8)', 'max_true_residual', maxval(abs(diag * (r%x - 1.0_dp)))
    print '(a,1x,es16.8)', 'seconds', elapsed_end - elapsed_start
    if (.not. r%success .or. r%nfev /= nfev_calls) error stop 'benchmark failed'
contains
    subroutine residual(x, fvec)
        real(dp), intent(in) :: x(:) !! nonlinear iterate
        real(dp), intent(out) :: fvec(:) !! residual
        nfev_calls = nfev_calls + 1
        fvec = diag * (x - 1.0_dp)
    end subroutine

    subroutine jvp(x, vec, product, ok)
        real(dp), intent(in) :: x(:) !! linearization point
        real(dp), intent(in) :: vec(:) !! Krylov direction
        real(dp), intent(out) :: product(:) !! exact Jacobian-vector action
        logical, intent(out) :: ok !! success
        jvp_calls = jvp_calls + 1
        product = diag * vec
        ok = .true.
    end subroutine
end program benchmark_krylov_restarted_gmres
