! SPDX-License-Identifier: MIT
! Compare exact and deterministic Walsh-probed automatic diagonal estimation.
program benchmark_krylov_diagonal_probes
    use scifort_kinds, only : dp
    use scifort_optimize_root, only : root, vector_root_result
    implicit none

    integer, parameter :: dimensions(3) = [128, 512, 2048]
    integer :: k, run, n, j, jv_calls
    real(dp), allocatable :: diagonal(:), x0(:), target(:)
    real(dp) :: tstart, tend, error_max
    type(vector_root_result) :: res

    print '(a)', 'n,probes,success,nit,nfev,jv_calls,max_error,cpu_seconds'
    do k = 1, size(dimensions)
        n = dimensions(k)
        allocate(diagonal(n), x0(n), target(n))
        do j = 1, n
            diagonal(j) = 10.0_dp**(-5.0_dp * real(j-1, dp) / real(n-1, dp))
            target(j) = 0.2_dp + 0.1_dp * sin(0.05_dp * real(j, dp))
        end do
        x0 = target + 0.25_dp
        do run = 1, 3
            jv_calls = 0
            call cpu_time(tstart)
            if (run == 1) then
                res = root(fun, x0, method='krylov', inner_method='lgmres', &
                    inner_m_diag_auto=.true., inner_m_diag_update_interval=0, &
                    jacobian_product=jvp, inner_maxiter=48, inner_restart=8, &
                    outer_k=0, line_search='none', fatol=1.0e-11_dp, maxiter=15)
            else if (run == 2) then
                res = root(fun, x0, method='krylov', inner_method='lgmres', &
                    inner_m_diag_auto=.true., inner_m_diag_probes=8, &
                    inner_m_diag_update_interval=0, &
                    jacobian_product=jvp, inner_maxiter=48, inner_restart=8, &
                    outer_k=0, line_search='none', fatol=1.0e-11_dp, maxiter=15)
            else
                res = root(fun, x0, method='krylov', inner_method='lgmres', &
                    inner_m_diag_auto=.true., inner_m_diag_probes=16, &
                    inner_m_diag_update_interval=0, &
                    jacobian_product=jvp, inner_maxiter=48, inner_restart=8, &
                    outer_k=0, line_search='none', fatol=1.0e-11_dp, maxiter=15)
            end if
            call cpu_time(tend)
            error_max = maxval(abs(res%x - target))
            print '(i0,a,i0,a,l1,a,i0,a,i0,a,i0,a,es12.4,a,f10.5)', &
                n, ',', merge(0, merge(8,16,run==2),run==1), ',', res%success, ',', &
                res%nit, ',', res%nfev, ',', jv_calls, ',', error_max, ',', tend-tstart
        end do
        deallocate(diagonal, x0, target)
    end do
contains
    subroutine apply_jac(v, prod, x)
        real(dp), intent(in) :: v(:), x(:)
        real(dp), intent(out) :: prod(:)
        integer :: i
        prod = diagonal * (1.0_dp + 0.2_dp * (x-target)) * v
        do i = 1, size(v)-1
            prod(i) = prod(i) + 0.14_dp * diagonal(i) * v(i+1)
            prod(i+1) = prod(i+1) + 0.14_dp * diagonal(i+1) * v(i)
        end do
    end subroutine apply_jac

    subroutine fun(x, f)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: f(:)
        real(dp) :: d(size(x))
        integer :: i
        d = x-target
        f = diagonal * (d + 0.1_dp * d*d)
        do i = 1, size(x)-1
            f(i) = f(i) + 0.14_dp * diagonal(i) * d(i+1)
            f(i+1) = f(i+1) + 0.14_dp * diagonal(i+1) * d(i)
        end do
    end subroutine fun

    subroutine jvp(x, v, prod, ok)
        real(dp), intent(in) :: x(:), v(:)
        real(dp), intent(out) :: prod(:)
        logical, intent(out) :: ok
        jv_calls = jv_calls + 1
        call apply_jac(v, prod, x)
        ok = .true.
    end subroutine jvp
end program benchmark_krylov_diagonal_probes
