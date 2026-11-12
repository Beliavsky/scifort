! SPDX-License-Identifier: MIT
! Estimated inverse diagonals on an irregularly coupled Jacobian.
! This is a probe-quality benchmark, not an end-to-end LGMRES timing test.
program benchmark_krylov_random_diagonal
    use scifort_kinds, only : dp
    use scifort_optimize_root, only : root, vector_root_result, &
        vector_root_linear_operator
    implicit none

    integer, parameter :: dimensions(3) = [128, 512, 2048]
    integer, parameter :: variants = 5
    integer :: trial, n, i, jv_calls, f_calls
    real(dp), allocatable :: diagonal(:), observed(:), x0(:)
    real(dp) :: tstart, tfinish, error_max, error_mean
    character(len=16) :: mode
    character(len=16), parameter :: names(variants) = [character(len=16) :: &
        'exact', 'walsh', 'rademacher', 'walsh-check', 'random-check']
    type(vector_root_result) :: res

    print '(a)', 'n,method,jv_calls,nfev,max_rel_error,mean_rel_error,cpu_seconds'
    do trial = 1, size(dimensions)
        n = dimensions(trial)
        allocate(diagonal(n), observed(n), x0(n))
        do i = 1, n
            diagonal(i) = 10.0_dp**(-4.0_dp * real(i-1,dp) / real(n-1,dp))
        end do
        x0 = 0.2_dp
        do i = 1, variants
            jv_calls = 0
            f_calls = 0
            mode = names(i)
            call cpu_time(tstart)
            select case (trim(mode))
            case ('exact')
                res = root(residual, x0, method='krylov', &
                    inner_m_diag_auto=.true., jacobian_product=jvp, &
                    inner_solver=inspect)
            case ('walsh')
                res = root(residual, x0, method='krylov', &
                    inner_m_diag_auto=.true., inner_m_diag_probes=8, &
                    jacobian_product=jvp, inner_solver=inspect)
            case ('rademacher')
                res = root(residual, x0, method='krylov', &
                    inner_m_diag_auto=.true., inner_m_diag_probes=8, &
                    inner_m_diag_probe_scheme='rademacher', inner_m_diag_probe_seed=17, &
                    jacobian_product=jvp, inner_solver=inspect)
            case ('walsh-check')
                res = root(residual, x0, method='krylov', &
                    inner_m_diag_auto=.true., inner_m_diag_probes=8, &
                    inner_m_diag_verify_samples=4, inner_m_diag_verify_rtol=0.15_dp, &
                    jacobian_product=jvp, inner_solver=inspect)
            case ('random-check')
                res = root(residual, x0, method='krylov', &
                    inner_m_diag_auto=.true., inner_m_diag_probes=8, &
                    inner_m_diag_probe_scheme='rademacher', inner_m_diag_probe_seed=17, &
                    inner_m_diag_verify_samples=4, inner_m_diag_verify_rtol=0.15_dp, &
                    jacobian_product=jvp, inner_solver=inspect)
            end select
            call cpu_time(tfinish)
            error_max = maxval(abs(diagonal * observed - 1.0_dp))
            error_mean = sum(abs(diagonal * observed - 1.0_dp)) / real(n,dp)
            print '(i0,a,a,a,i0,a,i0,a,es14.6,a,es14.6,a,es14.6)', &
                n, ',', trim(mode), ',', jv_calls, ',', f_calls, ',', &
                error_max, ',', error_mean, ',', tfinish-tstart
        end do
        deallocate(diagonal, observed, x0)
    end do
contains
    subroutine apply_jacobian(v, product)
        real(dp), intent(in) :: v(:)
        real(dp), intent(out) :: product(:)
        integer :: k
        product = diagonal * v
        do k = 1, n-1
            product(k) = product(k) + 0.1_dp * diagonal(k) * v(k+1)
        end do
        do k = 1, n-8
            product(k) = product(k) + 0.45_dp * diagonal(k) * v(k+8)
        end do
        do k = 1, n-19
            product(k) = product(k) + 0.1_dp * diagonal(k) * v(k+19)
        end do
    end subroutine apply_jacobian

    subroutine residual(x, fvec)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: fvec(:)
        f_calls = f_calls + 1
        call apply_jacobian(x, fvec)
    end subroutine residual

    subroutine jvp(x, vec, product, ok)
        real(dp), intent(in) :: x(:), vec(:)
        real(dp), intent(out) :: product(:)
        logical, intent(out) :: ok
        jv_calls = jv_calls + 1
        call apply_jacobian(vec, product)
        ok = .true.
    end subroutine jvp

    subroutine inspect(matvec, preconditioner, rhs, rtol, maxiter, solution, ok)
        procedure(vector_root_linear_operator) :: matvec, preconditioner
        real(dp), intent(in) :: rhs(:), rtol
        integer, intent(in) :: maxiter
        real(dp), intent(out) :: solution(:)
        logical, intent(out) :: ok
        call preconditioner(spread(1.0_dp, 1, size(rhs)), observed, ok)
        if (.not. ok) error stop 'invalid diagonal preconditioner in benchmark'
        solution = 0.0_dp
        ok = .false.
    end subroutine inspect
end program benchmark_krylov_random_diagonal
