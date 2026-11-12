! SPDX-License-Identifier: MIT
! Reproducible dense-coupling diagonal estimator comparison.
program benchmark_krylov_progressive_diagonal
    use scifort_kinds, only : dp
    use scifort_optimize_root, only : root, vector_root_linear_operator, vector_root_result
    implicit none

    real(dp), allocatable :: diagonal(:), u(:), v(:), inverse_estimate(:), x0(:)
    real(dp) :: strength, scale
    integer :: n, j, jv_calls, f_calls
    type(vector_root_result) :: result

    print '(a)', 'n,strength,method,jv_calls,residual_calls,mean_relative_diagonal_error,max_relative_diagonal_error'
    do n = 64, 1024, 192
        allocate(diagonal(n), u(n), v(n), inverse_estimate(n), x0(n))
        do j = 1, n
            diagonal(j) = 1.0_dp + real(j,dp) / real(n,dp)
            u(j) = sin(0.113_dp * real(j, dp)) + 0.5_dp * cos(0.29_dp * real(j, dp))
            v(j) = sin(0.263_dp * real(j, dp)) + 0.2_dp * cos(0.471_dp * real(j, dp))
        end do
        x0 = 0.2_dp
        do j = 1, 2
            if (j == 1) then
                strength = 0.2_dp
            else
                strength = 0.5_dp
            end if
            call run('exact', 0)
            call run('fixed-8', 1)
            call run('progressive-8-to-32', 2)
            call run('fixed-32', 3)
        end do
        deallocate(diagonal, u, v, inverse_estimate, x0)
    end do
contains
    subroutine run(name, mode)
        character(len=*), intent(in) :: name
        integer, intent(in) :: mode
        real(dp) :: diagtrue(size(diagonal)), err(size(diagonal))
        integer :: j

        jv_calls = 0
        f_calls = 0
        select case (mode)
        case (0)
            result = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
                jacobian_product=jvp, inner_solver=inspect)
        case (1)
            result = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
                inner_m_diag_probes=8, inner_m_diag_probe_scheme='rademacher', &
                inner_m_diag_probe_seed=17, inner_m_diag_verify_samples=4, &
                inner_m_diag_verify_rtol=0.25_dp, &
                jacobian_product=jvp, inner_solver=inspect)
        case (2)
            result = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
                inner_m_diag_probes=8, inner_m_diag_probe_scheme='rademacher', &
                inner_m_diag_probe_seed=17, inner_m_diag_verify_samples=4, &
                inner_m_diag_verify_rtol=0.25_dp, inner_m_diag_probe_max=32, &
                jacobian_product=jvp, inner_solver=inspect)
        case (3)
            result = root(residual, x0, method='krylov', inner_m_diag_auto=.true., &
                inner_m_diag_probes=32, inner_m_diag_probe_scheme='rademacher', &
                inner_m_diag_probe_seed=17, inner_m_diag_verify_samples=4, &
                inner_m_diag_verify_rtol=0.25_dp, &
                jacobian_product=jvp, inner_solver=inspect)
        end select
        scale = strength / sqrt(real(size(diagonal), dp))
        diagtrue = diagonal * (1.0_dp + scale * u * v)
        err = abs(inverse_estimate * diagtrue - 1.0_dp)
        print '(i0,",",f4.2,",",a,",",i0,",",i0,",",es12.4,",",es12.4)', &
            size(diagonal), strength, name, jv_calls, f_calls, &
            sum(err) / real(size(err), dp), maxval(err)
    end subroutine run

    subroutine apply_jacobian(input, output)
        real(dp), intent(in) :: input(:)
        real(dp), intent(out) :: output(:)
        real(dp) :: coupling
        coupling = strength * dot_product(v, input) / sqrt(real(size(input), dp))
        output = diagonal * (input + u * coupling)
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
        call apply_jacobian(vec, product)
        ok = .true.
    end subroutine jvp

    subroutine inspect(matvec, preconditioner, rhs, rtol, maxiter, solution, ok)
        procedure(vector_root_linear_operator) :: matvec, preconditioner
        real(dp), intent(in) :: rhs(:), rtol
        integer, intent(in) :: maxiter
        real(dp), intent(out) :: solution(:)
        logical, intent(out) :: ok
        call preconditioner(spread(1.0_dp, 1, size(rhs)), inverse_estimate, ok)
        solution = 0.0_dp
        ok = .false.
    end subroutine inspect
end program benchmark_krylov_progressive_diagonal
