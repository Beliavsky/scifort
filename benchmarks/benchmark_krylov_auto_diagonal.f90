! SPDX-License-Identifier: MIT
! Reproducible coupled nonlinear systems for automatic Jacobian-diagonal preconditioning.
!
! Compile with SciFort's normal FPM build or link its optimize root module.
! The output is CSV, suitable for a table or plotting tool.
program benchmark_krylov_auto_diagonal
    use scifort_kinds, only : dp
    use scifort_optimize_root, only : root, vector_root_result
    implicit none

    real(dp), allocatable :: scales(:), targets(:), start(:), fixed_inv(:)
    real(dp) :: tic, toc, error
    integer :: n, i, mode, f_calls, jv_calls
    type(vector_root_result) :: fit
    character(len=24) :: label

    print '(a)', 'size,policy,success,iterations,nfev,jv_calls,cpu_seconds,max_solution_error'
    do mode = 1, 2
        n = 32
        if (mode == 2) n = 128
        allocate(scales(n), targets(n), start(n), fixed_inv(n))
        do i = 1, n
            scales(i) = 10.0_dp**(-6.0_dp * real(i-1, dp) / real(n-1, dp))
            targets(i) = 0.6_dp + 0.1_dp * cos(real(i, dp))
        end do
        start = targets + 0.6_dp
        fixed_inv = 1.0_dp / (scales * (1.0_dp + 0.4_dp*(start-targets)))
        do i = 1, 4
            f_calls = 0
            jv_calls = 0
            call cpu_time(tic)
            select case(i)
            case(1)
                label = 'none'
                fit = root(residual, start, method='krylov', inner_method='lgmres', &
                    inner_maxiter=64, inner_restart=8, outer_k=2, &
                    jacobian_product=jvp, line_search='none', &
                    tol_norm=scaled_max_norm, fatol=1.0e-9_dp, maxiter=18)
            case(2)
                label = 'static_initial'
                fit = root(residual, start, method='krylov', inner_method='lgmres', &
                    inner_maxiter=64, inner_restart=8, outer_k=2, &
                    inner_m_diag=fixed_inv, jacobian_product=jvp, line_search='none', &
                    tol_norm=scaled_max_norm, fatol=1.0e-9_dp, maxiter=18)
            case(3)
                label = 'auto_initial'
                fit = root(residual, start, method='krylov', inner_method='lgmres', &
                    inner_maxiter=64, inner_restart=8, outer_k=2, &
                    inner_m_diag_auto=.true., inner_m_diag_update_interval=0, &
                    jacobian_product=jvp, line_search='none', &
                    tol_norm=scaled_max_norm, fatol=1.0e-9_dp, maxiter=18)
            case(4)
                label = 'auto_every_step'
                fit = root(residual, start, method='krylov', inner_method='lgmres', &
                    inner_maxiter=64, inner_restart=8, outer_k=2, &
                    inner_m_diag_auto=.true., inner_m_diag_update_interval=1, &
                    jacobian_product=jvp, line_search='none', &
                    tol_norm=scaled_max_norm, fatol=1.0e-9_dp, maxiter=18)
            end select
            call cpu_time(toc)
            error = maxval(abs(fit%x-targets))
            if (fit%nfev /= f_calls) error stop 'residual callback accounting mismatch'
            print '(i0,a,a,a,l1,a,i0,a,i0,a,i0,a,es13.5,a,es13.5)', &
                n, ',', trim(label), ',', fit%success, ',', fit%nit, ',', &
                fit%nfev, ',', jv_calls, ',', toc-tic, ',', error
        end do
        deallocate(scales, targets, start, fixed_inv)
    end do
contains
    subroutine residual(x, fvec)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: fvec(:)
        real(dp) :: e(size(x))
        integer :: j, next
        f_calls = f_calls + 1
        e = x - targets
        do j = 1, size(x)
            next = mod(j, size(x)) + 1
            fvec(j) = scales(j)*(e(j) + 0.2_dp*e(j)**2 + 0.06_dp*e(next))
        end do
    end subroutine residual

    subroutine jvp(x, vec, product, ok)
        real(dp), intent(in) :: x(:)
        real(dp), intent(in) :: vec(:)
        real(dp), intent(out) :: product(:)
        logical, intent(out) :: ok
        integer :: j, next
        jv_calls = jv_calls + 1
        do j = 1, size(x)
            next = mod(j, size(x)) + 1
            product(j) = scales(j)*((1.0_dp+0.4_dp*(x(j)-targets(j)))*vec(j) + &
                0.06_dp*vec(next))
        end do
        ok = .true.
    end subroutine jvp

    function scaled_max_norm(vec) result(norm_value)
        real(dp), intent(in) :: vec(:)
        real(dp) :: norm_value
        norm_value = maxval(abs(vec)/scales)
    end function scaled_max_norm
end program benchmark_krylov_auto_diagonal
