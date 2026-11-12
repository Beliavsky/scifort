! SPDX-License-Identifier: MIT
! Compare exact-diagonal, dense, streamed, and low-product rank-two discovery.
! A transpose product costs c forward products in the modeled work columns.
! Timings are measured but machine-dependent; modeled work is deterministic.
program benchmark_krylov_rank_two_discovery
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_optimize_root, only : root, vector_root_result
    implicit none

    integer, parameter :: policies = 7
    character(len=20), parameter :: names(policies) = [character(len=20) :: &
        'diagonal', 'dense', 'streamed', 'low_product', 'low_product_checked', &
        'adaptive_1x', 'adaptive_16x']
    real(dp), allocatable :: d(:), u(:,:), v(:,:), target(:), x0(:)
    real(dp) :: start_time, end_time, scaled_error, residual_error, exponent, theta
    integer :: n, i, isize, icond, model, policy, f_calls, jv_calls, jtv_calls
    type(vector_root_result) :: fit

    print '(a)', 'n,log10_scale,rank,policy,success,nit,nfev,njvp,njtv,'// &
        'work_1x,work_4x,work_16x,cpu_seconds,scaled_residual,solution_error'
    do isize = 1, 2
        n = 32
        if (isize == 2) n = 128
        allocate(d(n), u(n,3), v(n,3), target(n), x0(n))
        do icond = 1, 2
            exponent = 2.0_dp
            if (icond == 2) exponent = 6.0_dp
            do i = 1, n
                theta = real(i,dp)
                d(i) = 10.0_dp**(-exponent*real(i-1,dp)/real(n-1,dp))
                u(i,1) = d(i)*(1.2_dp+0.21_dp*sin(0.17_dp*theta))
                u(i,2) = d(i)*(-0.8_dp+0.23_dp*cos(0.29_dp*theta))
                u(i,3) = d(i)*(0.9_dp+0.27_dp*sin(0.53_dp*theta))
                v(i,1) = (0.6_dp+0.35_dp*cos(0.31_dp*theta))/real(n,dp)
                v(i,2) = (0.4_dp+0.43_dp*sin(0.41_dp*theta))/real(n,dp)
                v(i,3) = (0.8_dp+0.37_dp*cos(0.67_dp*theta))/real(n,dp)
                target(i) = 0.05_dp*cos(0.13_dp*theta)
            end do
            x0 = target + 0.3_dp
            do model = 2, 3
                do policy = 1, policies
                    f_calls = 0
                    jv_calls = 0
                    jtv_calls = 0
                    call cpu_time(start_time)
                    select case(policy)
                    case(1)
                        fit = root(residual, x0, method='krylov', inner_method='lgmres', &
                            jacobian_product=jvp, inner_maxiter=60, inner_restart=12, &
                            outer_k=0, inner_m_diag_auto=.true., &
                            tol_norm=scaled_norm, fatol=1.0e-9_dp, maxiter=7, &
                            line_search='none')
                    case(2)
                        fit = root(residual, x0, method='krylov', inner_method='lgmres', &
                            jacobian_product=jvp, inner_maxiter=60, inner_restart=12, &
                            outer_k=0, inner_m_diag_auto=.true., &
                            inner_m_diag_rank_two=.true., &
                            inner_m_diag_rank_two_discover=.true., &
                            tol_norm=scaled_norm, fatol=1.0e-9_dp, maxiter=7, &
                            line_search='none')
                    case(3)
                        fit = root(residual, x0, method='krylov', inner_method='lgmres', &
                            jacobian_product=jvp, inner_maxiter=60, inner_restart=12, &
                            outer_k=0, inner_m_diag_auto=.true., &
                            inner_m_diag_rank_two=.true., &
                            inner_m_diag_rank_two_discover=.true., &
                            inner_m_diag_rank_two_stream=.true., &
                            tol_norm=scaled_norm, fatol=1.0e-9_dp, maxiter=7, &
                            line_search='none')
                    case(4,5)
                        fit = root(residual, x0, method='krylov', inner_method='lgmres', &
                            jacobian_product=jvp, jacobian_transpose_product=jtv, &
                            inner_maxiter=60, inner_restart=12, &
                            outer_k=0, inner_m_diag_auto=.true., &
                            inner_m_diag_rank_two=.true., &
                            inner_m_diag_rank_two_discover=.true., &
                            inner_m_diag_rank_two_nonsymmetric_base=d, &
                            inner_m_diag_rank_two_nonsymmetric_verify_transpose=policy==5, &
                            tol_norm=scaled_norm, fatol=1.0e-9_dp, maxiter=7, &
                            line_search='none')
                    case(6,7)
                        fit = root(residual, x0, method='krylov', inner_method='lgmres', &
                            jacobian_product=jvp, jacobian_transpose_product=jtv, &
                            inner_maxiter=60, inner_restart=12, &
                            outer_k=0, inner_m_diag_auto=.true., &
                            inner_m_diag_rank_two=.true., &
                            inner_m_diag_rank_two_discover=.true., &
                            inner_m_diag_rank_two_nonsymmetric_base=d, &
                            inner_m_diag_rank_two_nonsymmetric_adaptive_transpose_cost= &
                                real(1+15*(policy-6),dp), &
                            tol_norm=scaled_norm, fatol=1.0e-9_dp, maxiter=7, &
                            line_search='none')
                    end select
                    call cpu_time(end_time)
                    if (fit%nfev /= f_calls) error stop 'incorrect evaluation count'
                    if (.not. all(ieee_is_finite(fit%x))) error stop 'nonfinite solution'
                    scaled_error = maxval(abs(fit%x-target))
                    residual_error = scaled_norm(fit%fun)
                    print '(*(g0,:,","))', &
                        n, exponent, model, trim(names(policy)), &
                        fit%success, fit%nit, fit%nfev, &
                        jv_calls, jtv_calls, jv_calls+jtv_calls, &
                        jv_calls+4*jtv_calls, jv_calls+16*jtv_calls, &
                        end_time-start_time, residual_error, scaled_error
                end do
            end do
        end do
        deallocate(d, u, v, target, x0)
    end do
contains
    subroutine residual(x, fvec)
        real(dp), intent(in) :: x(:) !! current iterate
        real(dp), intent(out) :: fvec(:) !! linear root residual
        f_calls = f_calls+1
        call action(x-target,fvec)
    end subroutine residual

    subroutine action(vec, product)
        real(dp), intent(in) :: vec(:) !! forward direction
        real(dp), intent(out) :: product(:) !! exact Jv product
        integer :: k
        product = d*vec
        do k = 1, model
            product = product + u(:,k)*dot_product(v(:,k),vec)
        end do
    end subroutine action

    subroutine jvp(x, vec, product, ok)
        real(dp), intent(in) :: x(:) !! current linearization point
        real(dp), intent(in) :: vec(:) !! forward multiplication direction
        real(dp), intent(out) :: product(:) !! Jacobian-vector product
        logical, intent(out) :: ok !! product is finite
        jv_calls = jv_calls + 1
        call action(vec,product)
        ok = all(ieee_is_finite(product))
    end subroutine jvp

    subroutine jtv(x, vec, product, ok)
        real(dp), intent(in) :: x(:) !! current linearization point
        real(dp), intent(in) :: vec(:) !! transpose multiplication direction
        real(dp), intent(out) :: product(:) !! transpose Jacobian-vector product
        logical, intent(out) :: ok !! product is finite
        integer :: k
        jtv_calls = jtv_calls + 1
        product = d*vec
        do k = 1, model
            product = product + v(:,k)*dot_product(u(:,k),vec)
        end do
        ok = all(ieee_is_finite(product))
    end subroutine jtv

    function scaled_norm(fvec) result(value)
        real(dp), intent(in) :: fvec(:) !! diagonal-scaled residual
        real(dp) :: value !! maximum row-normalized norm
        value = maxval(abs(fvec)/d)
    end function scaled_norm
end program benchmark_krylov_rank_two_discovery
