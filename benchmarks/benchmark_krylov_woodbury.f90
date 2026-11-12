! SPDX-License-Identifier: MIT
! Benchmark unpreconditioned, diagonal, and checked low-rank Krylov solves.
! Structured Jacobians are D + u*v^T or D + U*V^T, with known roots.
! CSV work counts are deterministic; CPU times depend on the environment.
program benchmark_krylov_woodbury
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_optimize_root, only : root, vector_root_result
    implicit none

    integer, parameter :: policies = 5
    character(len=16), parameter :: names(policies) = [character(len=16) :: &
        'none', 'diagonal', 'rank_one', 'rank_two', 'rank_two_sketch']
    real(dp), allocatable :: scales(:), target(:), start(:), right(:,:), left(:,:)
    real(dp) :: start_time, end_time, error, residual_norm, coupling_two
    real(dp) :: exponent
    integer :: n, i, p, model, isize, icond, f_calls, jv_calls
    type(vector_root_result) :: fit

    print '(a)', &
        'n,log10_scale,system,policy,success,nit,nfev,njvp,cpu_seconds,scaled_residual,solution_error'
    do isize = 1, 2
        n = 32
        if (isize == 2) n = 128
        allocate(scales(n), target(n), start(n), right(n,2), left(n,2))
        do i = 1, n
            target(i) = 0.3_dp + 0.1_dp*sin(real(i,dp))
            right(i,1) = 1.0_dp
            right(i,2) = cos(0.31_dp*real(i,dp))
            left(i,1) = 0.9_dp + 0.1_dp*cos(real(i,dp))
            left(i,2) = sin(0.23_dp*real(i,dp))
        end do
        start = target + 0.4_dp
        do icond = 1, 2
            exponent = 2.0_dp
            if (icond == 2) exponent = 6.0_dp
            do i = 1, n
                scales(i) = 10.0_dp**(-exponent*real(i-1,dp)/real(n-1,dp))
            end do
            do model = 1, 2
                coupling_two = 0.0_dp
                if (model == 2) coupling_two = 4.0_dp
                do p = 1, policies
                    f_calls = 0
                    jv_calls = 0
                    call cpu_time(start_time)
                    select case(p)
                    case(1)
                        fit = root(residual, start, method='krylov', &
                            inner_method='lgmres', jacobian_product=jvp, &
                            inner_maxiter=60, inner_restart=8, outer_k=0, &
                            line_search='none', tol_norm=scaled_norm, fatol=1.0e-9_dp, &
                            maxiter=12)
                    case(2)
                        fit = root(residual, start, method='krylov', &
                            inner_method='lgmres', jacobian_product=jvp, &
                            inner_maxiter=60, inner_restart=8, outer_k=0, &
                            inner_m_diag_auto=.true., line_search='none', &
                            tol_norm=scaled_norm, fatol=1.0e-9_dp, maxiter=12)
                    case(3)
                        fit = root(residual, start, method='krylov', &
                            inner_method='lgmres', jacobian_product=jvp, &
                            inner_maxiter=60, inner_restart=8, outer_k=0, &
                            inner_m_diag_auto=.true., inner_m_diag_rank_one=.true., &
                            inner_m_diag_rank_one_right=right(:,1), &
                            line_search='none', tol_norm=scaled_norm, &
                            fatol=1.0e-9_dp, maxiter=12)
                    case(4)
                        fit = root(residual, start, method='krylov', &
                            inner_method='lgmres', jacobian_product=jvp, &
                            inner_maxiter=60, inner_restart=8, outer_k=0, &
                            inner_m_diag_auto=.true., inner_m_diag_rank_two=.true., &
                            inner_m_diag_rank_two_right=right, line_search='none', &
                            tol_norm=scaled_norm, fatol=1.0e-9_dp, maxiter=12)
                    case(5)
                        fit = root(residual, start, method='krylov', &
                            inner_method='lgmres', jacobian_product=jvp, &
                            inner_maxiter=60, inner_restart=8, outer_k=0, &
                            inner_m_diag_auto=.true., inner_m_diag_rank_two=.true., &
                            inner_m_diag_rank_two_right=right, &
                            inner_m_diag_rank_two_sketch=.true., line_search='none', &
                            tol_norm=scaled_norm, fatol=1.0e-9_dp, maxiter=12)
                    end select
                    call cpu_time(end_time)
                    if (fit%nfev /= f_calls) error stop 'incorrect residual call accounting'
                    if (.not. all(ieee_is_finite(fit%x))) error stop 'nonfinite solution'
                    residual_norm = scaled_norm(fit%fun)
                    error = maxval(abs(fit%x-target))
                    print '(i0,a,f3.0,a,i0,a,a,a,l1,a,i0,a,i0,a,i0,a,es12.4,a,es12.4,a,es12.4)', &
                        n, ',', exponent, ',', model, ',', trim(names(p)), ',', &
                        fit%success, ',', fit%nit, ',', fit%nfev, ',', jv_calls, ',', &
                        end_time-start_time, ',', residual_norm, ',', error
                end do
            end do
        end do
        deallocate(scales, target, start, right, left)
    end do
contains
    subroutine residual(x, fvec)
        real(dp), intent(in) :: x(:) !! current iterate
        real(dp), intent(out) :: fvec(:) !! scaled coupled residual
        real(dp) :: delta(size(x))
        real(dp) :: moments(2)

        f_calls = f_calls+1
        delta = x-target
        moments = matmul(transpose(right),delta)/real(size(x),dp)
        fvec = scales*(delta + 0.05_dp*delta**2 + &
            8.0_dp*left(:,1)*moments(1) + coupling_two*left(:,2)*moments(2))
    end subroutine residual

    subroutine jvp(x, vec, product, ok)
        real(dp), intent(in) :: x(:) !! linearization point
        real(dp), intent(in) :: vec(:) !! Jacobian direction
        real(dp), intent(out) :: product(:) !! exact Jacobian-direction product
        logical, intent(out) :: ok !! true if the product is finite
        real(dp) :: moments(2)

        jv_calls = jv_calls+1
        moments = matmul(transpose(right),vec)/real(size(x),dp)
        product = scales*((1.0_dp + 0.1_dp*(x-target))*vec + &
            8.0_dp*left(:,1)*moments(1) + coupling_two*left(:,2)*moments(2))
        ok = all(ieee_is_finite(product))
    end subroutine jvp

    function scaled_norm(vec) result(value)
        real(dp), intent(in) :: vec(:) !! row-scaled residual
        real(dp) :: value !! maximum normalized residual
        value = maxval(abs(vec)/scales)
    end function scaled_norm
end program benchmark_krylov_woodbury
