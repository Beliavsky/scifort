! SPDX-License-Identifier: MIT
! Matrix-free, ill-conditioned nonlinear LGMRES restart regression.
! At the exact root the Jacobian is D*Q, where Q is orthogonal, and
! cond_2(D*Q) equals the prescribed diagonal ratio.
program test_optimize_root_krylov_restart_scaling
    use scifort_kinds, only : dp
    use scifort_optimize_root, only : root, vector_root_result
    implicit none

    integer :: n, nhalf, i, nfev_calls, jvp_calls, precond_calls
    integer :: failures
    real(dp) :: condition, vscale
    real(dp), allocatable :: diagonal(:), target(:), x0(:)
    real(dp), allocatable :: vec(:), product(:), plus(:), minus(:)
    type(vector_root_result) :: rraw, rscaled, radaptive
    logical :: ok

    failures = 0
    do n = 64, 256, 192
        nhalf = n / 2
        allocate(diagonal(n), target(n), x0(n), vec(n), product(n), plus(n), minus(n))
        condition = 1.0e6_dp
        call setup()

        vec = [(sin(real(i, dp)), i=1,n)]
        call jvp(x0, vec, product, ok)
        vscale = 1.0e-5_dp
        call model_residual(x0 + vscale*vec, plus)
        call model_residual(x0 - vscale*vec, minus)
        call check(ok, 'analytic Jacobian product is defined')
        call check(maxval(abs(((plus-minus)/(2.0_dp*vscale) - product) / diagonal)) < &
            1.0e-8_dp, 'analytic Jacobian product matches centered differences')

        call reset_counts()
        rraw = root(residual, x0, method='krylov', inner_method='lgmres', &
            inner_maxiter=48, outer_k=2, inner_m=preconditioner, &
            jacobian_product=jvp, line_search='none', fatol=1.0e-8_dp, maxiter=12)
        call check(rraw%success, 'full-width LGMRES meets unscaled tolerance')
        call check(rraw%nfev == nfev_calls, 'unscaled residual work count')
        call check(jvp_calls > 0 .and. precond_calls > 0, 'Jv and preconditioner dispatched')
        call check(maxval(abs(rraw%fun)) < 1.0e-8_dp, 'unscaled final residual')

        call reset_counts()
        rscaled = root(residual, x0, method='krylov', inner_method='lgmres', &
            inner_maxiter=48, outer_k=2, inner_m=preconditioner, &
            jacobian_product=jvp, tol_norm=scaled_norm, &
            line_search='none', fatol=1.0e-8_dp, maxiter=12)
        call check(rscaled%success, 'full-width LGMRES meets scaled tolerance')
        call check(rscaled%nfev == nfev_calls, 'scaled residual work count')
        call check(maxval(abs(rscaled%fun / diagonal)) < 1.0e-8_dp, &
            'full-width scaled residual is accurate')
        call check(maxval(abs(rscaled%x - target)) < 1.0e-6_dp, &
            'full-width solution is accurate after scaling')

        call reset_counts()
        radaptive = root(residual, x0, method='krylov', inner_method='lgmres', &
            inner_maxiter=48, inner_restart=6, inner_restart_adaptive=.true., &
            inner_restart_max=24, outer_k=2, inner_m=preconditioner, &
            jacobian_product=jvp, tol_norm=scaled_norm, &
            line_search='none', fatol=1.0e-8_dp, maxiter=12)
        call check(radaptive%success, 'adaptive LGMRES meets scaled tolerance')
        call check(radaptive%nfev == nfev_calls, 'adaptive residual work count')
        call check(maxval(abs(radaptive%fun / diagonal)) < 1.0e-8_dp, &
            'adaptive scaled residual is accurate')
        call check(maxval(abs(radaptive%x - target)) < 1.0e-6_dp, &
            'adaptive solution is accurate after scaling')
        deallocate(diagonal, target, x0, vec, product, plus, minus)
    end do

    if (failures /= 0) then
        print '(a,i0)', 'test_optimize_root_krylov_restart_scaling: FAIL ', failures
        stop 1
    end if
    print '(a)', 'test_optimize_root_krylov_restart_scaling: PASS'
contains
    subroutine check(valid, message)
        logical, intent(in) :: valid
        character(len=*), intent(in) :: message
        if (.not. valid) then
            failures = failures + 1
            print '(a)', 'FAIL: '//message
        end if
    end subroutine

    subroutine reset_counts()
        nfev_calls = 0
        jvp_calls = 0
        precond_calls = 0
    end subroutine

    subroutine setup()
        integer :: k
        do k = 1, n
            diagonal(k) = condition**(-real(k-1, dp) / real(n-1, dp))
            target(k) = 0.4_dp * sin(0.025_dp * real(k, dp))
            x0(k) = target(k) + 0.12_dp * cos(0.031_dp * real(k, dp))
        end do
    end subroutine

    subroutine model_residual(x, fvec)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: fvec(:)
        integer :: k
        real(dp) :: v1, v2
        do k = 1, nhalf
            v1 = x(k) - target(k)
            v2 = x(k+nhalf) - target(k+nhalf)
            fvec(k) = diagonal(k) * (cos(0.65_dp)*v1 - sin(0.65_dp)*v2 + 0.04_dp*v1*v1)
            fvec(k+nhalf) = diagonal(k+nhalf) * &
                (sin(0.65_dp)*v1 + cos(0.65_dp)*v2 + 0.04_dp*v2*v2)
        end do
    end subroutine

    subroutine residual(x, fvec)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: fvec(:)
        nfev_calls = nfev_calls + 1
        call model_residual(x, fvec)
    end subroutine

    subroutine jvp(x, vec, product, ok)
        real(dp), intent(in) :: x(:)
        real(dp), intent(in) :: vec(:)
        real(dp), intent(out) :: product(:)
        logical, intent(out) :: ok
        integer :: k
        jvp_calls = jvp_calls + 1
        do k = 1, nhalf
            product(k) = diagonal(k) * &
                (cos(0.65_dp)*vec(k) - sin(0.65_dp)*vec(k+nhalf) + &
                0.08_dp*(x(k)-target(k))*vec(k))
            product(k+nhalf) = diagonal(k+nhalf) * &
                (sin(0.65_dp)*vec(k) + cos(0.65_dp)*vec(k+nhalf) + &
                0.08_dp*(x(k+nhalf)-target(k+nhalf))*vec(k+nhalf))
        end do
        ok = .true.
    end subroutine

    subroutine preconditioner(vec, product, ok)
        real(dp), intent(in) :: vec(:)
        real(dp), intent(out) :: product(:)
        logical, intent(out) :: ok
        precond_calls = precond_calls + 1
        product = vec / sqrt(diagonal)
        ok = .true.
    end subroutine

    function scaled_norm(values) result(value)
        real(dp), intent(in) :: values(:)
        real(dp) :: value
        value = maxval(abs(values / diagonal))
    end function
end program test_optimize_root_krylov_restart_scaling
