! SPDX-License-Identifier: MIT
! Check growing Arnoldi storage, large unused outer capacity, and preconditioning.
program test_optimize_root_krylov_lazy_workspace
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_optimize_root, only : root, vector_root_result
    implicit none

    integer, parameter :: n = 64
    integer, parameter :: short_capacity = 8
    integer, parameter :: long_capacity = 96
    real(dp) :: start(n), target(n), diagonal(n), coupling(n)
    type(vector_root_result) :: short_run, long_run, gmres_run
    integer :: j, failures, nfev_calls, jvp_calls, precond_calls
    integer :: short_jvp, short_precond, short_nfev, mode
    logical :: filter_images

    failures = 0
    do j = 1, n
        target(j) = 1.0_dp + 0.13_dp * sin(0.33_dp * real(j, dp))
        start(j) = target(j) + 0.75_dp + 0.2_dp * cos(0.17_dp * real(j, dp))
        diagonal(j) = 0.4_dp + real(mod(j, 29), dp) / 60.0_dp
        coupling(j) = sin(0.25_dp * real(j, dp))
    end do

    do mode = 1, 2
        filter_images = mode == 2
        call reset_counts()
        short_run = root(residual, start, method='krylov', inner_method='lgmres', &
            inner_m=preconditioner, jacobian_product=jvp, &
            inner_maxiter=8, outer_k=short_capacity, &
            outer_image_prune_auto=filter_images, outer_image_reuse=filter_images, &
            outer_candidate_image_reuse=filter_images, &
            inner_orthogonalization='twice', line_search='none', &
            fatol=1.0e-11_dp, maxiter=30)
        call verify(short_run, 'short outer capacity')
        short_jvp = jvp_calls
        short_precond = precond_calls
        short_nfev = nfev_calls
        call check(short_run%nit >= 5, 'several steps increase the Arnoldi basis')

        call reset_counts()
        long_run = root(residual, start, method='krylov', inner_method='lgmres', &
            inner_m=preconditioner, jacobian_product=jvp, &
            inner_maxiter=8, outer_k=long_capacity, &
            outer_image_prune_auto=filter_images, outer_image_reuse=filter_images, &
            outer_candidate_image_reuse=filter_images, &
            inner_orthogonalization='twice', line_search='none', &
            fatol=1.0e-11_dp, maxiter=30)
        call verify(long_run, 'large unused outer capacity')
        call check(long_run%nit == short_run%nit, 'capacity does not change iterations')
        call check(long_run%nfev == short_run%nfev, 'capacity does not change nfev')
        call check(jvp_calls == short_jvp, 'capacity does not change Jv calls')
        call check(precond_calls == short_precond, 'capacity does not change preconditioner calls')
        call check(nfev_calls == short_nfev, 'capacity preserves residual calls')
        call check(maxval(abs(long_run%x - short_run%x)) < 2.0e-12_dp, &
            'large and small capacities produce the same solution')
    end do

    call reset_counts()
    gmres_run = root(residual, start, method='krylov', inner_method='gmres', &
        inner_m=preconditioner, jacobian_product=jvp, inner_maxiter=12, &
        inner_orthogonalization='selective', line_search='none', &
        fatol=1.0e-11_dp, maxiter=30)
    call verify(gmres_run, 'preconditioned GMRES')
    call check(precond_calls > 0, 'GMRES calls preconditioner')

    if (failures /= 0) then
        print '(a,i0)', 'test_optimize_root_krylov_lazy_workspace: FAIL ', failures
        stop 1
    end if
    print '(a)', 'test_optimize_root_krylov_lazy_workspace: PASS'
contains
    subroutine reset_counts()
        nfev_calls = 0
        jvp_calls = 0
        precond_calls = 0
    end subroutine reset_counts

    subroutine check(valid, label)
        logical, intent(in) :: valid !! check result
        character(len=*), intent(in) :: label !! failure label
        if (.not. valid) then
            failures = failures + 1
            print '(a)', 'FAIL: '//label
        end if
    end subroutine check

    subroutine verify(result, label)
        type(vector_root_result), intent(in) :: result !! root solver result
        character(len=*), intent(in) :: label !! run description
        real(dp) :: check_residual(n)

        call check(result%success, trim(label)//' converged')
        call check(result%nfev == nfev_calls, trim(label)//' nfev accounting')
        call check(jvp_calls > 0, trim(label)//' called Jv')
        call check(precond_calls > 0, trim(label)//' called preconditioner')
        call check(all(ieee_is_finite(result%x)), trim(label)//' finite iterate')
        call model_residual(result%x, check_residual)
        call check(maxval(abs(check_residual)) < 2.0e-11_dp, &
            trim(label)//' true nonlinear residual')
        call check(maxval(abs(result%x - target)) < 2.0e-9_dp, &
            trim(label)//' expected root')
    end subroutine verify

    pure subroutine model_residual(x, fvec)
        real(dp), intent(in) :: x(:) !! nonlinear variables
        real(dp), intent(out) :: fvec(:) !! residual values
        real(dp) :: delta(n), projection

        delta = x - target
        projection = dot_product(coupling, delta) / real(n, dp)
        fvec = diagonal * (x*x - target*target) + 0.03_dp * projection * delta
    end subroutine model_residual

    subroutine residual(x, fvec)
        real(dp), intent(in) :: x(:) !! nonlinear variables
        real(dp), intent(out) :: fvec(:) !! residual values
        nfev_calls = nfev_calls + 1
        call model_residual(x, fvec)
    end subroutine residual

    subroutine jvp(x, vec, product, ok)
        real(dp), intent(in) :: x(:) !! current point
        real(dp), intent(in) :: vec(:) !! multiplication direction
        real(dp), intent(out) :: product(:) !! Jacobian-vector product
        logical, intent(out) :: ok !! success flag
        real(dp) :: delta(n), projection, vprojection

        jvp_calls = jvp_calls + 1
        delta = x - target
        projection = dot_product(coupling, delta) / real(n, dp)
        vprojection = dot_product(coupling, vec) / real(n, dp)
        product = 2.0_dp * diagonal * x * vec + &
            0.03_dp * (projection * vec + vprojection * delta)
        ok = all(ieee_is_finite(product))
    end subroutine jvp

    subroutine preconditioner(vec, product, ok)
        real(dp), intent(in) :: vec(:) !! vector to precondition
        real(dp), intent(out) :: product(:) !! inverse-Jacobian approximation
        logical, intent(out) :: ok !! success flag

        precond_calls = precond_calls + 1
        product = vec / (2.0_dp * diagonal * start)
        ok = all(ieee_is_finite(product))
    end subroutine preconditioner
end program test_optimize_root_krylov_lazy_workspace
