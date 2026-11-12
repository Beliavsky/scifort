! SPDX-License-Identifier: MIT
! Regress persistent Arnoldi workspaces and in-place matvec scratch storage.
program test_optimize_root_krylov_workspace_reuse
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_optimize_root, only : root, vector_root_result
    implicit none

    real(dp), allocatable :: start(:), target(:), diagonal(:)
    type(vector_root_result) :: standard, cached, short_history, gmres
    integer :: failures, nfev_actual, jvp_actual, n, j, k
    integer :: standard_jvp, cached_jvp

    failures = 0
    do k = 1, 2
        if (k == 1) then
            n = 24
        else
            n = 64
        end if
        allocate(start(n), target(n), diagonal(n))
        do j = 1, n
            target(j) = 0.8_dp + 0.2_dp * sin(real(j, dp) * 0.31_dp)
            diagonal(j) = 0.5_dp + real(mod(j, 19), dp) / 30.0_dp
            start(j) = target(j) + 0.8_dp + 0.15_dp * cos(real(j, dp) * 0.21_dp)
        end do

        call reset_counts()
        standard = root(residual, start, method='krylov', inner_method='lgmres', &
            jacobian_product=jvp, inner_m=preconditioner, inner_maxiter=8, &
            outer_k=4, inner_orthogonalization='twice', &
            line_search='none', fatol=1.0e-11_dp, maxiter=35)
        standard_jvp = jvp_actual
        call verify(standard, 'uncached LGMRES')
        call check(standard%nit >= 3, 'several iterations reuse the basis workspace')

        call reset_counts()
        cached = root(residual, start, method='krylov', inner_method='lgmres', &
            jacobian_product=jvp, inner_m=preconditioner, inner_maxiter=8, &
            outer_k=4, outer_image_prune_auto=.true., &
            outer_image_reuse=.true., outer_candidate_image_reuse=.true., &
            inner_orthogonalization='twice', line_search='none', &
            fatol=1.0e-11_dp, maxiter=35)
        cached_jvp = jvp_actual
        call verify(cached, 'cached LGMRES')
        call check(cached_jvp > 0 .and. standard_jvp > 0, &
            'both algorithms apply Jacobian products')
        call check(maxval(abs(standard%x - cached%x)) < 2.0e-9_dp, &
            'caching and repeated workspaces give consistent roots')

        call reset_counts()
        short_history = root(residual, start, method='krylov', inner_method='lgmres', &
            jacobian_product=jvp, inner_m=preconditioner, inner_maxiter=8, &
            outer_k=1, outer_image_prune_auto=.true., &
            outer_image_reuse=.true., outer_candidate_image_reuse=.true., &
            line_search='none', fatol=1.0e-11_dp, maxiter=35)
        call verify(short_history, 'one-vector history')

        call reset_counts()
        gmres = root(residual, start, method='krylov', inner_method='gmres', &
            jacobian_product=jvp, inner_m=preconditioner, inner_maxiter=8, &
            inner_orthogonalization='selective', line_search='none', &
            fatol=1.0e-11_dp, maxiter=35)
        call verify(gmres, 'GMRES residual scratch reuse')
        deallocate(start, target, diagonal)
    end do

    if (failures /= 0) then
        print '(a,i0)', 'test_optimize_root_krylov_workspace_reuse: FAIL ', failures
        stop 1
    end if
    print '(a)', 'test_optimize_root_krylov_workspace_reuse: PASS'

contains
    subroutine check(valid, label)
        logical, intent(in) :: valid !! test outcome
        character(len=*), intent(in) :: label !! assertion name
        if (.not. valid) then
            failures = failures + 1
            print '(a)', 'FAIL: '//label
        end if
    end subroutine check

    subroutine reset_counts()
        nfev_actual = 0
        jvp_actual = 0
    end subroutine reset_counts

    subroutine verify(result, label)
        type(vector_root_result), intent(in) :: result !! solver result
        character(len=*), intent(in) :: label !! assertion prefix
        real(dp) :: exact_residual(size(result%x))

        call check(result%success, trim(label)//': convergence')
        call check(result%nfev == nfev_actual, trim(label)//': residual calls')
        call check(jvp_actual > 0, trim(label)//': Jacobian calls')
        call check(all(ieee_is_finite(result%x)), trim(label)//': finite root')
        exact_residual = diagonal * (result%x*result%x - target*target)
        call check(maxval(abs(exact_residual)) < 2.0e-11_dp, trim(label)//': true residual')
        call check(maxval(abs(result%x - target)) < 2.0e-9_dp, trim(label)//': root accuracy')
    end subroutine verify

    subroutine residual(x, fvec)
        real(dp), intent(in) :: x(:) !! model parameters
        real(dp), intent(out) :: fvec(:) !! nonlinear residual
        nfev_actual = nfev_actual + 1
        fvec = diagonal * (x*x - target*target)
    end subroutine residual

    subroutine jvp(x, vec, product, ok)
        real(dp), intent(in) :: x(:) !! current iterate
        real(dp), intent(in) :: vec(:) !! multiplication direction
        real(dp), intent(out) :: product(:) !! exact Jacobian action
        logical, intent(out) :: ok !! success status
        jvp_actual = jvp_actual + 1
        product = 2.0_dp * diagonal * x * vec
        ok = all(ieee_is_finite(product))
    end subroutine jvp

    subroutine preconditioner(vec, product, ok)
        real(dp), intent(in) :: vec(:) !! vector to transform
        real(dp), intent(out) :: product(:) !! approximate inverse action
        logical, intent(out) :: ok !! finite transformed vector
        product = vec / (2.0_dp * diagonal * start)
        ok = all(ieee_is_finite(product))
    end subroutine preconditioner
end program test_optimize_root_krylov_workspace_reuse
