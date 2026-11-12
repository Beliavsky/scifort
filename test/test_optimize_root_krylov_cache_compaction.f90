! SPDX-License-Identifier: MIT
! Exercise both image caches with retained-vector compaction and preconditioning.
program test_optimize_root_krylov_cache_compaction
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_optimize_root, only : root, vector_root_result
    implicit none

    integer, parameter :: n = 64
    real(dp) :: initial(n), target(n), diagonal(n), coupling(n)
    type(vector_root_result) :: base, retained, candidate, both, small_history
    integer :: failures, calls_base, calls_retained, calls_candidate
    integer :: calls_both, nfev_actual, jvp_calls, j

    failures = 0
    do j = 1, n
        target(j) = 1.0_dp + 0.12_dp * sin(0.5_dp * real(j, dp))
        initial(j) = target(j) + 0.7_dp + 0.2_dp * sin(0.17_dp * real(j, dp))
        diagonal(j) = 0.3_dp + real(mod(j, 23), dp) / 55.0_dp
        coupling(j) = cos(0.27_dp * real(j, dp))
    end do

    call reset()
    base = root(residual, initial, method='krylov', inner_method='lgmres', &
        inner_m=preconditioner, jacobian_product=jvp, inner_maxiter=8, &
        outer_k=4, outer_image_prune_auto=.true., &
        inner_orthogonalization='twice', fatol=1.0e-11_dp, maxiter=50)
    calls_base = jvp_calls
    call verify(base, 'uncached reference')

    call reset()
    retained = root(residual, initial, method='krylov', inner_method='lgmres', &
        inner_m=preconditioner, jacobian_product=jvp, inner_maxiter=8, &
        outer_k=4, outer_image_prune_auto=.true., outer_image_reuse=.true., &
        inner_orthogonalization='twice', fatol=1.0e-11_dp, maxiter=50)
    calls_retained = jvp_calls
    call verify(retained, 'retained-image reuse')

    call reset()
    candidate = root(residual, initial, method='krylov', inner_method='lgmres', &
        inner_m=preconditioner, jacobian_product=jvp, inner_maxiter=8, &
        outer_k=4, outer_image_prune_auto=.true., &
        outer_candidate_image_reuse=.true., &
        inner_orthogonalization='twice', fatol=1.0e-11_dp, maxiter=50)
    calls_candidate = jvp_calls
    call verify(candidate, 'candidate-image reuse')

    call reset()
    both = root(residual, initial, method='krylov', inner_method='lgmres', &
        inner_m=preconditioner, jacobian_product=jvp, inner_maxiter=8, &
        outer_k=4, outer_image_prune_auto=.true., outer_image_reuse=.true., &
        outer_candidate_image_reuse=.true., &
        inner_orthogonalization='twice', fatol=1.0e-11_dp, maxiter=50)
    calls_both = jvp_calls
    call verify(both, 'combined image reuse')

    call check(base%success .and. retained%success .and. &
        candidate%success .and. both%success, 'all policies converge')
    call check(base%nit > 2 .and. both%nit > 2, &
        'more than two nonlinear steps exercise retained image compaction')
    call check(calls_retained < calls_base, 'retained cache saves Jacobian calls')
    call check(calls_candidate < calls_base, 'candidate cache saves Jacobian calls')
    call check(calls_both <= calls_retained .and. calls_both <= calls_candidate, &
        'both caches reuse the existing operator images')
    call check(maxval(abs(both%x - base%x)) < 2.0e-9_dp, &
        'combined cache agrees with uncached solution')

    ! outer_k=1 exercises a full history and FIFO replacement on every
    ! independent insertion. Both cached paths must handle the one-column
    ! image cache without assuming extra retained storage.
    call reset()
    small_history = root(residual, initial, method='krylov', inner_method='lgmres', &
        inner_m=preconditioner, jacobian_product=jvp, inner_maxiter=8, &
        outer_k=1, outer_image_prune_tol=1.0e-8_dp, &
        outer_image_reuse=.true., outer_candidate_image_reuse=.true., &
        inner_orthogonalization='selective', fatol=1.0e-11_dp, maxiter=50)
    call verify(small_history, 'one retained column')
    call check(small_history%success, 'one-column FIFO history converges')

    if (failures > 0) then
        print '(a,i0)', 'test_optimize_root_krylov_cache_compaction: FAIL ', failures
        stop 1
    end if
    print '(a)', 'test_optimize_root_krylov_cache_compaction: PASS'
    print '(a,4(i0,1x))', 'analytic Jv calls: ', &
        calls_base, calls_retained, calls_candidate, calls_both

contains
    subroutine reset()
        nfev_actual = 0
        jvp_calls = 0
    end subroutine reset

    subroutine check(valid, label)
        logical, intent(in) :: valid !! assertion outcome
        character(len=*), intent(in) :: label !! assertion name
        if (.not. valid) then
            failures = failures + 1
            print '(a)', 'FAIL: '//label
        end if
    end subroutine check

    subroutine verify(result, label)
        type(vector_root_result), intent(in) :: result !! nonlinear solver result
        character(len=*), intent(in) :: label !! reference label
        real(dp) :: value(n)

        call check(result%success, trim(label)//' converged')
        call check(result%nfev == nfev_actual, trim(label)//' callback accounting')
        call check(all(ieee_is_finite(result%x)), trim(label)//' finite solution')
        call model_residual(result%x, value)
        call check(maxval(abs(value)) < 2.0e-11_dp, trim(label)//' true residual')
        call check(maxval(abs(result%x - target)) < 1.0e-9_dp, &
            trim(label)//' correct solution')
    end subroutine verify

    pure subroutine model_residual(x, fvec)
        real(dp), intent(in) :: x(:) !! nonlinear variables
        real(dp), intent(out) :: fvec(:) !! nonlinear residual
        real(dp) :: delta(n)
        real(dp) :: projection

        delta = x - target
        projection = dot_product(coupling, delta) / real(n, dp)
        fvec = diagonal * (x*x - target*target) + 0.04_dp * projection * delta
    end subroutine model_residual

    subroutine residual(x, fvec)
        real(dp), intent(in) :: x(:) !! nonlinear variables
        real(dp), intent(out) :: fvec(:) !! nonlinear residual
        nfev_actual = nfev_actual + 1
        call model_residual(x, fvec)
    end subroutine residual

    subroutine jvp(x, vec, product, ok)
        real(dp), intent(in) :: x(:) !! nonlinear iterate
        real(dp), intent(in) :: vec(:) !! Krylov direction
        real(dp), intent(out) :: product(:) !! analytic Jacobian-vector action
        logical, intent(out) :: ok !! finite output
        real(dp) :: delta(n), projection, vprojection

        jvp_calls = jvp_calls + 1
        delta = x - target
        projection = dot_product(coupling, delta) / real(n, dp)
        vprojection = dot_product(coupling, vec) / real(n, dp)
        product = 2.0_dp * diagonal * x * vec + &
            0.04_dp * (projection * vec + vprojection * delta)
        ok = all(ieee_is_finite(product))
    end subroutine jvp

    subroutine preconditioner(vec, product, ok)
        real(dp), intent(in) :: vec(:) !! matrix-free direction
        real(dp), intent(out) :: product(:) !! approximate inverse action
        logical, intent(out) :: ok !! finite output

        product = vec / (2.0_dp * diagonal * initial)
        ok = all(ieee_is_finite(product))
    end subroutine preconditioner
end program test_optimize_root_krylov_cache_compaction
