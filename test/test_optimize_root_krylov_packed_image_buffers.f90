! SPDX-License-Identifier: MIT
! Exercise short and growing image histories with borrowed Arnoldi scratch.
program test_optimize_root_krylov_packed_image_buffers
    use scifort_kinds, only : dp
    use scifort_optimize_root, only : root, vector_root_result
    implicit none

    integer, parameter :: n = 8
    real(dp) :: target(n), x0(n)
    type(vector_root_result) :: reference, cached, finite_difference
    integer :: i, k, failures, jvp_calls, reference_calls, cached_calls
    logical :: auto_mode

    failures = 0
    do i = 1, n
        target(i) = 0.7_dp + 0.2_dp * real(i, dp)
        x0(i) = target(i) + 0.45_dp + 0.01_dp * real(i, dp)
    end do

    do k = 1, 3
        ! outer_k=1 exercises both the initial zero-column image array and
        ! FIFO replacement. Higher capacities exercise growing live images.
        jvp_calls = 0
        reference = root(residual, x0, method='krylov', &
            inner_method='lgmres', inner_maxiter=8, outer_k=2**(k-1), &
            outer_image_prune_auto=.true., jacobian_product=jvp, &
            fatol=1.0e-11_dp, line_search='none', maxiter=40)
        reference_calls = jvp_calls
        jvp_calls = 0
        cached = root(residual, x0, method='krylov', &
            inner_method='lgmres', inner_maxiter=8, outer_k=2**(k-1), &
            outer_image_prune_auto=.true., outer_image_reuse=.true., &
            outer_candidate_image_reuse=.true., &
            outer_image_basis_reuse=.true., jacobian_product=jvp, &
            fatol=1.0e-11_dp, line_search='none', maxiter=40)
        cached_calls = jvp_calls
        call check(reference%success .and. cached%success, 'analytic convergence')
        call check(reference%nit == cached%nit .and. &
            reference%nfev == cached%nfev, 'analytic work counts')
        call check(all(reference%x == cached%x), 'analytic solution bitwise identity')
        call check(cached_calls <= reference_calls, 'analytic Jv cache saves calls')
        call check(maxval(abs(cached%fun)) < 1.0e-11_dp, 'analytic residual')
    end do

    do k = 1, 3
        auto_mode = k == 3
        finite_difference = root(residual, x0, method='krylov', &
            inner_method='lgmres', inner_maxiter=8, outer_k=k, &
            outer_image_prune_tol=merge(0.0_dp, 1.0e-8_dp, auto_mode), &
            outer_image_prune_auto=auto_mode, outer_image_reuse=.true., &
            outer_image_basis_reuse=.true., &
            jacobian_diff_scheme='3-point', &
            fatol=1.0e-9_dp, line_search='none', maxiter=40)
        call check(finite_difference%success, 'finite-difference convergence')
        call check(maxval(abs(finite_difference%fun)) < 1.0e-9_dp, &
            'finite-difference residual')
    end do

    if (failures /= 0) then
        print '(a,i0)', 'test_optimize_root_krylov_packed_image_buffers: FAIL ', failures
        stop 1
    end if
    print '(a)', 'test_optimize_root_krylov_packed_image_buffers: PASS'
contains
    subroutine residual(x, f)
        real(dp), intent(in) :: x(:) !! iterate
        real(dp), intent(out) :: f(:) !! residual
        f = x*x - target*target + 0.05_dp * sum(x - target)
    end subroutine residual

    subroutine jvp(x, v, result_vec, ok)
        real(dp), intent(in) :: x(:) !! current iterate
        real(dp), intent(in) :: v(:) !! multiplication direction
        real(dp), intent(out) :: result_vec(:) !! Jacobian image
        logical, intent(out) :: ok !! validity flag
        jvp_calls = jvp_calls + 1
        result_vec = 2.0_dp*x*v + 0.05_dp * sum(v)
        ok = .true.
    end subroutine jvp

    subroutine check(condition, description)
        logical, intent(in) :: condition !! required outcome
        character(len=*), intent(in) :: description !! failure label
        if (.not. condition) then
            failures = failures + 1
            print '(a)', 'FAIL: '//description
        end if
    end subroutine check
end program test_optimize_root_krylov_packed_image_buffers
