! SPDX-License-Identifier: MIT
! Exercise the opt-in full-history LGMRES QR-downdate path through public APIs.
program test_optimize_root_krylov_fifo_qr
    use scifort_kinds, only : dp
    use scifort_optimize_root, only : root, newton_krylov, vector_root_result
    implicit none

    integer, parameter :: n = 8
    real(dp) :: target(n)
    real(dp) :: x0(n)
    real(dp) :: check_f(n)
    type(vector_root_result) :: old_result
    type(vector_root_result) :: fast_result
    type(vector_root_result) :: wrapper_result
    type(vector_root_result) :: bad_result
    integer :: failures
    integer :: i
    integer :: calls

    failures = 0
    calls = 0
    do i = 1, n
        target(i) = 0.7_dp + 0.2_dp * real(i, dp)
        x0(i) = target(i) + 0.4_dp + 0.05_dp * real(i, dp)
    end do
    old_result = root(residual, x0, method='krylov', inner_method='lgmres', &
        inner_maxiter=8, outer_k=2, outer_prune_tol=1.0e-8_dp, &
        jacobian_product=jvp, fatol=1.0e-12_dp, line_search='none', &
        maxiter=40)
    fast_result = root(residual, x0, method='krylov', inner_method='lgmres', &
        inner_maxiter=8, outer_k=2, outer_prune_tol=1.0e-8_dp, &
        outer_fifo_prune_fast=.true., jacobian_product=jvp, &
        fatol=1.0e-12_dp, line_search='none', maxiter=40)
    wrapper_result = newton_krylov(residual, x0, inner_method='lgmres', &
        inner_maxiter=8, outer_k=2, outer_prune_tol=1.0e-8_dp, &
        outer_fifo_prune_fast=.true., jacobian_product=jvp, &
        f_tol=1.0e-12_dp, line_search='none', maxiter=40)
    call check(old_result%success, 'legacy public root converges')
    call check(fast_result%success, 'fast FIFO public root converges')
    call check(wrapper_result%success, 'newton_krylov forwards fast option')
    call check(fast_result%nit >= 4, 'multiple nonlinear steps exercise FIFO')
    call model(fast_result%x, check_f)
    call check(maxval(abs(check_f)) < 1.0e-11_dp, 'verified fast root residual')
    call check(maxval(abs(fast_result%x - old_result%x)) < 1.0e-9_dp, &
        'fast and baseline converge to same root')
    call check(maxval(abs(wrapper_result%x - fast_result%x)) < 1.0e-9_dp, &
        'public wrappers agree')

    ! All invalid combinations must fail before calling the residual.
    calls = 0
    bad_result = root(residual, x0, method='hybr', outer_fifo_prune_fast=.true.)
    call check(.not. bad_result%success .and. calls == 0, 'non-Krylov method rejected')
    bad_result = root(residual, x0, method='krylov', outer_fifo_prune_fast=.true.)
    call check(.not. bad_result%success .and. calls == 0, 'missing pruning tolerance rejected')
    bad_result = root(residual, x0, method='krylov', inner_method='gmres', &
        outer_fifo_prune_fast=.true., outer_prune_tol=1.0e-8_dp)
    call check(.not. bad_result%success .and. calls == 0, 'non-LGMRES method rejected')
    bad_result = root(residual, x0, method='krylov', inner_method='lgmres', &
        outer_k=1, outer_fifo_prune_fast=.true., outer_prune_tol=1.0e-8_dp)
    call check(.not. bad_result%success .and. calls == 0, 'single-vector FIFO rejected')

    if (failures /= 0) then
        print '(a,i0)', 'test_optimize_root_krylov_fifo_qr: FAIL ', failures
        stop 1
    end if
    print '(a)', 'test_optimize_root_krylov_fifo_qr: PASS'
contains
    subroutine model(x, f)
        real(dp), intent(in) :: x(:) !! trial root
        real(dp), intent(out) :: f(:) !! residual
        f = x*x - target*target + 0.05_dp * sum(x - target)
    end subroutine model

    subroutine residual(x, f)
        real(dp), intent(in) :: x(:) !! trial root
        real(dp), intent(out) :: f(:) !! residual
        calls = calls + 1
        call model(x, f)
    end subroutine residual

    subroutine jvp(x, vec, result_vec, ok)
        real(dp), intent(in) :: x(:) !! linearization point
        real(dp), intent(in) :: vec(:) !! tangent direction
        real(dp), intent(out) :: result_vec(:) !! Jacobian-vector product
        logical, intent(out) :: ok !! valid product
        result_vec = 2.0_dp*x*vec + 0.05_dp * sum(vec)
        ok = .true.
    end subroutine jvp

    subroutine check(condition, description)
        logical, intent(in) :: condition !! regression requirement
        character(len=*), intent(in) :: description !! failure label
        if (.not. condition) then
            failures = failures + 1
            print '(a)', 'FAIL: '//description
        end if
    end subroutine check
end program test_optimize_root_krylov_fifo_qr
