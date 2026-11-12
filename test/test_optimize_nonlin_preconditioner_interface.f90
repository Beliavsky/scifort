! SPDX-License-Identifier: MIT
! Regression for calling an optional inverse-Jacobian callback from the
! internal Krylov preconditioner, with an explicit procedure interface.
program test_optimize_nonlin_preconditioner_interface
    use scifort_kinds, only : dp
    use scifort_optimize_nonlin, only : solve_krylov_root
    implicit none

    real(dp), parameter :: diagonal(4) = [1.0_dp, 10.0_dp, 100.0_dp, 1000.0_dp]
    real(dp), parameter :: rhs(4) = [2.0_dp, -3.0_dp, 4.0_dp, -5.0_dp]
    real(dp) :: x(4), residual(4)
    integer :: nfev, nit, status, callback_calls
    logical :: success
    character(len=160) :: message

    callback_calls = 0
    call solve_krylov_root(linear_system, 0.0_dp * rhs, x, residual, nfev, nit, &
        status, success, message, inner_m=inverse_diagonal, inner_method='gmres', &
        inner_maxiter=1, outer_k=0, rdiff=1.0e-6_dp, line_search='none', &
        fatol=1.0e-7_dp)
    if (.not. success .or. status /= 1) error stop 'preconditioned Krylov solve'
    if (maxval(abs(x-rhs/diagonal)) > 1.0e-7_dp) error stop 'Krylov solution'
    if (callback_calls <= 0) error stop 'optional operator was not invoked'

    callback_calls = 0
    call solve_krylov_root(linear_system, 0.0_dp * rhs, x, residual, nfev, nit, &
        status, success, message, inner_m=failing_operator, inner_method='gmres', &
        inner_maxiter=1, outer_k=0, line_search='none', fatol=1.0e-7_dp)
    if (success .or. callback_calls <= 0) error stop 'failed callback was ignored'
    if (index(message, 'inversion failed') <= 0) error stop 'callback failure message'

    call solve_krylov_root(linear_system, 0.0_dp * rhs, x, residual, nfev, nit, &
        status, success, message, inner_method='gmres', inner_maxiter=4, &
        outer_k=0, rdiff=1.0e-6_dp, line_search='none', fatol=1.0e-7_dp)
    if (.not. success) error stop 'identity fallback Krylov solve'
    if (maxval(abs(x-rhs/diagonal)) > 1.0e-7_dp) error stop 'identity solution'

    print '(a)', 'test_optimize_nonlin_preconditioner_interface passed'
contains
    subroutine linear_system(values, fvec)
        real(dp), intent(in) :: values(:) !! current nonlinear iterate
        real(dp), intent(out) :: fvec(:) !! linear residual
        fvec = diagonal * values - rhs
    end subroutine linear_system

    subroutine inverse_diagonal(vec, product, ok)
        real(dp), intent(in) :: vec(:) !! vector to precondition
        real(dp), intent(out) :: product(:) !! inverse diagonal product
        logical, intent(out) :: ok !! true if the operator succeeded
        callback_calls = callback_calls + 1
        product = vec / diagonal
        ok = .true.
    end subroutine inverse_diagonal

    subroutine failing_operator(vec, product, ok)
        real(dp), intent(in) :: vec(:) !! vector to precondition
        real(dp), intent(out) :: product(:) !! zeroed output on failure
        logical, intent(out) :: ok !! false, intentionally indicating failure
        callback_calls = callback_calls + 1
        product = 0.0_dp * vec
        ok = .false.
    end subroutine failing_operator
end program test_optimize_nonlin_preconditioner_interface
