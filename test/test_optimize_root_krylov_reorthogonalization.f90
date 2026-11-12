! SPDX-License-Identifier: BSD-3-Clause
! Check opt-in GMRES/LGMRES modified Gram-Schmidt reorthogonalization.
program test_optimize_root_krylov_reorthogonalization
    use scifort_kinds, only : dp
    use scifort_optimize_root, only : root, newton_krylov, vector_root_result, &
        vector_root_linear_operator
    implicit none

    integer, parameter :: n = 48
    integer, parameter :: inner_dim = 36
    character(len=9), parameter :: modes(3) = &
        [character(len=9) :: 'once', 'selective', 'twice']
    real(dp) :: diagonal(n)
    real(dp) :: rhs(n)
    real(dp) :: saved_basis(n, inner_dim)
    real(dp) :: overlap(3)
    real(dp) :: solution_once(n)
    real(dp) :: easy_target(4)
    real(dp) :: easy_diag(4)
    type(vector_root_result) :: res
    integer :: failures
    integer :: calls
    integer :: calls_by_mode(3)
    integer :: quadratic_calls
    integer :: quadratic_counts(3)
    integer :: i
    integer :: j
    integer :: mode

    failures = 0
    calls = 0
    do i = 1, n
        diagonal(i) = 10.0_dp**(-12.0_dp * real(i-1, dp) / real(n-1, dp))
        rhs(i) = sin(0.41_dp * real(i, dp))
    end do
    easy_diag = [1.0_dp, 2.0_dp, 4.0_dp, 8.0_dp]
    easy_target = [2.0_dp, -1.0_dp, 0.5_dp, 3.0_dp]

    ! A diagonally ill-conditioned matrix leads to nearly dependent Arnoldi
    ! vectors. Record the directions received by the exact J*v callback:
    ! each Arnoldi direction must remain orthogonal to earlier directions.
    ! We deliberately stop after one nonlinear step so every run starts with
    ! the same normalized right-hand side and builds the same number of
    ! inner Krylov vectors.
    do mode = 1, size(modes)
        calls = 0
        saved_basis = 0.0_dp
        res = root(ill_residual, spread(0.0_dp, 1, n), method='krylov', &
            inner_method='gmres', inner_maxiter=inner_dim, maxiter=1, &
            inner_orthogonalization=trim(modes(mode)), &
            jacobian_product=ill_jvp, line_search='none', fatol=1.0e-14_dp)
        calls_by_mode(mode) = calls
        overlap(mode) = max_basis_overlap(min(calls, inner_dim))
        call check(res%nfev == 2 .and. calls >= inner_dim, &
            trim(modes(mode))//' finite-difference-free Arnoldi work')
        call check(all(abs(res%x) < huge(1.0_dp)) .and. &
            all(abs(res%fun) < huge(1.0_dp)), &
            trim(modes(mode))//' finite GMRES iterate')
        if (mode == 1) solution_once = res%x
    end do
    call check(all(calls_by_mode == calls_by_mode(1)), &
        'reorthogonalization adds no Jacobian products')
    call check(overlap(2) < 1.0e-11_dp .and. overlap(3) < 1.0e-11_dp, &
        'DGKS and two-pass modes preserve Arnoldi orthogonality')
    call check(overlap(2) <= overlap(1) + 1.0e-13_dp .and. &
        overlap(3) <= overlap(1) + 1.0e-13_dp, &
        'reorthogonalization does not degrade the basis')

    ! The default must agree with explicitly requesting the original
    ! one-pass Arnoldi projection, including unchanged solver work counts.
    calls = 0
    res = root(ill_residual, spread(0.0_dp, 1, n), method='krylov', &
        inner_method='gmres', inner_maxiter=inner_dim, maxiter=1, &
        jacobian_product=ill_jvp, line_search='none', fatol=1.0e-14_dp)
    call check(calls == calls_by_mode(1) .and. res%nfev == 2 .and. &
        all(res%x == solution_once), 'implicit one-pass default is unchanged')

    ! Exercise the public LGMRES, preconditioner, and newton_krylov paths.
    do mode = 1, size(modes)
        res = root(easy_residual, [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], &
            method='krylov', inner_method='lgmres', outer_k=3, &
            inner_maxiter=12, inner_orthogonalization=trim(modes(mode)), &
            inner_m=easy_preconditioner, jacobian_product=easy_jvp, &
            line_search='none', fatol=1.0e-11_dp)
        call check(res%success .and. maxval(abs(res%x-easy_target)) < 1.0e-10_dp, &
            trim(modes(mode))//' preconditioned LGMRES convergence')
    end do
    res = newton_krylov(easy_residual, [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], &
        inner_method='gmres', inner_orthogonalization='TWICE', &
        jacobian_product=easy_jvp, line_search='none', f_tol=1.0e-11_dp)
    call check(res%success .and. maxval(abs(res%x-easy_target)) < 1.0e-10_dp, &
        'newton_krylov forwards case-insensitive orthogonalization option')

    ! Force several nonlinear iterations so LGMRES actually consumes and
    ! updates its retained outer augmentation vectors. A nonlinear diagonal
    ! quadratic system also tests that the Jacobian changes between steps.
    do mode = 1, size(modes)
        quadratic_calls = 0
        res = root(quadratic_residual, [4.0_dp, 4.0_dp, 4.0_dp, 4.0_dp], &
            method='krylov', inner_method='lgmres', inner_maxiter=4, &
            outer_k=3, inner_orthogonalization=trim(modes(mode)), &
            jacobian_product=quadratic_jvp, maxiter=25, &
            fatol=1.0e-10_dp, line_search='none')
        quadratic_counts(mode) = quadratic_calls
        call check(res%success .and. res%nit > 3 .and. &
            maxval(abs(res%x - [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp])) < 1.0e-10_dp, &
            trim(modes(mode))//' LGMRES augmented nonlinear convergence')
    end do
    call check(all(quadratic_counts == quadratic_counts(1)), &
        'reorthogonalization preserves augmented LGMRES Jacobian work')

    ! Reject unsupported combinations without evaluating the residual.
    res = root(easy_residual, [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], &
        method='krylov', inner_orthogonalization='third')
    call check(.not. res%success .and. res%nfev == 0 .and. &
        index(res%message, 'once, selective, or twice') /= 0, &
        'invalid projection option rejected')
    res = root(easy_residual, [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], &
        method='krylov', inner_method='tfqmr', inner_orthogonalization='twice')
    call check(.not. res%success .and. res%nfev == 0 .and. &
        index(res%message, 'requires built-in gmres or lgmres') /= 0, &
        'non-Arnoldi inner method rejects option')
    res = root(easy_residual, [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], &
        method='krylov', inner_solver=identity_solver, &
        inner_orthogonalization='selective')
    call check(.not. res%success .and. res%nfev == 0 .and. &
        index(res%message, 'requires built-in gmres or lgmres') /= 0, &
        'custom inner solver rejects option')
    res = root(easy_residual, [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], &
        method='broyden1', inner_orthogonalization='twice')
    call check(.not. res%success .and. res%nfev == 0 .and. &
        index(res%message, 'requires method=krylov') /= 0, &
        'non-Krylov root method rejects option')

    if (failures /= 0) then
        print '(a,i0)', 'test_optimize_root_krylov_reorthogonalization: FAIL ', failures
        stop 1
    end if
    print '(a)', 'test_optimize_root_krylov_reorthogonalization: PASS'

contains

    subroutine check(condition, description)
        logical, intent(in) :: condition !! expected condition
        character(len=*), intent(in) :: description !! case name

        if (.not. condition) then
            failures = failures + 1
            print '(a)', 'FAIL: '//description
        end if
    end subroutine check

    function max_basis_overlap(num_columns) result(value)
        integer, intent(in) :: num_columns !! recorded Krylov directions
        real(dp) :: value !! maximum absolute off-diagonal Gram matrix entry
        integer :: col
        integer :: prev

        value = 0.0_dp
        do col = 2, num_columns
            do prev = 1, col - 1
                value = max(value, abs(dot_product(saved_basis(:, col), &
                    saved_basis(:, prev))))
            end do
        end do
    end function max_basis_overlap

    subroutine ill_residual(x, f)
        real(dp), intent(in) :: x(:) !! nonlinear point
        real(dp), intent(out) :: f(:) !! diagonal-system residual

        f = diagonal * x - rhs
    end subroutine ill_residual

    subroutine ill_jvp(x, vec, product, ok)
        real(dp), intent(in) :: x(:) !! nonlinear point
        real(dp), intent(in) :: vec(:) !! Arnoldi direction
        real(dp), intent(out) :: product(:) !! J times direction
        logical, intent(out) :: ok !! callback success

        calls = calls + 1
        if (calls <= inner_dim) saved_basis(:, calls) = vec
        product = diagonal * vec
        ok = .true.
    end subroutine ill_jvp

    subroutine easy_residual(x, f)
        real(dp), intent(in) :: x(:) !! nonlinear point
        real(dp), intent(out) :: f(:) !! linear residual

        f = easy_diag * (x-easy_target)
    end subroutine easy_residual

    subroutine easy_jvp(x, vec, product, ok)
        real(dp), intent(in) :: x(:) !! nonlinear point
        real(dp), intent(in) :: vec(:) !! Jacobian-vector direction
        real(dp), intent(out) :: product(:) !! J times direction
        logical, intent(out) :: ok !! callback success

        product = easy_diag * vec
        ok = .true.
    end subroutine easy_jvp

    subroutine easy_preconditioner(vec, product, ok)
        real(dp), intent(in) :: vec(:) !! Krylov direction
        real(dp), intent(out) :: product(:) !! approximate inverse action
        logical, intent(out) :: ok !! callback success

        product = vec / easy_diag
        ok = .true.
    end subroutine easy_preconditioner

    subroutine quadratic_residual(x, f)
        real(dp), intent(in) :: x(:) !! nonlinear iterate
        real(dp), intent(out) :: f(:) !! quadratic system residual

        f = x * x - [1.0_dp, 4.0_dp, 9.0_dp, 16.0_dp]
    end subroutine quadratic_residual

    subroutine quadratic_jvp(x, vec, product, ok)
        real(dp), intent(in) :: x(:) !! current nonlinear iterate
        real(dp), intent(in) :: vec(:) !! multiplication direction
        real(dp), intent(out) :: product(:) !! exact Jacobian-vector product
        logical, intent(out) :: ok !! callback success

        quadratic_calls = quadratic_calls + 1
        product = 2.0_dp * x * vec
        ok = .true.
    end subroutine quadratic_jvp

    subroutine identity_solver(matvec, preconditioner, rhs_vec, rtol, maxiter, solution, ok)
        procedure(vector_root_linear_operator) :: matvec !! Jacobian action
        procedure(vector_root_linear_operator) :: preconditioner !! inverse approximation
        real(dp), intent(in) :: rhs_vec(:) !! right-hand side
        real(dp), intent(in) :: rtol !! inner relative tolerance
        integer, intent(in) :: maxiter !! maximum inner iterations
        real(dp), intent(out) :: solution(:) !! returned step
        logical, intent(out) :: ok !! callback success

        solution = rhs_vec
        ok = .true.
    end subroutine identity_solver

end program test_optimize_root_krylov_reorthogonalization
