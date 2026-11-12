! SPDX-License-Identifier: BSD-3-Clause
! Public LGMRES regressions with dependent/ill-conditioned augmentation data.
program test_optimize_root_krylov_rank_ls
    use scifort_kinds, only : dp
    use scifort_optimize_root, only : root, newton_krylov, vector_root_result
    implicit none

    real(dp), parameter :: x0(6) = [3.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp, 7.0_dp]
    character(len=9), parameter :: modes(3) = &
        [character(len=9) :: 'once', 'selective', 'twice']
    type(vector_root_result) :: res
    integer :: failures
    integer :: j
    integer :: jvp_count

    failures = 0
    do j = 1, size(modes)
        jvp_count = 0
        res = root(singular_residual, x0, method='krylov', &
            inner_method='lgmres', inner_maxiter=6, outer_k=4, &
            jacobian_product=singular_jvp, inner_orthogonalization=trim(modes(j)), &
            maxiter=25, fatol=1.0e-11_dp, line_search='none')
        call check(res%success .and. maxval(abs(res%fun)) < 1.0e-11_dp, &
            trim(modes(j))//' rank-deficient nonlinear convergence')
        call check(res%nit >= 3 .and. jvp_count > 6, &
            trim(modes(j))//' reuses LGMRES augmented vectors')
        call check(abs(res%x(1) + res%x(2) - 2.0_dp) < 1.0e-10_dp .and. &
            maxval(abs(res%x(3:6) - [3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp])) < 1.0e-10_dp, &
            trim(modes(j))//' correct rank-deficient root')
    end do

    ! A two-variable diagonal system with singular values 1 and 1e-10
    ! forces the small upper-triangular LGMRES system through the rank-
    ! revealing fallback. The normal-equations fallback drops the small
    ! direction and needs an extra nonlinear iteration; direct SVD resolves
    ! it within the available three termination checks.
    res = root(ill_residual, [0.0_dp, 0.0_dp], method='krylov', &
        inner_method='lgmres', inner_maxiter=2, outer_k=0, &
        jacobian_product=ill_jvp, line_search='none', &
        fatol=1.0e-8_dp, maxiter=3)
    call check(res%success .and. res%nfev <= 3, &
        'ill-conditioned LGMRES rank-revealing fallback converges promptly')
    call check(abs(res%x(1) - 1.0_dp) < 1.0e-7_dp .and. &
        abs(1.0e-10_dp * res%x(2) - 1.0_dp) < 1.0e-7_dp, &
        'ill-conditioned LGMRES preserves the small singular direction')

    res = newton_krylov(singular_residual, x0, &
        inner_method='lgmres', inner_maxiter=6, outer_k=4, &
        jacobian_product=singular_jvp, inner_orthogonalization='twice', &
        maxiter=25, f_tol=1.0e-11_dp, line_search='none')
    call check(res%success .and. maxval(abs(res%fun)) < 1.0e-11_dp, &
        'newton_krylov dependent augmented basis')

    if (failures > 0) then
        print '(a,i0)', 'test_optimize_root_krylov_rank_ls: FAIL ', failures
        stop 1
    end if
    print '(a)', 'test_optimize_root_krylov_rank_ls: PASS'

contains

    subroutine check(condition, label)
        logical, intent(in) :: condition !! expected outcome
        character(len=*), intent(in) :: label !! regression description

        if (.not. condition) then
            failures = failures + 1
            print '(a)', 'FAIL: '//label
        end if
    end subroutine check

    subroutine ill_residual(x, values)
        real(dp), intent(in) :: x(:) !! poorly scaled linear state
        real(dp), intent(out) :: values(:) !! diagonal residual with small singular value

        values = [x(1) - 1.0_dp, 1.0e-10_dp * x(2) - 1.0_dp]
    end subroutine ill_residual

    subroutine ill_jvp(x, vec, product, ok)
        real(dp), intent(in) :: x(:) !! nonlinear point
        real(dp), intent(in) :: vec(:) !! Jacobian direction
        real(dp), intent(out) :: product(:) !! Jacobian-vector product
        logical, intent(out) :: ok !! validity flag

        product = [vec(1), 1.0e-10_dp * vec(2)]
        ok = .true.
    end subroutine ill_jvp

    subroutine singular_residual(x, values)
        real(dp), intent(in) :: x(:) !! nonlinear state
        real(dp), intent(out) :: values(:) !! rank-deficient nonlinear residual

        values(1:2) = (x(1) + x(2))**2 - 4.0_dp
        values(3:6) = x(3:6)**2 - [9.0_dp, 16.0_dp, 25.0_dp, 36.0_dp]
    end subroutine singular_residual

    subroutine singular_jvp(x, vec, product, ok)
        real(dp), intent(in) :: x(:) !! current nonlinear point
        real(dp), intent(in) :: vec(:) !! matrix-free direction
        real(dp), intent(out) :: product(:) !! Jacobian-vector product
        logical, intent(out) :: ok !! true when computed successfully

        jvp_count = jvp_count + 1
        product(1:2) = 2.0_dp * (x(1) + x(2)) * (vec(1) + vec(2))
        product(3:6) = 2.0_dp * x(3:6) * vec(3:6)
        ok = .true.
    end subroutine singular_jvp

end program test_optimize_root_krylov_rank_ls
