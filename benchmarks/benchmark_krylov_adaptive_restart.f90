! SPDX-License-Identifier: MIT
! Deterministic fixed-versus-adaptive restarted LGMRES work comparison.
!
! Build against SciFort with FPM or the offline root test dependencies
! described in docs/validation.md. Output is CSV on standard output.
program benchmark_krylov_adaptive_restart
    use scifort_kinds, only : dp
    use scifort_optimize_root, only : root, vector_root_result
    implicit none

    integer, parameter :: n = 12
    real(dp) :: diagonal(n), target(n), x0(n)
    type(vector_root_result) :: result
    integer :: i
    integer :: problem
    integer :: mode
    integer :: jvp_count
    integer :: residual_count
    character(len=12) :: problem_label
    character(len=12) :: mode_label

    x0 = 0.0_dp
    target = 1.0_dp
    do i = 1, n
        diagonal(i) = 2.0_dp**(0.6_dp * real(i - 1, dp))
    end do
    print '(a)', 'problem,mode,success,nit,nfev,njvp,max_residual'
    do problem = 1, 2
        if (problem == 1) then
            problem_label = 'skew2'
        else
            problem_label = 'scaled12'
        end if
        do mode = 1, 3
            residual_count = 0
            jvp_count = 0
            select case (mode)
            case (1)
                mode_label = 'fixed'
                if (problem == 1) then
                    result = root(residual, x0(:2), method='krylov', &
                        inner_method='lgmres', inner_restart=1, inner_maxiter=8, &
                        outer_k=0, jacobian_product=jvp, line_search='none', &
                        fatol=1.0e-12_dp, maxiter=5)
                else
                    result = root(residual, x0, method='krylov', &
                        inner_method='lgmres', inner_restart=1, inner_maxiter=36, &
                        outer_k=0, jacobian_product=jvp, line_search='none', &
                        fatol=1.0e-8_dp, maxiter=5)
                end if
            case (2:3)
                mode_label = 'adaptive'
                if (mode == 3) mode_label = 'capped'
                if (problem == 1) then
                    result = root(residual, x0(:2), method='krylov', &
                        inner_method='lgmres', inner_restart=1, inner_maxiter=8, &
                        inner_restart_adaptive=.true., inner_restart_max=merge(1, 2, mode == 3), &
                        outer_k=0, jacobian_product=jvp, line_search='none', &
                        fatol=1.0e-12_dp, maxiter=5)
                else
                    result = root(residual, x0, method='krylov', &
                        inner_method='lgmres', inner_restart=1, inner_maxiter=36, &
                        inner_restart_adaptive=.true., inner_restart_max=merge(1, 4, mode == 3), &
                        outer_k=0, jacobian_product=jvp, line_search='none', &
                        fatol=1.0e-8_dp, maxiter=5)
                end if
            end select
            if (result%nfev /= residual_count) error stop 'evaluation count mismatch'
            print '(a,",",a,",",l1,",",i0,",",i0,",",i0,",",es16.8)', &
                trim(problem_label), trim(mode_label), result%success, result%nit, &
                result%nfev, jvp_count, maxval(abs(result%fun))
        end do
    end do
contains
    subroutine residual(x, fvec)
        real(dp), intent(in) :: x(:) !! current nonlinear point
        real(dp), intent(out) :: fvec(:) !! residual vector
        residual_count = residual_count + 1
        if (problem == 1) then
            fvec(1) = -(x(2) - 2.0_dp)
            fvec(2) = x(1) - 1.0_dp
        else
            fvec = diagonal * (x - target)
        end if
    end subroutine

    subroutine jvp(x, vec, product, ok)
        real(dp), intent(in) :: x(:) !! current linearization point
        real(dp), intent(in) :: vec(:) !! direction to multiply
        real(dp), intent(out) :: product(:) !! analytic Jacobian-vector product
        logical, intent(out) :: ok !! true for valid product
        jvp_count = jvp_count + 1
        if (problem == 1) then
            product(1) = -vec(2)
            product(2) = vec(1)
        else
            product = diagonal * vec
        end if
        ok = .true.
    end subroutine
end program benchmark_krylov_adaptive_restart
