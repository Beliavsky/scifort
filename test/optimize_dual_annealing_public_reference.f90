! SPDX-License-Identifier: MIT
! These contract tests link against a mock optimizer, not the real L-BFGS-B.
program test_optimize_dual_annealing_public_reference
    use scifort_kinds, only : dp
    use scifort_optimize_dual_annealing, only : da_maxiter, da_invalid, dual_annealing_result
    use scifort_optimize_dual_annealing_minimize, only : dual_annealing
    implicit none
    type(dual_annealing_result) :: res
    integer :: checks, i
    real(dp) :: bounds(2,2)
    character(len=12) :: methods(2)

    checks = 0
    bounds(:,1) = -4.0_dp
    bounds(:,2) = 4.0_dp
    methods = [character(len=12) :: 'l-bfgs-b', 'powell']
    do i = 1, size(methods)
        res = dual_annealing(objective, bounds, maxiter=5, seed=8, method=trim(methods(i)), &
            local_interval=1, local_maxiter=80)
        call must(res%status == da_maxiter .and. res%success, 'supported local method')
        call must(res%nlocal > 0, 'local optimizer invoked')
        call must(res%fun < 1.0e-8_dp, 'local optimizer finds quadratic minimum')
        call must(res%nfev > 21, 'include local objective evaluations')
    end do
    res = dual_annealing(objective, bounds, maxiter=3, seed=8, no_local_search=.true.)
    call must(res%nlocal == 0 .and. res%status == da_maxiter, 'public no-local-search option')
    res = dual_annealing(objective, bounds, method='invalid')
    call must(res%status == da_invalid, 'unknown local method rejected')
    res = dual_annealing(objective, bounds, local_maxiter=-2)
    call must(res%status == da_invalid, 'negative local iterations rejected')
    print '(a,1x,i0)', 'DUAL_ANNEALING_PUBLIC_CHECKS', checks
contains
    function objective(x) result(f)
        real(dp), intent(in) :: x(:)
        real(dp) :: f
        f = (x(1)-0.5_dp)**2 + (x(2)+0.75_dp)**2
    end function objective
    subroutine must(condition, label)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: label
        if (.not. condition) then
            print '(a)', 'FAIL: '//label
            error stop 1
        end if
        checks = checks + 1
    end subroutine must
end program test_optimize_dual_annealing_public_reference
