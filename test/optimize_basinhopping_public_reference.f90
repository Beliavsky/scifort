! SPDX-License-Identifier: MIT
program test_optimize_basinhopping_public_reference
    use scifort_kinds, only : dp
    use scifort_optimize_basinhopping_minimize, only : basinhopping
    use scifort_optimize_basinhopping, only : basinhopping_result, bh_completed, bh_invalid
    implicit none
    type(basinhopping_result) :: res
    integer :: checks, i
    character(len=12) :: methods(3)

    checks = 0
    methods = [character(len=12) :: 'nelder-mead', 'powell', 'bfgs']
    do i = 1, size(methods)
        res = basinhopping(quadratic, [3.0_dp, -3.0_dp], niter=0, &
            method=trim(methods(i)), seed=58)
        call must(res%success .and. res%status == bh_completed, 'local method available')
        call must(res%nlocal == 1 .and. res%nfev > 0, 'initial local solver called')
        call must(res%fun < 1.0e-5_dp, 'local minimum agrees with known solution')
    end do
    res = basinhopping(quadratic, [3.0_dp, -3.0_dp], niter=5, seed=114, &
        temperature=0.0_dp, local_maxiter=1000)
    call must(res%success .and. res%nit == 5, 'basin hopping calls SciFort local minimizer')
    call must(res%fun < 1.0e-5_dp, 'global best and local result agree')
    res = basinhopping(quadratic, [3.0_dp, -3.0_dp], method='not-a-method')
    call must(res%status == bh_invalid, 'invalid local method rejected')
    res = basinhopping(quadratic, [3.0_dp, -3.0_dp], local_maxiter=-1)
    call must(res%status == bh_invalid, 'invalid local iteration budget rejected')
    print '(a,1x,i0)', 'BASINHOPPING_PUBLIC_CHECKS', checks
contains
    function quadratic(x) result(fx)
        real(dp), intent(in) :: x(:)
        real(dp) :: fx
        fx = (x(1) - 1.25_dp)**2 + (x(2) + 0.75_dp)**2
    end function quadratic
    subroutine must(condition, label)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: label
        if (.not.condition) then
            print '(a)', 'FAIL: '//label
            error stop 1
        end if
        checks = checks + 1
    end subroutine must
end program test_optimize_basinhopping_public_reference
