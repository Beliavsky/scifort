! SPDX-License-Identifier: MIT
! Deterministic global search without a local optimizer or random generator.
program optimize_direct_example
    use scifort_kinds, only : dp
    use scifort_optimize_direct, only : direct, direct_result
    implicit none
    type(direct_result) :: res
    real(dp) :: domain(2,2)
    domain(:,1) = [-5.0_dp, -5.0_dp]
    domain(:,2) = [ 5.0_dp,  5.0_dp]
    res = direct(rastrigin, domain, maxfun=1200, maxiter=100)
    if (.not. allocated(res%x)) error stop 'DIRECT returned no solution'
    if (res%fun > 1.0e-6_dp) error stop 'DIRECT failed Rastrigin example'
    print '(a,2(1x,f12.7))', 'best point:', res%x
    print '(a,es14.6)', 'objective: ', res%fun
    print '(a,i0,a,i0)', 'evaluations: ', res%nfev, ', iterations: ', res%nit
contains
    function rastrigin(x) result(f)
        real(dp), intent(in) :: x(:) !! real argument, length two
        real(dp) :: f !! nonnegative Rastrigin objective
        real(dp) :: pi
        pi = acos(-1.0_dp)
        f = 20.0_dp + sum(x*x - 10.0_dp*cos(2.0_dp*pi*x))
    end function rastrigin
end program optimize_direct_example
