! SPDX-License-Identifier: MIT
! Opt-in sampled homology-birth-count stopping for two symmetric wells.
program optimize_shgo_hgr_example
    use scifort_kinds, only : dp
    use scifort_optimize, only : shgo, shgo_result, shgo_hgr_stable
    implicit none
    real(dp) :: bounds(1,2)
    type(shgo_result) :: early, full

    bounds(1,:) = [-1.0_dp,1.0_dp]
    early = shgo(two_wells,bounds,iters=5,maxfun=1000,local_iter=2, &
        hgr_stability_iters=2,hgr_min_iters=3,hgr_f_tol=1.0e-8_dp)
    full = shgo(two_wells,bounds,iters=5,maxfun=1000,local_iter=2)
    print '(a,i0,a,i0)', 'HGR result: status=',early%status,', rounds=',early%nit
    print '(a,i0,a,i0)', 'graph minima=',early%hgr_estimate, &
        ', unchanged comparisons=',early%hgr_streak
    print '(a,i0,a,i0)', 'evaluations with HGR=',early%nfev, &
        ', without HGR=',full%nfev
    print '(a,f10.6,a,es12.4)', 'best x=',early%x(1),', f=',early%fun
    if (early%status == shgo_hgr_stable) then
        print '(a)', 'Terminated by heuristic HGR stability, not certified optimality.'
    end if
contains
    function two_wells(x) result(f)
        real(dp), intent(in) :: x(:) !! one-dimensional bounded coordinate
        real(dp) :: f !! two equal minima at -0.5 and +0.5
        f = (x(1)**2-0.25_dp)**2
    end function two_wells
end program optimize_shgo_hgr_example
