! SPDX-License-Identifier: MIT
! Example of selective dyadic refinement in the public SciFort SHGO API.
program optimize_shgo_adaptive_example
    use scifort_kinds, only : dp
    use scifort_optimize, only : shgo, shgo_result
    implicit none
    type(shgo_result) :: adaptive, uniform
    real(dp) :: bounds(2,2)

    bounds(:,1) = -6.0_dp
    bounds(:,2) = 6.0_dp
    ! Comparing the same number of rounds, with bounded local polishing.
    adaptive = shgo(himmelblau, bounds, iters=4, maxfun=4000, &
        adaptive=.true., adaptive_fraction=0.25_dp, local_iter=3, &
        method='nelder-mead')
    uniform = shgo(himmelblau, bounds, iters=4, maxfun=4000, &
        local_iter=3, method='nelder-mead')
    print '(a,l1,2(a,i0))', 'adaptive success=',adaptive%success, &
        ', nfev=',adaptive%nfev,', nleaf=',adaptive%nleaf
    print '(a,2f14.7,a,es15.7)', 'adaptive best=',adaptive%x, &
        ', f=',adaptive%fun
    print '(a,l1,a,i0)', 'uniform success=',uniform%success, &
        ', nfev=',uniform%nfev
    print '(a,2f14.7,a,es15.7)', 'uniform best =',uniform%x, &
        ', f=',uniform%fun
contains
    function himmelblau(x) result(f)
        real(dp), intent(in) :: x(:) !! two coordinates
        real(dp) :: f !! Himmelblau objective
        f = (x(1)**2+x(2)-11.0_dp)**2 + (x(1)+x(2)**2-7.0_dp)**2
    end function himmelblau
end program optimize_shgo_adaptive_example
