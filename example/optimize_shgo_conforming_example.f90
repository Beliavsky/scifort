! SPDX-License-Identifier: MIT
! Two-dimensional conforming adaptive SHGO using Nelder-Mead polishing.
program optimize_shgo_conforming_example
    use scifort_kinds, only : dp
    use scifort_optimize, only : shgo, shgo_result, shgo_completed
    implicit none
    type(shgo_result) :: result
    real(dp) :: bounds(2,2)

    bounds(:,1) = [-6.0_dp,-6.0_dp]
    bounds(:,2) = [ 6.0_dp, 6.0_dp]
    result = shgo(himmelblau,bounds,iters=4,maxfun=5000, &
        adaptive=.true.,conforming=.true.,adaptive_fraction=0.25_dp, &
        method='nelder-mead',local_iter=4)
    print '(a,i0)', 'status: ',result%status
    print '(a,2f13.7)', 'best point: ',result%x
    print '(a,es13.5)', 'best objective: ',result%fun
    print '(a,i0)', 'evaluations: ',result%nfev
    print '(a,i0)', 'adaptive leaves: ',result%nleaf
    print '(a,i0)', 'conforming triangles: ',result%ntriangles
    print '(a,l1)', 'triangulation validated: ',result%conforming_valid
    if (result%status /= shgo_completed) print '(a)',trim(result%message)
contains
    function himmelblau(x) result(value)
        real(dp), intent(in) :: x(:) !! two-dimensional search coordinate
        real(dp) :: value !! Himmelblau multimodal objective
        value = (x(1)*x(1)+x(2)-11.0_dp)**2 &
            +(x(1)+x(2)*x(2)-7.0_dp)**2
    end function himmelblau
end program optimize_shgo_conforming_example
