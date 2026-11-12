! SPDX-License-Identifier: MIT
! Run with: fpm run --example optimize_shgo_conforming_3d_example
program optimize_shgo_conforming_3d_example
    use scifort_kinds, only : dp
    use scifort_optimize_shgo, only : shgo_core, shgo_result
    implicit none
    type(shgo_result) :: res
    real(dp) :: bounds(3,2)

    bounds(:,1) = [-1.0_dp,-1.0_dp,-1.0_dp]
    bounds(:,2) = [ 1.0_dp, 1.0_dp, 1.0_dp]
    res = shgo_core(shifted_sphere,bounds,iters=5,maxfun=10000, &
        adaptive=.true.,conforming=.true.,adaptive_fraction=0.1_dp, &
        local_iter=0)
    print '(a,l1)', 'mesh valid: ',res%conforming_valid
    print '(a,4(i0,1x))', 'vertices, leaves, splits, tetrahedra: ', &
        res%nmesh,res%nleaf,res%nrefined,res%ntetrahedra
    print '(a,es14.6)', 'sampled objective: ',res%fun
    print '(a,3(f10.6,1x))', 'best point: ',res%x
contains
    function shifted_sphere(x) result(value)
        real(dp), intent(in) :: x(:) !! bounded three-dimensional vector
        real(dp) :: value !! squared distance to known optimum
        value = sum((x-[0.3_dp,-0.2_dp,0.15_dp])**2)
    end function shifted_sphere
end program optimize_shgo_conforming_3d_example
