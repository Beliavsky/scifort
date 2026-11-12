! SPDX-License-Identifier: MIT
! Independent SciPy objective-quality comparison for HGR-terminated searches.
program optimize_shgo_hgr_driver
    use scifort_kinds, only : dp
    use scifort_optimize_shgo, only : shgo_core, shgo_result
    implicit none
    integer :: case_id, dim, j, streak_requirement
    real(dp) :: bounds(2,2)
    type(shgo_result) :: result
    bounds(:,1) = -1.0_dp
    bounds(:,2) = 1.0_dp
    do case_id = 1, 4
        dim = 1
        if (case_id >= 3) dim = 2
        streak_requirement = 2
        if (case_id == 2) streak_requirement = 1
        result = shgo_core(objective,bounds(1:dim,:),iters=5, &
            hgr_stability_iters=streak_requirement,hgr_f_tol=0.0_dp,local_iter=0)
        write (*,'(5(i0,1x),es24.16,1x,3(i0,1x))',advance='no') &
            case_id,dim,result%status,result%nit,result%nfev,result%fun, &
            result%hgr_estimate,result%hgr_streak,result%nmesh
        do j = 1, dim
            write (*,'(es24.16,1x)',advance='no') result%x(j)
        end do
        write (*,*)
    end do
contains
    function objective(x) result(f)
        real(dp), intent(in) :: x(:) !! bounded search coordinates
        real(dp) :: f !! objective selected by case identifier
        real(dp) :: first, second
        select case (case_id)
        case (1)
            f = (x(1)**2-0.25_dp)**2
        case (2)
            f = (x(1)-0.3_dp)**2
        case (3)
            f = sum(x**2)
        case (4)
            first = (x(1)-0.5_dp)**2+(x(2)+0.5_dp)**2
            second = (x(1)+0.5_dp)**2+(x(2)-0.5_dp)**2
            f = min(first,second)
        end select
    end function objective
end program optimize_shgo_hgr_driver
