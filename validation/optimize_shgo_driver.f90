! SPDX-License-Identifier: MIT
! Stand-alone SciFort SHGO/SciPy benchmark exchange driver.
program optimize_shgo_driver
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_optimize_shgo, only : shgo_core, shgo_result, shgo_objective
    implicit none
    integer :: cid, levels, dim, j
    real(dp) :: bounds(3,2)
    type(shgo_result) :: res

    do cid = 1, 11
        dim = 2
        if (cid == 1 .or. cid == 7 .or. cid == 10) dim = 1
        if (cid == 8) dim = 3
        bounds(:,1) = -2.0_dp
        bounds(:,2) = 2.0_dp
        select case (cid)
        case (4)
            bounds = reshape([-5.12_dp,-5.12_dp,-5.12_dp, &
                5.12_dp,5.12_dp,5.12_dp],[3,2])
        case (5)
            bounds(:,1) = -5.0_dp
            bounds(:,2) = 5.0_dp
        case (6)
            bounds(:,1) = -6.0_dp
            bounds(:,2) = 6.0_dp
        end select
        do levels = 1, 3
            res = shgo_core(objective,bounds(1:dim,:),iters=levels, &
                maxfun=10000,local_iter=20,minimizer=coordinate_polish)
            write (*,'(3(i0,1x),es24.16,1x,5(i0,1x))',advance='no') &
                cid,levels,res%status,res%fun,res%nfev,res%nlfev,res%nit,res%nlocal,res%nminima
            do j = 1, dim
                write (*,'(es24.16,1x)',advance='no') res%x(j)
            end do
            write (*,*)
        end do
    end do
contains
    function objective(x) result(f)
        real(dp), intent(in) :: x(:) !! parameter vector of length dim
        real(dp) :: f !! one of eleven independently comparable objectives
        real(dp) :: a, b, m
        f = huge(1.0_dp)
        select case (cid)
        case (1)
            f = (x(1)-0.317_dp)**2
        case (2)
            f = (x(1)-0.43_dp)**2 + 3.0_dp*(x(2)+0.61_dp)**2
        case (3)
            f = (x(1)+x(2)-0.3_dp)**2 + 0.1_dp*(x(1)-x(2)-1.2_dp)**2
        case (4)
            f = 20.0_dp + sum(x**2) - 10.0_dp*sum(cos(2.0_dp*acos(-1.0_dp)*x))
        case (5)
            m = real(size(x), dp)
            f = -20.0_dp*exp(-0.2_dp*sqrt(sum(x**2)/m)) &
                - exp(sum(cos(2.0_dp*acos(-1.0_dp)*x))/m) + 20.0_dp + exp(1.0_dp)
        case (6)
            f = (x(1)**2+x(2)-11.0_dp)**2 + (x(1)+x(2)**2-7.0_dp)**2
        case (7)
            f = (x(1)**2-1.0_dp)**2 + 0.1_dp*x(1)
        case (8)
            f = sum((x-[0.25_dp,-0.5_dp,0.75_dp])**2)
        case (9)
            f = sin(5.0_dp*x(1))**2 + (x(1)-0.71_dp)**2 &
                + sin(3.0_dp*x(2))**2+(x(2)+0.27_dp)**2
        case (10)
            f = abs(x(1)-0.123_dp)
        case (11)
            a = (x(1)-1.0_dp)**2 + (x(2)-1.0_dp)**2
            b = (x(1)+1.0_dp)**2 + (x(2)+1.0_dp)**2 + 0.25_dp
            f = min(a,b)
        end select
    end function objective

    subroutine coordinate_polish(fun,lower,upper,start,maxfev,x,value,nfev,success)
        procedure(shgo_objective) :: fun !! budget-guarded objective function
        real(dp), intent(in) :: lower(:) !! inclusive lower bounds
        real(dp), intent(in) :: upper(:) !! inclusive upper bounds
        real(dp), intent(in) :: start(:) !! initial feasible seed
        integer, intent(in) :: maxfev !! objective evaluation budget
        real(dp), intent(out) :: x(:) !! locally polished point
        real(dp), intent(out) :: value !! locally polished objective
        integer, intent(out) :: nfev !! number of objective calls
        logical, intent(out) :: success !! whether a finite point is available
        real(dp) :: step(size(x)), trial(size(x)), ftrial
        integer :: iter, k, dir
        logical :: improved
        x = start
        nfev = 0
        success = .false.
        value = huge(1.0_dp)
        if (maxfev < 1) return
        value = fun(x)
        nfev = 1
        if (.not. ieee_is_finite(value)) return
        step = 0.15_dp*(upper-lower)
        do iter = 1, 250
            improved = .false.
            do k = 1, size(x)
                do dir = -1, 1, 2
                    trial = x
                    trial(k) = max(lower(k),min(upper(k),x(k)+real(dir,dp)*step(k)))
                    if (trial(k) == x(k)) cycle
                    if (nfev >= maxfev) exit
                    ftrial = fun(trial)
                    nfev = nfev + 1
                    if (ieee_is_finite(ftrial) .and. ftrial < value) then
                        value = ftrial
                        x = trial
                        improved = .true.
                    end if
                end do
                if (nfev >= maxfev) exit
            end do
            if (nfev >= maxfev) exit
            if (.not. improved) step = step * 0.5_dp
            if (maxval(step/(upper-lower)) < 1.0e-9_dp) exit
        end do
        success = .true.
    end subroutine coordinate_polish
end program optimize_shgo_driver
