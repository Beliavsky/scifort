! SPDX-License-Identifier: MIT
! Output full solutions for independent SciPy validation.
program optimize_direct_driver
    use scifort_kinds, only : dp
    use scifort_optimize_direct, only : direct, direct_result
    implicit none
    integer :: cid, biased, n, budget, ibudget
    real(dp), allocatable :: bounds(:,:)
    type(direct_result) :: r

    do cid = 1, 10
        call setup(cid, bounds)
        n = size(bounds, 1)
        do ibudget = 1, 3
            select case(ibudget)
            case(1)
                budget = 120
            case(2)
                budget = 480
            case default
                budget = 1600
            end select
            do biased = 0, 1
                r = direct(objective, bounds, maxfun=budget, maxiter=180, &
                    locally_biased=(biased==1), vol_tol=0.0_dp, len_tol=0.0_dp)
                if (.not. allocated(r%x)) error stop 'direct failed unexpectedly'
                write (*,'(4(i0,1x),es24.16,1x,2(i0,1x),*(es24.16,1x))') &
                    cid,biased,budget,r%status,r%fun,r%nfev,r%nit,r%x
            end do
        end do
    end do
contains
    subroutine setup(id, range)
        integer, intent(in) :: id !! benchmark identifier
        real(dp), allocatable, intent(out) :: range(:,:) !! bounds for each variable
        integer :: num
        select case(id)
        case(1)
            num=1
        case(8)
            num=4
        case(9)
            num=3
        case default
            num=2
        end select
        allocate(range(num,2))
        select case(id)
        case(1)
            range(:,1)=-1.0_dp; range(:,2)=1.0_dp
        case(4,5)
            range(:,1)=-5.0_dp; range(:,2)=5.0_dp
        case(6)
            range(:,1)=-2.0_dp; range(:,2)=2.0_dp
        case(7)
            range(:,1)=-6.0_dp; range(:,2)=6.0_dp
        case default
            range(:,1)=-2.0_dp; range(:,2)=2.0_dp
        end select
    end subroutine setup
    function objective(x) result(y)
        real(dp), intent(in) :: x(:) !! test point
        real(dp) :: y !! scalar objective
        real(dp) :: pi
        pi = acos(-1.0_dp)
        y = huge(1.0_dp)
        select case (cid)
        case(1)
            y = (x(1)-0.317_dp)**2
        case(2)
            y = (x(1)-0.4_dp)**2 + 3.0_dp*(x(2)+0.7_dp)**2
        case(3)
            y = (x(1)+x(2)-0.3_dp)**2 + 0.1_dp*(x(1)-x(2)-1.2_dp)**2
        case(4)
            y = 20.0_dp + x(1)**2 + x(2)**2 &
                - 10.0_dp*(cos(2.0_dp*pi*x(1))+cos(2.0_dp*pi*x(2)))
        case(5)
            y = -20.0_dp*exp(-0.2_dp*sqrt(sum(x*x)/2.0_dp)) &
                - exp(sum(cos(2.0_dp*pi*x))/2.0_dp) + 20.0_dp + exp(1.0_dp)
        case(6)
            y = 100.0_dp*(x(2)-x(1)**2)**2 + (1.0_dp-x(1))**2
        case(7)
            y = (x(1)**2+x(2)-11.0_dp)**2 + (x(1)+x(2)**2-7.0_dp)**2
        case(8)
            y = sum((x - [0.25_dp,-0.5_dp,0.75_dp,-0.25_dp])**2)
        case(9)
            y = sum(abs(x-[0.18_dp,-0.71_dp,0.34_dp]))
        case(10)
            y = sin(5.0_dp*x(1))**2 + (x(1)-0.71_dp)**2 &
                + sin(3.0_dp*x(2))**2 + (x(2)+0.27_dp)**2
        end select
    end function objective
end program optimize_direct_driver
