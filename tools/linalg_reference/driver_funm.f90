! SPDX-License-Identifier: MIT
! Batch input/output driver for independent scipy.linalg.funm comparison.
program driver_funm
    use scifort_kinds, only : dp
    use scifort_linalg_funm, only : funm
    implicit none
    real(dp), allocatable :: real_input(:, :)
    complex(dp), allocatable :: a(:, :), x(:, :)
    real(dp) :: re, im, commutator
    integer :: case_count, kind_input, n, i, j, k, status, function_id, with_derivative

    read(*,*) case_count
    do k = 1, case_count
        read(*,*) kind_input, n, function_id, with_derivative
        allocate(a(n,n), x(n,n), real_input(n,n))
        do j = 1, n
            do i = 1, n
                read(*,*) re, im
                a(i,j) = cmplx(re, im, dp)
            end do
        end do
        if (kind_input == 1) then
            real_input = real(a,dp)
            if (with_derivative == 2) then
                call funm(real_input, scalar_function, x, status, &
                    commutator_residual=commutator, separation_tol=0.02_dp, &
                    derivative=scalar_derivative, taylor_max_order=96)
            else if (with_derivative == 1) then
                call funm(real_input, scalar_function, x, status, &
                    commutator_residual=commutator, derivative=scalar_derivative)
            else
                call funm(real_input, scalar_function, x, status, commutator_residual=commutator)
            end if
        else
            if (with_derivative == 2) then
                call funm(a, scalar_function, x, status, &
                    commutator_residual=commutator, separation_tol=0.02_dp, &
                    derivative=scalar_derivative, taylor_max_order=96)
            else if (with_derivative == 1) then
                call funm(a, scalar_function, x, status, &
                    commutator_residual=commutator, derivative=scalar_derivative)
            else
                call funm(a, scalar_function, x, status, commutator_residual=commutator)
            end if
        end if
        write(*,'(i0,1x,es29.17e3)') status, commutator
        do j = 1, n
            do i = 1, n
                write(*,'(es29.17e3,1x,es29.17e3)') real(x(i,j),dp), aimag(x(i,j))
            end do
        end do
        deallocate(a,x,real_input)
    end do

contains
    function scalar_function(z) result(value)
        complex(dp), intent(in) :: z !! complex scalar function argument
        complex(dp) :: value
        select case (function_id)
        case (1)
            value = exp(z)
        case (2)
            value = log(z)
        case (3)
            value = sin(z)
        case (4)
            value = z*z
        case (5)
            value = 1.0_dp/(1.0_dp+z)
        case (6)
            value = cos(z)
        case (7)
            value = exp(-z*z)
        case (8)
            value = sqrt(z)
        case (9)
            value = sinh(z)
        case default
            value = z
        end select
    end function scalar_function

    function scalar_derivative(z, order) result(value)
        complex(dp), intent(in) :: z !! derivative evaluation point
        integer, intent(in) :: order !! positive order of derivative
        complex(dp) :: value
        integer :: i
        select case(function_id)
        case (1)
            value = exp(z)
        case (3)
            select case (mod(order,4))
            case (0)
                value = sin(z)
            case (1)
                value = cos(z)
            case (2)
                value = -sin(z)
            case default
                value = -cos(z)
            end select
        case (6)
            select case (mod(order,4))
            case (0)
                value = cos(z)
            case (1)
                value = -sin(z)
            case (2)
                value = -cos(z)
            case default
                value = sin(z)
            end select
        case (9)
            if (mod(order,2) == 0) then
                value = sinh(z)
            else
                value = cosh(z)
            end if
        case (4)
            value = cmplx(0.0_dp,0.0_dp,dp)
            if (order == 1) value = 2.0_dp*z
            if (order == 2) value = cmplx(2.0_dp,0.0_dp,dp)
        case (2)
            value = 1.0_dp/z**order
            do i = 2, order
                value = -real(i-1,dp)*value
            end do
        case default
            value = cmplx(0.0_dp,0.0_dp,dp)
        end select
    end function scalar_derivative
end program driver_funm
