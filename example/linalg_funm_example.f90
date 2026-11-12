! SPDX-License-Identifier: MIT
! Demonstrates block funm: Jordan blocks need derivatives; semisimple blocks do not.
program linalg_funm_example
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success
    use scifort_linalg_funm, only : funm
    implicit none
    real(dp) :: a(2,2), b(3,3)
    complex(dp) :: x(2,2), y(3,3)
    integer :: status

    a = 0.0_dp
    a(1,1) = 0.5_dp
    a(2,2) = 0.5_dp
    a(1,2) = 1.0_dp
    call funm(a, scalar_exp, x, status, derivative=exp_derivatives)
    if (status /= linalg_status_success) error stop 'funm failed'
    print '(a,2(1x,f11.7))', 'exp(A) first row:', real(x(1,1)), real(x(1,2))

    ! Repeated eigenvalues in positions 1 and 3 are not adjacent.
    b = 0.0_dp
    b(1,1) = 0.5_dp
    b(2,2) = 1.8_dp
    b(3,3) = 0.5_dp
    b(1,2) = 0.4_dp
    b(2,3) = 0.2_dp
    b(1,3) = -0.1_dp
    call funm(b, scalar_exp, y, status, derivative=exp_derivatives)
    if (status /= linalg_status_success) error stop 'nonadjacent funm failed'
    print '(a,3(1x,f11.7))', 'exp(B) first row:', real(y(1,1:3))

    ! Remove the induced Jordan coupling while retaining eigenvalue repeats
    ! at positions 1 and 3.  The result is diagonalizable: no derivative.
    b(1,3) = b(1,2)*b(2,3)/(b(2,2)-b(1,1))
    call funm(b, scalar_exp, y, status)
    if (status /= linalg_status_success) error stop 'semisimple funm failed'
    print '(a,3(1x,f11.7))', 'exp(C) first row:', real(y(1,1:3))

    ! Close but different eigenvalues need more than the Jordan polynomial.
    ! Here the callback supports derivatives of every positive order.
    a(1,1) = 0.5_dp
    a(2,2) = 0.503_dp
    a(1,2) = 0.8_dp
    call funm(a, scalar_exp, x, status, separation_tol=0.02_dp, &
        derivative=exp_derivatives, taylor_max_order=96)
    if (status /= linalg_status_success) error stop 'clustered funm failed'
    print '(a,2(1x,f11.7))', 'exp(D) first row:', real(x(1,1:2))

contains
    function scalar_exp(z) result(value)
        complex(dp), intent(in) :: z !! scalar argument of the exponential
        complex(dp) :: value
        value = exp(z)
    end function scalar_exp

    function exp_derivatives(z, order) result(value)
        complex(dp), intent(in) :: z !! evaluation point for k-th derivative
        integer, intent(in) :: order !! strictly positive order
        complex(dp) :: value
        value = exp(z)
    end function exp_derivatives
end program linalg_funm_example
