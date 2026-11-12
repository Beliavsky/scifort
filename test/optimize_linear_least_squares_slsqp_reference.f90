! SPDX-License-Identifier: MIT
! These tests require the actual SLSQP 1.6.0 FPM dependency.
program optimize_linear_least_squares_slsqp_reference
    use scifort_kinds, only : dp
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf
    use scifort_optimize_linear_least_squares_slsqp
    use scifort_optimize_linear_least_squares, only: linear_lsq_result
    implicit none
    real(dp) :: a(3,2),b(3),x(2),rnorm,inf
    type(linear_lsq_result) :: result
    integer :: status
    a(:,1)=[1.0_dp,1.0_dp,0.0_dp]
    a(:,2)=[0.0_dp,0.0_dp,1.0_dp]
    b=[2.0_dp,1.0_dp,1.0_dp]
    inf=ieee_value(1.0_dp,ieee_positive_inf)
    call nnls_slsqp(a,b,x,rnorm,status)
    if (status /= 0 .or. maxval(abs(x-[1.5_dp,1.0_dp])) > 1.0e-10_dp) &
        error stop 'SLSQP nnls mismatch'
    if (abs(rnorm-sqrt(0.5_dp)) > 1.0e-10_dp) error stop 'SLSQP nnls residual'
    result=lsq_linear_slsqp(a,b,[0.0_dp,0.0_dp],[1.0_dp,inf])
    if (.not.result%success) error stop 'SLSQP BVLS failure'
    if (maxval(abs(result%x-[1.0_dp,1.0_dp]))>1.0e-10_dp) error stop 'SLSQP BVLS vector'
    if (result%optimality>1.0e-10_dp) error stop 'SLSQP BVLS KKT'
    print *, 'SLSQP BVLS bridge: passed'
end program optimize_linear_least_squares_slsqp_reference
