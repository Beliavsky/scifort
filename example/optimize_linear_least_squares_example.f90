! SPDX-License-Identifier: MIT
program optimize_linear_least_squares_example
    use scifort_kinds, only : dp
    use scifort_optimize_linear_least_squares, only : nnls, lsq_linear, linear_lsq_result
    implicit none
    real(dp) :: a(3,2), b(3), x(2), residual_norm
    integer :: status
    type(linear_lsq_result) :: result
    a(:,1)=[1.0_dp,1.0_dp,0.0_dp]
    a(:,2)=[0.0_dp,0.0_dp,1.0_dp]
    b=[2.0_dp,1.0_dp,1.0_dp]
    call nnls(a,b,x,residual_norm,status)
    print '(a,2f12.6)', 'nnls solution: ',x
    print '(a,f12.6,1x,i0)', 'residual norm and status: ',residual_norm,status
    result=lsq_linear(a,b,lower=[0.0_dp,0.0_dp],upper=[1.0_dp,2.0_dp])
    print '(a,2f12.6)', 'bvls solution: ',result%x
    print '(a,es12.4,1x,i0)', 'cost and status: ',result%cost,result%status
end program optimize_linear_least_squares_example
