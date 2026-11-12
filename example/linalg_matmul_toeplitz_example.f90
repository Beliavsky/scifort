! SPDX-License-Identifier: MIT
program linalg_matmul_toeplitz_example
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success
    use scifort_linalg_matmul_toeplitz, only : matmul_toeplitz
    implicit none
    real(dp) :: c(3)
    real(dp) :: r(4)
    real(dp) :: x(4)
    real(dp) :: y(3)
    real(dp) :: rhs(4,2)
    real(dp) :: products(3,2)
    complex(dp) :: cc(3)
    complex(dp) :: cx(3)
    complex(dp) :: cy(3)
    integer :: status
    integer :: i

    c = [1.0_dp,2.0_dp,3.0_dp]
    r = [123.0_dp,4.0_dp,5.0_dp,6.0_dp] ! r(1) is ignored.
    x = [1.0_dp,2.0_dp,-1.0_dp,3.0_dp]
    call matmul_toeplitz(c,x,y,status,r)
    if (status /= linalg_status_success) error stop 'real Toeplitz multiplication failed'
    if (maxval(abs(y-[22.0_dp,15.0_dp,18.0_dp])) > 1.0e-12_dp) &
        error stop 'unexpected real Toeplitz product'
    print '(a,*(f12.5,1x))', 'Rectangular Toeplitz * vector: ',y

    rhs(:,1) = x
    rhs(:,2) = [0.0_dp,1.0_dp,2.0_dp,3.0_dp]
    call matmul_toeplitz(c,rhs,products,status,r)
    if (status /= linalg_status_success) error stop 'matrix multiplication failed'
    print '(a,*(f12.5,1x))', 'First output column: ',products(:,1)
    print '(a,*(f12.5,1x))', 'Second output column:',products(:,2)

    cc = [cmplx(2.0_dp,1.0_dp,dp), cmplx(0.2_dp,0.7_dp,dp), &
        cmplx(-0.5_dp,0.3_dp,dp)]
    cx = [cmplx(1.0_dp,2.0_dp,dp),cmplx(2.0_dp,0.0_dp,dp), &
        cmplx(-0.1_dp,1.0_dp,dp)]
    call matmul_toeplitz(cc,cx,cy,status)
    if (status /= linalg_status_success) error stop 'complex Toeplitz multiplication failed'
    print '(a)', 'Complex product with conjugate first row:'
    do i=1,size(cy)
        print '(i2,2(1x,f12.6))', i,real(cy(i),dp),aimag(cy(i))
    end do
end program linalg_matmul_toeplitz_example
