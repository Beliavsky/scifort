! SPDX-License-Identifier: MIT
program linalg_dense_cholesky_example
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success
    use scifort_linalg_dense, only : cholesky, cho_factor, cho_solve
    implicit none
    real(dp) :: a(3,3), f(3,3), b(3,2), x(3,2), truth(3,2)
    complex(dp) :: ac(3,3), fc(3,3), bc(3), xc(3), expected(3)
    integer :: status, i

    a(1,:)=[5.0_dp,1.0_dp,0.0_dp]
    a(2,:)=[1.0_dp,4.0_dp,1.0_dp]
    a(3,:)=[0.0_dp,1.0_dp,3.0_dp]
    truth(:,1)=[1.0_dp,-1.0_dp,2.0_dp]
    truth(:,2)=[0.5_dp,2.0_dp,-3.0_dp]
    b=matmul(a,truth)
    call cho_factor(a,f,status,lower=.true.)
    if (status /= linalg_status_success) error stop 'real factorization failed'
    call cho_solve(f,b,x,status,lower=.true.)
    if (status /= linalg_status_success) error stop 'real solve failed'
    print '(A,ES12.4)', 'real Cholesky solve error: ',maxval(abs(x-truth))

    ac=cmplx(a,0.0_dp,dp)
    ac(1,2)=cmplx(1.0_dp,0.5_dp,dp)
    ac(2,1)=conjg(ac(1,2))
    expected=[cmplx(1.0_dp,0.5_dp,dp),cmplx(-0.5_dp,1.0_dp,dp), &
        cmplx(0.25_dp,-1.0_dp,dp)]
    bc=matmul(ac,expected)
    call cholesky(ac,fc,status) ! upper-triangular factor U, default lower=.false.
    if (status /= linalg_status_success) error stop 'complex factorization failed'
    call cho_solve(fc,bc,xc,status)
    if (status /= linalg_status_success) error stop 'complex solve failed'
    print '(A,ES12.4)', 'complex Cholesky solve error: ',maxval(abs(xc-expected))

    ! One factorization can be reused for any number of right-hand sides.
    do i=1,3
        bc(i)=cmplx(real(i,dp),-0.5_dp,dp)
    end do
    call cho_solve(fc,bc,xc,status)
    if (status /= linalg_status_success) error stop 'complex factor reuse failed'
    print '(A,ES12.4)', 'complex reuse residual: ',maxval(abs(matmul(ac,xc)-bc))
end program linalg_dense_cholesky_example
