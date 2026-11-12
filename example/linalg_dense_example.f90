! SPDX-License-Identifier: MIT
program linalg_dense_example
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success
    use scifort_linalg_dense, only : solve_triangular, lu_factor, lu_solve
    implicit none
    real(dp) :: a(3,3), lu(3,3), b(3,2), x(3,2), exact(3,2)
    complex(dp) :: t(3,3), bc(3), xc(3), expected(3)
    integer :: piv(3), status

    ! Factor once and reuse for multiple right-hand sides.
    a(1,:)=[0.0_dp,1.0_dp,1.0_dp]
    a(2,:)=[5.0_dp,4.0_dp,1.0_dp]
    a(3,:)=[1.0_dp,1.0_dp,3.0_dp]
    exact(:,1)=[1.0_dp,2.0_dp,-1.0_dp]
    exact(:,2)=[2.0_dp,-1.0_dp,0.5_dp]
    b=matmul(a,exact)
    call lu_factor(a,lu,piv,status)
    if (status /= linalg_status_success) error stop 'factorization failed'
    call lu_solve(lu,piv,b,x,status)
    if (status /= linalg_status_success) error stop 'solve failed'
    print '(A,ES12.4)', 'dense LU residual: ',maxval(abs(x-exact))

    ! Reuse factors to solve A^T x = b.
    b=matmul(transpose(a),exact)
    call lu_solve(lu,piv,b,x,status,trans='T')
    if (status /= linalg_status_success) error stop 'transpose solve failed'
    print '(A,ES12.4)', 'transpose LU residual: ',maxval(abs(x-exact))

    ! Upper-triangular, complex conjugate transpose with unit diagonal.
    t=cmplx(0.0_dp,0.0_dp,dp)
    t(1,2)=cmplx(0.5_dp,0.2_dp,dp)
    t(1,3)=cmplx(-0.3_dp,0.1_dp,dp)
    t(2,3)=cmplx(0.7_dp,-0.1_dp,dp)
    expected=[cmplx(1.0_dp,0.1_dp,dp),cmplx(-1.0_dp,0.5_dp,dp), &
        cmplx(0.4_dp,-0.2_dp,dp)]
    t(1,1)=cmplx(1.0_dp,0.0_dp,dp)
    t(2,2)=cmplx(1.0_dp,0.0_dp,dp)
    t(3,3)=cmplx(1.0_dp,0.0_dp,dp)
    bc=matmul(transpose(conjg(t)),expected)
    t(1,1)=cmplx(0.0_dp,0.0_dp,dp)
    t(2,2)=cmplx(0.0_dp,0.0_dp,dp)
    t(3,3)=cmplx(0.0_dp,0.0_dp,dp)
    call solve_triangular(t,bc,xc,status,trans='C',unit_diagonal=.true.)
    if (status /= linalg_status_success) error stop 'triangular solve failed'
    print '(A,ES12.4)', 'triangular adjoint residual: ',maxval(abs(xc-expected))
end program linalg_dense_example
