! SPDX-License-Identifier: MIT
program linalg_toeplitz_example
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success
    use scifort_linalg_toeplitz, only : toeplitz, solve_toeplitz
    implicit none
    real(dp) :: c(4), r(4), b(4), x(4), t(4,4)
    complex(dp) :: cc(3), zz(3,3), xx(3,2), rhs(3,2), exact(3,2)
    integer :: i, status

    ! Example from SciPy's documented solve_toeplitz.
    c = [1.0_dp,3.0_dp,6.0_dp,10.0_dp]
    r = [1.0_dp,-1.0_dp,-2.0_dp,-3.0_dp]
    b = [1.0_dp,2.0_dp,2.0_dp,5.0_dp]
    call solve_toeplitz(c,b,x,status,r)
    if (status /= linalg_status_success) error stop 'real Toeplitz solve failed'
    print '(a,*(f12.7,1x))', 'real x: ',x
    call toeplitz(c,t,status,r)
    if (status /= linalg_status_success) error stop 'real Toeplitz constructor failed'
    print '(a,es12.4)', 'real residual: ',maxval(abs(matmul(t,x)-b))

    ! No r: the default row is conjg(c), with a preserved complex diagonal.
    cc = [cmplx(4.0_dp,0.5_dp,dp),cmplx(0.2_dp,0.6_dp,dp), &
        cmplx(-0.3_dp,0.1_dp,dp)]
    call toeplitz(cc,zz,status)
    if (status /= linalg_status_success) error stop 'complex constructor failed'
    do i = 1, 3
        exact(i,1) = cmplx(real(i,dp),-0.2_dp*real(i,dp),dp)
        exact(i,2) = cmplx(0.3_dp*real(i,dp),real(i,dp),dp)
    end do
    rhs = matmul(zz,exact)
    call solve_toeplitz(cc,rhs,xx,status)
    if (status /= linalg_status_success) error stop 'complex Toeplitz solve failed'
    print '(a,es12.4)', 'complex matrix RHS residual: ', &
        maxval(abs(matmul(zz,xx)-rhs))
end program linalg_toeplitz_example
