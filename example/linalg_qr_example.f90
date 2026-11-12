! SPDX-License-Identifier: MIT
program linalg_qr_example
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success
    use scifort_linalg_qr, only : qr_factor, qr_unpack, qr_solve, lstsq
    implicit none
    integer, parameter :: m=6, n=3
    real(dp) :: a(m,n), af(m,n), tau(n), q(m,n), r(n,n)
    real(dp) :: b(m,2), x(n,2), x2(n,2), truth(n,2)
    complex(dp) :: wide(3,5), rhs(3), minimum_norm(5), known(5)
    integer :: pivots(n), rank, status, i, j

    do j = 1, n
        do i = 1, m
            a(i,j) = sin(real(2*i*j+j,dp)) + 0.2_dp*cos(real(i*j,dp))
        end do
        truth(j,1) = real(j,dp)/3.0_dp
        truth(j,2) = -real(j,dp)/7.0_dp
    end do
    b = matmul(a,truth)
    call qr_factor(a,af,tau,status,pivots)
    if (status /= linalg_status_success) error stop 'QR factorization failed'
    call qr_unpack(af,tau,q,r,status)
    if (status /= linalg_status_success) error stop 'QR unpack failed'
    print '(A,ES12.4)', 'real QR reconstruction: ', maxval(abs(matmul(q,r)-a(:,pivots)))
    call qr_solve(af,tau,pivots,b,x,status,rank)
    if (status /= linalg_status_success) error stop 'QR solve failed'
    print '(A,I0,A,ES12.4)', 'real QR rank=',rank,', solution error=',maxval(abs(x-truth))
    call qr_solve(af,tau,pivots,2.0_dp*b,x2,status)
    if (status /= linalg_status_success) error stop 'second QR solve failed'
    print '(A,ES12.4)', 'factor reuse error: ',maxval(abs(x2-2.0_dp*truth))

    do j = 1, 5
        do i = 1, 3
            wide(i,j) = cmplx(sin(real(i*j,dp)),cos(real(i+2*j,dp)),dp)
        end do
    end do
    ! A^H*y is in the row space of A; it is the minimum-norm solution.
    known = matmul(transpose(conjg(wide)), &
        [cmplx(1.0_dp,0.5_dp,dp),cmplx(-0.2_dp,0.3_dp,dp),cmplx(0.7_dp,-0.1_dp,dp)])
    rhs = matmul(wide,known)
    call lstsq(wide,rhs,minimum_norm,status,rank)
    if (status /= linalg_status_success) error stop 'wide least squares failed'
    print '(A,I0,A,ES12.4)', 'complex wide rank=',rank,', min-norm error=', &
        maxval(abs(minimum_norm-known))
end program linalg_qr_example
