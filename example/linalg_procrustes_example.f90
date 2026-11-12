! SPDX-License-Identifier: MIT
program linalg_procrustes_example
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success
    use scifort_linalg_procrustes, only : orthogonal_procrustes
    implicit none
    real(dp) :: a(4,3), b(4,3), q(3,3), r(3,3), scale_factor
    complex(dp) :: bc(4,3), qc(3,3), rc(3,3)
    integer :: i, status

    a = 0.0_dp
    do i = 1, 3
        a(i,i) = real(i,dp)
    end do
    q = 0.0_dp
    q(1,2) = 1.0_dp
    q(2,3) = 1.0_dp
    q(3,1) = 1.0_dp
    b = matmul(a,q)
    call orthogonal_procrustes(a,b,r,scale_factor,status)
    if (status /= linalg_status_success) error stop 'real Procrustes failed'
    print '(A,ES12.4)', 'real alignment error: ', maxval(abs(matmul(a,r)-b))
    print '(A,F6.2)', 'real singular-value sum: ', scale_factor

    qc = cmplx(q,0.0_dp,dp)
    qc(1,2) = cmplx(0.0_dp,1.0_dp,dp)
    bc = matmul(cmplx(a,0.0_dp,dp),qc)
    call orthogonal_procrustes(a,bc,rc,scale_factor,status)
    if (status /= linalg_status_success) error stop 'complex Procrustes failed'
    print '(A,ES12.4)', 'mixed alignment error: ', &
        maxval(abs(matmul(cmplx(a,0.0_dp,dp),rc)-bc))
    print '(A,F6.2)', 'complex singular-value sum: ', scale_factor
end program linalg_procrustes_example
