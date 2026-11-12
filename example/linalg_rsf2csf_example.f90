! SPDX-License-Identifier: MIT
! Compute the real Schur form, then convert it without another Schur QR run.
program linalg_rsf2csf_example
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success
    use scifort_linalg_schur, only : schur, rsf2csf
    implicit none

    real(dp) :: a(4,4), tr(4,4), zr(4,4), relerr
    complex(dp) :: tc(4,4), zc(4,4)
    integer :: i, status

    a = 0.0_dp
    a(1,1) = -2.0_dp
    a(2,2) = 0.5_dp
    a(2,3) = -3.0_dp
    a(3,2) = 3.0_dp
    a(3,3) = 0.5_dp
    a(4,4) = 4.0_dp
    a(1,2) = 0.6_dp
    a(1,3) = -0.2_dp
    a(2,4) = 0.7_dp
    a(3,4) = -0.4_dp

    call schur(a, tr, zr, status)
    if (status /= linalg_status_success) error stop 'real Schur failed'
    call rsf2csf(tr, zr, tc, zc, status)
    if (status /= linalg_status_success) error stop 'rsf2csf failed'

    print '(a)', 'Eigenvalues (diagonal of complex Schur factor):'
    do i = 1, size(a,1)
        print '(i2,2es17.8)', i, real(tc(i,i), dp), aimag(tc(i,i))
    end do
    relerr = maxval(abs(matmul(matmul(zc,tc),transpose(conjg(zc)))-cmplx(a,0.0_dp,dp))) &
        / max(1.0_dp,maxval(abs(a)))
    print '(a,es13.5)', 'Maximum normalized reconstruction error: ', relerr
    if (relerr > 1.0e-10_dp) error stop 'reconstruction failed'
end program linalg_rsf2csf_example
