! SPDX-License-Identifier: MIT
program linalg_qr_multiply_example
    use scifort_kinds, only : dp
    use scifort_linalg_qr, only : qr_factor
    use scifort_linalg_qr_multiply, only : qr_multiply, qr_apply
    implicit none
    real(dp) :: a(4,3), c(3,2), product(4,2), r(3,3)
    real(dp) :: af(4,3), tau(3), reused(4,2)
    complex(dp) :: ac(4,3), cc(2,4), complex_product(2,3), rc(3,3)
    integer :: status, pivots(3)

    a(:,1) = [1.0_dp,2.0_dp,3.0_dp,4.0_dp]
    a(:,2) = [2.0_dp,0.0_dp,-1.0_dp,3.0_dp]
    a(:,3) = [0.0_dp,2.0_dp,3.0_dp,1.0_dp]
    c(:,1) = [1.0_dp,2.0_dp,0.0_dp]
    c(:,2) = [0.0_dp,1.0_dp,1.0_dp]

    call qr_multiply(a,c,product,r,status,mode='left')
    if (status /= 0) error stop 'qr_multiply failed'
    print '(a)', 'Q @ C without forming Q:'
    print '(3f14.7)', product

    call qr_factor(a,af,tau,status)
    if (status /= 0) error stop 'qr_factor failed'
    call qr_apply(af,tau,c,reused,status,mode='left')
    if (status /= 0) error stop 'qr_apply failed'
    print '(a,es12.4)', 'Reused-factor discrepancy: ',maxval(abs(reused-product))

    ac = cmplx(a,0.25_dp*a,dp)
    cc = cmplx(1.0_dp,2.0_dp,dp)
    call qr_multiply(ac,cc,complex_product,rc,status,mode='right', &
        pivots=pivots,conjugate=.true.)
    if (status /= 0) error stop 'complex qr_multiply failed'
    print '(a,3(i0,1x))', 'Pivot order (one-based): ',pivots
    print '(a,es12.4)', 'Norm of complex C @ conjugate(Q): ', &
        sqrt(sum(abs(complex_product)**2))
end program linalg_qr_multiply_example
