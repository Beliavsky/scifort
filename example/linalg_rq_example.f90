! SPDX-License-Identifier: MIT
program linalg_rq_example
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success
    use scifort_linalg_rq, only : rq, rq_r
    implicit none
    real(dp) :: a(2,3), r(2,3), q(3,3), re(2,2), qe(2,3), rr(2,3)
    complex(dp) :: ac(2,3), rc(2,2), qc(2,3)
    integer :: status

    a(1,:) = [1.0_dp, 3.0_dp, 5.0_dp]
    a(2,:) = [2.0_dp, 4.0_dp, 6.0_dp]

    call rq(a,r,q,status)
    if (status /= linalg_status_success) error stop 'full RQ failed'
    print '(a,es12.3)', 'Full real RQ reconstruction: ',maxval(abs(a-matmul(r,q)))
    call rq(a,re,qe,status,mode='economic')
    if (status /= linalg_status_success) error stop 'economic RQ failed'
    print '(a,es12.3)', 'Economic real RQ reconstruction: ',maxval(abs(a-matmul(re,qe)))
    call rq_r(a,rr,status)
    if (status /= linalg_status_success) error stop 'R-only failed'
    print '(a,es12.3)', 'R-only difference: ',maxval(abs(r-rr))

    ac = cmplx(a,0.5_dp*a,dp)
    call rq(ac,rc,qc,status,mode='economic')
    if (status /= linalg_status_success) error stop 'complex RQ failed'
    print '(a,es12.3)', 'Economic complex RQ reconstruction: ',maxval(abs(ac-matmul(rc,qc)))
end program linalg_rq_example
