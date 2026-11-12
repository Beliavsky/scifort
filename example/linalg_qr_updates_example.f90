! SPDX-License-Identifier: MIT
program linalg_qr_updates_example
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success
    use scifort_linalg_qr, only : qr
    use scifort_linalg_qr_updates, only : qr_update, qr_insert, qr_delete
    implicit none
    real(dp) :: a(5,3), q(5,3), r(3,3), q1(5,3), r1(3,3)
    real(dp) :: q2(5,4), r2(4,4), q3(5,3), r3(3,3)
    real(dp) :: u(5), v(3), new_column(5), expected(5,3), temporary(5,4)
    complex(dp) :: ac(3,2), qc(3,2), rc(2,2), qc1(3,2), rc1(2,2)
    complex(dp) :: uc(3), vc(2), expectedc(3,2)
    integer :: i, j, status

    do j = 1, 3
        do i = 1, 5
            a(i,j) = sin(real(i*5+j*2,dp))
        end do
    end do
    do i = 1, 5
        u(i) = cos(real(i,dp))
        new_column(i) = sin(real(i*3,dp))
    end do
    v = [0.15_dp,-0.5_dp,0.4_dp]
    call qr(a,q,r,status)
    if (status /= linalg_status_success) error stop 'initial QR'
    call qr_update(q,r,u,v,q1,r1,status)
    if (status /= linalg_status_success) error stop 'real QR rank-one update'
    expected = a+matmul(reshape(u,[5,1]),reshape(v,[1,3]))
    print '(A,ES12.4)', 'rank-one QR update residual: ',maxval(abs(matmul(q1,r1)-expected))

    ! Insert a column before the current second column, then remove it.
    call qr_insert(q1,r1,new_column,2,q2,r2,status,which='col')
    if (status /= linalg_status_success) error stop 'QR column insert'
    temporary(:,1) = expected(:,1)
    temporary(:,2) = new_column
    temporary(:,3:4) = expected(:,2:3)
    print '(A,ES12.4)', 'insert column residual:    ',maxval(abs(matmul(q2,r2)-temporary))
    call qr_delete(q2,r2,2,q3,r3,status,which='col')
    if (status /= linalg_status_success) error stop 'QR column delete'
    print '(A,ES12.4)', 'delete column residual:    ',maxval(abs(matmul(q3,r3)-expected))

    do j = 1, 2
        do i = 1, 3
            ac(i,j) = cmplx(sin(real(i*2+j,dp)),cos(real(i+j*3,dp)),dp)
        end do
    end do
    uc = [cmplx(1.0_dp,0.2_dp,dp),cmplx(-0.2_dp,1.0_dp,dp),cmplx(0.4_dp,0.1_dp,dp)]
    vc = [cmplx(0.1_dp,0.7_dp,dp),cmplx(-0.3_dp,0.1_dp,dp)]
    call qr(ac,qc,rc,status)
    if (status /= linalg_status_success) error stop 'initial complex QR'
    call qr_update(qc,rc,uc,vc,qc1,rc1,status)
    if (status /= linalg_status_success) error stop 'complex QR update'
    expectedc = ac+matmul(reshape(uc,[3,1]),reshape(conjg(vc),[1,2]))
    print '(A,ES12.4)', 'complex QR update residual: ',maxval(abs(matmul(qc1,rc1)-expectedc))
end program linalg_qr_updates_example
