! SPDX-License-Identifier: MIT
program linalg_hessenberg_example
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success
    use scifort_linalg_hessenberg, only : hessenberg
    implicit none

    real(dp) :: a(4,4), h(4,4), q(4,4), residual
    complex(dp) :: ac(4,4), hc(4,4), qc(4,4)
    integer :: i, j, status

    do j = 1, 4
        do i = 1, 4
            a(i,j) = real(3*i-j*j,dp) / real(i+j,dp)
        end do
    end do
    call hessenberg(a,h,status,q)
    if (status /= linalg_status_success) error stop 'real Hessenberg reduction failed'
    residual = maxval(abs(a-matmul(matmul(q,h),transpose(q))))
    print '(A,ES12.4)', 'real Hessenberg residual: ',residual
    print '(A,ES12.4)', 'real structural zero maximum: ',maxval(abs(h(3:4,1)))

    ac = cmplx(a,0.05_dp*transpose(a),dp)
    call hessenberg(ac,hc,status,qc)
    if (status /= linalg_status_success) error stop 'complex Hessenberg reduction failed'
    residual = maxval(abs(ac-matmul(matmul(qc,hc),transpose(conjg(qc)))))
    print '(A,ES12.4)', 'complex Hessenberg residual: ',residual
    print '(A,ES12.4)', 'complex structural zero maximum: ',maxval(abs(hc(3:4,1)))
end program linalg_hessenberg_example
