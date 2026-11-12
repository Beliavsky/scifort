! SPDX-License-Identifier: MIT
program linalg_cdf2rdf_example
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success
    use scifort_linalg_general_eigen, only : eig
    use scifort_linalg_cdf2rdf, only : cdf2rdf
    implicit none

    real(dp) :: a(3,3), wr(3,3), vr(3,3), residual
    complex(dp) :: w(3), v(3,3)
    integer :: status, i

    a = 0.0_dp
    a(1,1) = 1.0_dp
    a(1,2) = 2.0_dp
    a(1,3) = 3.0_dp
    a(2,2) = 4.0_dp
    a(2,3) = 5.0_dp
    a(3,2) = -5.0_dp
    a(3,3) = 4.0_dp

    call eig(a,w,status,vr=v)
    if (status /= linalg_status_success) error stop 'eig failed'
    call cdf2rdf(w,v,wr,vr,status)
    if (status /= linalg_status_success) error stop 'cdf2rdf failed'
    residual = maxval(abs(matmul(a,vr)-matmul(vr,wr)))

    print '(a)', 'Real block eigenvalue matrix:'
    do i = 1, 3
        print '(3(f12.6,1x))', wr(i,:)
    end do
    print '(a)', 'Real invariant eigenvector basis:'
    do i = 1, 3
        print '(3(f12.6,1x))', vr(i,:)
    end do
    print '(a,es13.5)', 'A*vr - vr*wr max error: ',residual
end program linalg_cdf2rdf_example
