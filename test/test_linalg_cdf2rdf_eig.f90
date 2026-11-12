! SPDX-License-Identifier: MIT
! Check interoperability with SciFort's self-contained nonsymmetric eigensolver.
program test_linalg_cdf2rdf_eig
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success
    use scifort_linalg_general_eigen, only : eig
    use scifort_linalg_cdf2rdf, only : cdf2rdf
    implicit none

    real(dp) :: a(3,3), wr(3,3), vr(3,3), b(2,2), br(2,2), bv(2,2)
    complex(dp) :: w(3), v(3,3), bw(2), bvec(2,2)
    integer :: status

    a = 0.0_dp
    a(1,1) = 2.0_dp
    a(2,2) = 3.0_dp
    a(3,3) = 3.0_dp
    a(2,3) = 5.0_dp
    a(3,2) = -5.0_dp
    call eig(a,w,status,vr=v)
    if (status /= linalg_status_success) error stop '3x3 eig failed'
    call cdf2rdf(w,v,wr,vr,status)
    if (status /= linalg_status_success) error stop '3x3 cdf2rdf after eig failed'
    if (maxval(abs(matmul(a,vr)-matmul(vr,wr))) > 1.0e-11_dp) &
        error stop '3x3 eig/cdf2rdf reconstruction failed'

    b(1,:) = [3.0_dp,5.0_dp]
    b(2,:) = [-5.0_dp,3.0_dp]
    call eig(b,bw,status,vr=bvec)
    if (status /= linalg_status_success) error stop '2x2 eig failed'
    call cdf2rdf(bw,bvec,br,bv,status)
    if (status /= linalg_status_success) error stop '2x2 cdf2rdf after eig failed'
    if (maxval(abs(matmul(b,bv)-matmul(bv,br))) > 1.0e-11_dp) &
        error stop '2x2 eig/cdf2rdf reconstruction failed'

    print '(a)', 'PASS cdf2rdf SciFort eig interoperability (2x2 and 3x3)'
end program test_linalg_cdf2rdf_eig
