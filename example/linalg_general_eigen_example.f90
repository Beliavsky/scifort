! SPDX-License-Identifier: MIT
program linalg_general_eigen_example
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success
    use scifort_linalg_general_eigen, only : eig, eigvals
    implicit none

    real(dp) :: rotation(2,2)
    complex(dp) :: w(2), vr(2,2), vl(2,2), a(2,2)
    integer :: status, i

    rotation(1,:) = [0.0_dp,-1.0_dp]
    rotation(2,:) = [1.0_dp,0.0_dp]
    call eig(rotation,w,status,vr=vr,vl=vl)
    if (status /= linalg_status_success) error stop 'real eigenproblem failed'
    print '(a)', 'Eigenvalues of 90-degree planar rotation:'
    do i = 1, 2
        print '(2(es15.6,1x))', real(w(i),dp),aimag(w(i))
    end do

    a(1,:) = [cmplx(1.0_dp,2.0_dp,dp), cmplx(3.0_dp,1.0_dp,dp)]
    a(2,:) = [cmplx(0.0_dp,0.0_dp,dp), cmplx(4.0_dp,-1.0_dp,dp)]
    call eigvals(a,w,status)
    if (status /= linalg_status_success) error stop 'complex eigenvalues failed'
    print '(a)', 'Eigenvalues of complex triangular matrix:'
    do i = 1, 2
        print '(2(es15.6,1x))', real(w(i),dp),aimag(w(i))
    end do
end program linalg_general_eigen_example
