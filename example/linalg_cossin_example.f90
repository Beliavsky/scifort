! SPDX-License-Identifier: MIT
! Complete real CS decomposition of a 4x4 rotation across two partitions.
program linalg_cossin_example
    use scifort_kinds, only : dp
    use scifort_linalg_cossin, only : cossin, cossin_separate
    use scifort_linalg, only : linalg_status_success
    implicit none
    real(dp) :: x(4,4), u(4,4), cs(4,4), vt(4,4), theta(2)
    real(dp) :: u1(2,2),u2(2,2),v1t(2,2),v2t(2,2), t1, t2, error
    integer :: status

    t1 = 0.3_dp
    t2 = 0.7_dp
    x = 0.0_dp
    x(1,1) = cos(t1)
    x(2,2) = cos(t2)
    x(3,3) = cos(t1)
    x(4,4) = cos(t2)
    x(1,3) = -sin(t1)
    x(2,4) = -sin(t2)
    x(3,1) = sin(t1)
    x(4,2) = sin(t2)

    call cossin(x,2,2,u,cs,vt,status)
    if (status /= linalg_status_success) error stop 'cossin failed'
    error = maxval(abs(matmul(matmul(u,cs),vt)-x))
    print '(a,es12.4)', 'max |U*CS*V^T-X| = ', error

    call cossin_separate(x,2,2,u1,u2,theta,v1t,v2t,status)
    if (status /= linalg_status_success) error stop 'cossin_separate failed'
    print '(a,2(f8.5,1x))', 'CS angles (radians): ',theta
end program linalg_cossin_example
