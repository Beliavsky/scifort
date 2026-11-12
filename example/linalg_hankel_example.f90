! SPDX-License-Identifier: MIT
! Hankel construction and FFT-based Hankel multiplication.
program linalg_hankel_example
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success
    use scifort_linalg_hankel, only : hankel
    use scifort_linalg_matmul_hankel, only : matmul_hankel
    implicit none
    real(dp) :: c(4), r(5), a(4,5)
    real(dp) :: x(5), y(4)
    integer :: status, i

    c = [1.0_dp,2.0_dp,3.0_dp,4.0_dp]
    r = [-99.0_dp,7.0_dp,7.0_dp,8.0_dp,9.0_dp]
    x = [2.0_dp,3.0_dp,4.0_dp,5.0_dp,6.0_dp]

    call hankel(c,a,status,r)
    if (status /= linalg_status_success) error stop 'hankel construction failed'
    print '(a)', 'Hankel matrix:'
    do i=1,size(a,1)
        print '(*(f8.2,1x))', a(i,:)
    end do

    call matmul_hankel(c,x,y,status,r)
    if (status /= linalg_status_success) error stop 'FFT Hankel multiplication failed'
    print '(a)', 'H*x (FFT algorithm):'
    print '(*(f10.3,1x))', y
    if (maxval(abs(y-matmul(a,x))) > 1.0e-11_dp) error stop 'Hankel product differs from dense'
end program linalg_hankel_example
