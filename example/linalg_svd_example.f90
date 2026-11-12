! SPDX-License-Identifier: MIT
program linalg_svd_example
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success
    use scifort_linalg_svd, only : svd, pinv
    use scifort_linalg_qr, only : lstsq
    implicit none
    real(dp) :: a(3,2), b(3), u(3,2), s(2), vh(2,2), x(2), ap(2,3)
    integer :: status, rank

    a = 0.0_dp
    a(1,:) = [1.0_dp,1.0_dp]
    b = [2.0_dp,1.0_dp,0.0_dp]
    call svd(a,u,s,vh,status)
    if (status /= linalg_status_success) error stop 'SVD did not converge'
    print '(A,2(1X,ES12.4))', 'singular values:',s
    call pinv(a,ap,status,rank)
    if (status /= linalg_status_success) error stop 'pinv failed'
    print '(A,I0)', 'numerical rank: ',rank
    call lstsq(a,b,x,status,rank)
    if (status /= linalg_status_success) error stop 'minimum-norm lstsq failed'
    print '(A,2(1X,F12.6))', 'minimum-norm solution:',x
    print '(A,1X,ES12.4)', 'least-squares residual:',sum(abs(matmul(a,x)-b))
end program linalg_svd_example
