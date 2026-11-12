! SPDX-License-Identifier: MIT
program linalg_interpolative_svd_example
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success
    use scifort_linalg_interpolative_svd, only : interpolative_svd
    implicit none
    real(dp) :: a(6,5), b(6,2), c(2,5), u(6,2), s(2), v(5,2)
    complex(dp) :: z(6,5)
    complex(dp), allocatable :: uz(:, :), vz(:, :)
    real(dp), allocatable :: sz(:)
    real(dp) :: err
    integer :: i, j, k, status

    do i = 1,6
        b(i,1) = sin(real(2*i,dp))
        b(i,2) = cos(real(3*i,dp))
    end do
    do j = 1,5
        c(1,j) = cos(real(5*j,dp))
        c(2,j) = sin(real(7*j,dp))
    end do
    a = matmul(b,c)
    call interpolative_svd(a,2,u,s,v,status)
    if (status /= linalg_status_success) error stop 'fixed-rank interpolative SVD failed'
    err = sqrt(sum((a-matmul(u,matmul(diagonal(s),transpose(v))))**2)) / &
        sqrt(sum(a**2))
    print '(a,2(es13.5,1x))', 'real fixed-rank singular values: ',s
    print '(a,es13.5)', 'real reconstruction relative error: ',err
    if (err > 1.0e-9_dp) error stop 'unexpected real reconstruction error'

    do j = 1,5
        z(:,j) = cmplx(a(:,j),a(:,j)*0.5_dp,dp)
    end do
    call interpolative_svd(z,1.0e-8_dp,k,uz,sz,vz,status)
    if (status /= linalg_status_success) error stop 'adaptive complex SVD failed'
    err = sqrt(sum(abs(z-matmul(uz,matmul(diagonal(sz),transpose(conjg(vz)))))**2)) / &
        sqrt(sum(abs(z)**2))
    print '(a,i0)', 'complex adaptive selected rank: ',k
    print '(a,es13.5)', 'complex reconstruction relative error: ',err
    if (k /= 2 .or. err > 1.0e-9_dp) error stop 'unexpected complex SVD result'

contains

    function diagonal(x) result(d)
        real(dp), intent(in) :: x(:) !! singular values
        real(dp) :: d(size(x),size(x))
        integer :: t
        d = 0.0_dp
        do t=1,size(x)
            d(t,t) = x(t)
        end do
    end function diagonal
end program linalg_interpolative_svd_example
