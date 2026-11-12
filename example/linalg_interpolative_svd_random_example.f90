! SPDX-License-Identifier: MIT
program linalg_interpolative_svd_random_example
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success
    use scifort_linalg_interpolative_svd_random, only : randomized_interpolative_svd
    implicit none
    real(dp) :: a(40,30), left(40,3), right(3,30), u(40,3), s(3), v(30,3)
    complex(dp) :: z(40,30)
    complex(dp), allocatable :: uz(:, :), vz(:, :)
    real(dp), allocatable :: sz(:)
    real(dp) :: err
    integer :: i, j, k, status

    do i = 1,40
        left(i,1) = sin(real(i,dp))
        left(i,2) = cos(real(3*i,dp))
        left(i,3) = sin(real(7*i,dp))
    end do
    do j = 1,30
        right(1,j) = sin(real(2*j,dp))
        right(2,j) = cos(real(5*j,dp))
        right(3,j) = sin(real(11*j,dp))
    end do
    a = matmul(left,right)
    call randomized_interpolative_svd(a,3,u,s,v,status,seed=12345_int64)
    if (status /= linalg_status_success) error stop 'fixed-rank randomized SVD failed'
    err = sqrt(sum((a-matmul(u,matmul(diagonal(s),transpose(v))))**2))/sqrt(sum(a**2))
    print '(a,3(es13.5,1x))', 'real rank-three singular values: ',s
    print '(a,es13.5)', 'real reconstruction relative error: ',err
    if (err > 1.0e-10_dp) error stop 'unexpected real error'

    z = cmplx(a,0.5_dp*a,dp)
    call randomized_interpolative_svd(z,1.0e-8_dp,k,uz,sz,vz,status, &
        seed=12345_int64,oversampling=8)
    if (status /= linalg_status_success) error stop 'adaptive randomized SVD failed'
    err = sqrt(sum(abs(z-matmul(uz,matmul(diagonal(sz),transpose(conjg(vz)))))**2))/ &
        sqrt(sum(abs(z)**2))
    print '(a,i0)', 'complex adaptive selected rank: ',k
    print '(a,es13.5)', 'complex reconstruction relative error: ',err
    if (k /= 3 .or. err > 1.0e-8_dp) error stop 'unexpected complex rank or error'

contains
    function diagonal(values) result(x)
        real(dp), intent(in) :: values(:) !! singular values
        real(dp) :: x(size(values),size(values))
        integer :: t
        x = 0.0_dp
        do t = 1,size(values)
            x(t,t) = values(t)
        end do
    end function diagonal
end program linalg_interpolative_svd_random_example
