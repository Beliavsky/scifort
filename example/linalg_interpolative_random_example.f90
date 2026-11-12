! SPDX-License-Identifier: MIT
program linalg_interpolative_random_example
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success
    use scifort_linalg_interpolative_random, only : randomized_interp_decomp
    use scifort_linalg_interpolative, only : reconstruct_matrix_from_id
    implicit none
    real(dp) :: a(40,30), left(40,3), right(3,30), reconstructed(40,30)
    real(dp), allocatable :: proj(:, :)
    complex(dp) :: z(40,30), reconstructed_z(40,30)
    complex(dp), allocatable :: projz(:, :)
    integer, allocatable :: idx(:)
    real(dp) :: error
    integer :: i,j,k,status

    do i=1,40
        left(i,1) = sin(real(i,dp))
        left(i,2) = cos(real(3*i,dp))
        left(i,3) = sin(real(7*i,dp))
    end do
    do j=1,30
        right(1,j) = sin(real(2*j,dp))
        right(2,j) = cos(real(5*j,dp))
        right(3,j) = sin(real(11*j,dp))
    end do
    a = matmul(left,right)
    call randomized_interp_decomp(a,3,idx,proj,status,seed=12345_int64)
    if (status /= linalg_status_success) error stop 'randomized fixed-rank ID failed'
    call reconstruct_matrix_from_id(a(:,idx(1:3)),idx,proj,reconstructed,status)
    if (status /= linalg_status_success) error stop 'ID reconstruction failed'
    error = sqrt(sum((a-reconstructed)**2))/sqrt(sum(a**2))
    print '(a,3(i0,1x))', 'first three original columns chosen: ',idx(1:3)
    print '(a,es13.5)', 'real fixed-rank ID relative error: ',error
    if (error > 1.0e-9_dp) error stop 'unexpected reconstruction error'

    z = cmplx(a,0.5_dp*a,dp)
    call randomized_interp_decomp(z,1.0e-8_dp,k,idx,projz,status,seed=12345_int64)
    if (status /= linalg_status_success) error stop 'randomized adaptive ID failed'
    call reconstruct_matrix_from_id(z(:,idx(1:k)),idx,projz,reconstructed_z,status)
    if (status /= linalg_status_success) error stop 'complex ID reconstruction failed'
    error = sqrt(sum(abs(z-reconstructed_z)**2))/sqrt(sum(abs(z)**2))
    print '(a,i0)', 'complex adaptive rank: ',k
    print '(a,es13.5)', 'complex adaptive relative error: ',error
    if (k /= 3 .or. error > 1.0e-8_dp) error stop 'unexpected complex rank or error'
end program linalg_interpolative_random_example
