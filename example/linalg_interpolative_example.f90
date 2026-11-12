! SPDX-License-Identifier: MIT
program linalg_interpolative_example
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success
    use scifort_linalg_interpolative, only : interp_decomp, &
        reconstruct_skel_matrix, reconstruct_interp_matrix, &
        reconstruct_matrix_from_id, id_to_svd
    implicit none
    real(dp) :: a(4,5), b1(4), b2(4), error
    real(dp), allocatable :: b(:, :), p(:, :), approx(:, :)
    real(dp), allocatable :: u(:, :), s(:), v(:, :), proj(:, :)
    integer, allocatable :: idx(:)
    integer :: k, status

    b1 = [1.0_dp, 0.0_dp, 2.0_dp, 1.0_dp]
    b2 = [0.0_dp, 2.0_dp, 1.0_dp, -1.0_dp]
    a(:,1) = b1
    a(:,2) = b2
    a(:,3) = b1+b2
    a(:,4) = -b1+2.0_dp*b2
    a(:,5) = 3.0_dp*b1-b2

    ! Fixed-rank, deterministic interpolative decomposition.
    call interp_decomp(a,2,idx,proj,status)
    if (status /= linalg_status_success) error stop 'fixed-rank ID failed'
    k = size(proj,1)
    allocate(b(4,k),p(k,5),approx(4,5),u(4,k),s(k),v(5,k))
    call reconstruct_skel_matrix(a,k,idx,b,status)
    if (status /= linalg_status_success) error stop 'skeleton construction failed'
    call reconstruct_interp_matrix(idx,proj,p,status)
    if (status /= linalg_status_success) error stop 'interpolation construction failed'
    call reconstruct_matrix_from_id(b,idx,proj,approx,status)
    if (status /= linalg_status_success) error stop 'ID reconstruction failed'

    error = maxval(abs(approx-a))
    write(*,'(a,i0)') 'Fixed ID rank: ',k
    write(*,'(a,*(i0,1x))') 'One-based column order: ',idx
    write(*,'(a,es12.4)') 'Maximum ID reconstruction difference: ',error
    if (error > 1.0e-12_dp) error stop 'ID reconstruction inaccurate'

    call id_to_svd(b,idx,proj,u,s,v,status)
    if (status /= linalg_status_success) error stop 'ID to SVD failed'
    write(*,'(a,*(es12.4,1x))') 'Nonzero singular values: ',s

    ! Compute the numerical rank automatically from a Frobenius tolerance.
    call interp_decomp(a,1.0e-10_dp,k,idx,proj,status)
    if (status /= linalg_status_success) error stop 'tolerance ID failed'
    if (k /= 2) error stop 'expected numerical rank two'
    write(*,'(a,i0)') 'Tolerance-selected rank: ',k
end program linalg_interpolative_example
