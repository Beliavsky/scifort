! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program linalg_block_diag_khatri_rao_example
    use scifort_kinds, only : dp
    use scifort_linalg_block_diag, only : as_block, block_diag
    use scifort_linalg_khatri_rao, only : khatri_rao
    implicit none

    integer :: status, row
    real(dp) :: a(2,2), b(1,3), diagonal(3,5)
    real(dp) :: u(2,3), v(3,3), product(6,3)
    real(dp) :: batches(2,2,2), batch_diagonal(2,3,3)
    complex(dp) :: z(1,2), complex_diag(3,4)

    a(1,:) = [1.0_dp, 2.0_dp]
    a(2,:) = [3.0_dp, 4.0_dp]
    b(1,:) = [5.0_dp, 6.0_dp, 7.0_dp]
    call block_diag([as_block(a), as_block(b)], diagonal, status)
    if (status /= 0) error stop 'block_diag failed'
    print *, 'block_diag([2x2,1x3]) has shape (3,5):'
    do row=1,size(diagonal,1)
        print '(5f7.1)', diagonal(row,:)
    end do

    u(1,:) = [1.0_dp, 2.0_dp, 3.0_dp]
    u(2,:) = [4.0_dp, 5.0_dp, 6.0_dp]
    v(1,:) = [3.0_dp, 4.0_dp, 5.0_dp]
    v(2,:) = [6.0_dp, 7.0_dp, 8.0_dp]
    v(3,:) = [2.0_dp, 3.0_dp, 9.0_dp]
    call khatri_rao(u,v,product,status)
    if (status /= 0) error stop 'khatri_rao failed'
    print *, 'columnwise Kronecker product has shape (6,3):'
    do row=1,size(product,1)
        print '(3f8.1)', product(row,:)
    end do

    batches(1,:,:) = a
    batches(2,:,:) = 2.0_dp*a
    call block_diag([as_block(batches), as_block(9.0_dp)], batch_diagonal, status)
    if (status /= 0) error stop 'broadcast block_diag failed'
    print *, 'second batch bottom-right:', batch_diagonal(2,3,3)

    z(1,:) = [cmplx(1.0_dp,2.0_dp,dp), cmplx(3.0_dp,-1.0_dp,dp)]
    call block_diag([as_block(a), as_block(z)], complex_diag, status)
    if (status /= 0) error stop 'mixed-type block_diag failed'
    print *, 'complex block sample:', complex_diag(3,3)
end program linalg_block_diag_khatri_rao_example
