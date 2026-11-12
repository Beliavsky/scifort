! SPDX-License-Identifier: MIT
program linalg_lu_det_inv_example
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success
    use scifort_linalg_lu, only : lu, lu_indices, lu_permuted
    use scifort_linalg_det_inv, only : det, inv
    implicit none
    real(dp) :: a(3,2), p(3,3), l(3,2), u(2,2), pl(3,2)
    real(dp) :: b(2,2), b_inv(2,2), determinant
    complex(dp) :: z(2,2), z_inv(2,2), z_det
    integer :: perm(3), status

    a(1,:) = [2.0_dp, 7.0_dp]
    a(2,:) = [5.0_dp, 1.0_dp]
    a(3,:) = [3.0_dp, 4.0_dp]
    call lu(a,p,l,u,status)
    if (status /= linalg_status_success) error stop 'lu example'
    print '(a,es12.4)', 'Rectangular LU residual: ', maxval(abs(a-matmul(p,matmul(l,u))))
    call lu_indices(a,perm,l,u,status)
    if (status /= linalg_status_success) error stop 'lu_indices example'
    print '(a,3(1x,i0))', 'One-based row permutation:', perm
    call lu_permuted(a,pl,u,status)
    if (status /= linalg_status_success) error stop 'lu_permuted example'
    print '(a,es12.4)', 'Permuted-L residual:   ', maxval(abs(a-matmul(pl,u)))

    b(1,:) = [1.0_dp, 2.0_dp]
    b(2,:) = [3.0_dp, 4.0_dp]
    call det(b,determinant,status)
    if (status /= linalg_status_success) error stop 'det example'
    print '(a,f8.3)', 'Real determinant: ', determinant
    call inv(b,b_inv,status)
    if (status /= linalg_status_success) error stop 'inv example'
    print '(a,es12.4)', 'Real inverse residual: ', &
        maxval(abs(matmul(b,b_inv)-reshape([1.0_dp,0.0_dp,0.0_dp,1.0_dp],[2,2])))

    z(1,:) = [cmplx(1.0_dp,1.0_dp,dp),cmplx(2.0_dp,0.0_dp,dp)]
    z(2,:) = [cmplx(3.0_dp,0.0_dp,dp),cmplx(2.0_dp,-1.0_dp,dp)]
    call det(z,z_det,status)
    if (status /= linalg_status_success) error stop 'complex det example'
    print '(a,2(1x,f8.3))', 'Complex determinant (real, imag):', real(z_det,dp), aimag(z_det)
    call inv(z,z_inv,status)
    if (status /= linalg_status_success) error stop 'complex inv example'
    print '(a,es12.4)', 'Complex inverse residual: ', &
        maxval(abs(matmul(z,z_inv)-reshape([cmplx(1.0_dp,0.0_dp,dp), &
        cmplx(0.0_dp,0.0_dp,dp),cmplx(0.0_dp,0.0_dp,dp), &
        cmplx(1.0_dp,0.0_dp,dp)],[2,2])))
end program linalg_lu_det_inv_example
