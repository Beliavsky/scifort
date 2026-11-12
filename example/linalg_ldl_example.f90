! SPDX-License-Identifier: MIT
program linalg_ldl_example
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success
    use scifort_linalg_ldl, only : ldl
    implicit none
    real(dp) :: a(3,3), lu(3,3), d(3,3), recon(3,3)
    integer :: perm(3), status

    ! A has zero diagonal entries but is nonsingular and indefinite.
    a(1,:) = [0.0_dp, 2.0_dp, 1.0_dp]
    a(2,:) = [2.0_dp, 0.0_dp, 0.5_dp]
    a(3,:) = [1.0_dp, 0.5_dp, 3.0_dp]

    call ldl(a,lu,d,perm,status)
    if (status /= linalg_status_success) error stop 'LDL factorization failed'

    recon = matmul(lu,matmul(d,transpose(lu)))
    print '(a,3(i0,1x))', 'One-based permutation: ', perm
    print '(a,es13.5)', 'Maximum reconstruction error: ', maxval(abs(recon-a))
    print '(a)', 'Block diagonal factor D:'
    print '(3(f10.4,1x))', transpose(d)
end program linalg_ldl_example
