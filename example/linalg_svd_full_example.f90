! SPDX-License-Identifier: MIT
program linalg_svd_full_example
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success
    use scifort_linalg_svd, only : svd, svdvals, diagsvd
    implicit none
    real(dp) :: a(4,2), u(4,4), s(2), ss(2), vh(2,2), d(4,2)
    complex(dp) :: ac(2,4), uc(2,2), vhc(4,4)
    real(dp) :: sc(2)
    integer :: status

    a = reshape([1.0_dp,2.0_dp,0.0_dp,0.0_dp, &
                 1.0_dp,0.0_dp,3.0_dp,0.0_dp],shape(a))
    call svd(a,u,s,vh,status,full_matrices=.true.)
    if (status /= linalg_status_success) error stop 'real full SVD failed'
    call diagsvd(s,4,2,d,status)
    if (status /= linalg_status_success) error stop 'diagsvd failed'
    call svdvals(a,ss,status)
    if (status /= linalg_status_success) error stop 'svdvals failed'
    print '(A,ES12.4)', 'real 4x2 full SVD error: ', maxval(abs(matmul(matmul(u,d),vh)-a))
    print '(A,2(1X,ES12.4))', 'singular values: ', ss

    ac = cmplx(0.0_dp,0.0_dp,dp)
    ac(1,1) = cmplx(1.0_dp,1.0_dp,dp)
    ac(1,2) = cmplx(2.0_dp,-1.0_dp,dp)
    ac(2,3) = cmplx(0.0_dp,3.0_dp,dp)
    call svd(ac,uc,sc,vhc,status,full_matrices=.true.)
    if (status /= linalg_status_success) error stop 'complex full SVD failed'
    print '(A,ES12.4)', 'complex 2x4 full SVD error: ', &
        maxval(abs(matmul(uc*spread(sc,1,2),vhc(1:2,:))-ac))
    print '(A,2(I0,1X))', 'complete complex VH shape: ', shape(vhc)
end program linalg_svd_full_example
