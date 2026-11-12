! SPDX-License-Identifier: MIT
program linalg_ordqz_example
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success
    use scifort_linalg_ordqz, only : ordqz, ordqz_factors
    implicit none
    real(dp) :: a(4,4), b(4,4)
    complex(dp) :: aa(4,4), bb(4,4), q(4,4), z(4,4), alpha(4), beta(4)
    integer :: status, nselected, i

    a = 0.0_dp
    b = 0.0_dp
    do i=1,4
        b(i,i) = 1.0_dp
    end do
    a(1,1) = 2.0_dp
    a(2,2) = -3.0_dp
    a(3,3) = 0.5_dp
    a(4,4) = -1.0_dp
    a(1,3) = 0.2_dp
    a(2,4) = 0.3_dp
    call ordqz(a,b,aa,bb,q,z,status,alpha,beta,sort='lhp',sdim=nselected)
    if (status /= linalg_status_success) error stop 'ordqz failed'
    print '(a,i0)', 'Left-half-plane eigenvalues: ',nselected
    do i=1,4
        print '(a,i0,a,2(es15.6,1x))', 'lambda(',i,') = ', &
            real(alpha(i)/beta(i),dp),aimag(alpha(i)/beta(i))
    end do

    ! Reorder the computed factors without recomputing the generalized Schur form.
    call ordqz_factors(aa,bb,q,z,status,sort='iuc',sdim=nselected)
    if (status /= linalg_status_success) error stop 'factor reorder failed'
    print '(a,i0)', 'Eigenvalues inside the unit circle: ',nselected
end program linalg_ordqz_example
