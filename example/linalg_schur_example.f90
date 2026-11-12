! SPDX-License-Identifier: MIT
program linalg_schur_example
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success
    use scifort_linalg_schur, only : schur
    implicit none

    real(dp) :: a(3,3)
    complex(dp) :: t(3,3), z(3,3)
    real(dp) :: rt(3,3), rz(3,3)
    integer :: status, sdim, i

    ! Real input can have complex eigenvalues; the Schur output is complex.
    a(1,:) = [-2.0_dp, -1.0_dp, 0.0_dp]
    a(2,:) = [ 1.0_dp, -2.0_dp, 0.0_dp]
    a(3,:) = [ 0.2_dp,  0.1_dp, 3.0_dp]
    call schur(a, t, z, status, sort='lhp', sdim=sdim)
    if (status /= linalg_status_success) error stop 'Schur decomposition failed'
    if (sdim /= 2) error stop 'Unexpected selected invariant subspace size'
    print '(a,i0)', 'Eigenvalues selected in left half plane: ', sdim
    print '(a)', 'Complex Schur diagonal:'
    do i=1,3
        print '(2(es18.10,1x))', real(t(i,i),dp),aimag(t(i,i))
    end do
    print '(a,es12.4)', 'Max complex-Schur reconstruction error: ', &
        maxval(abs(matmul(matmul(z,t),transpose(conjg(z)))-cmplx(a,kind=dp)))

    call schur(a, rt, rz, status, sort='lhp', sdim=sdim)
    if (status /= linalg_status_success .or. sdim /= 2) error stop 'Real Schur failed'
    print '(a)', 'Real Schur form (2x2 block for complex pair):'
    do i = 1,3
        print '(3(f12.6,1x))', rt(i,:)
    end do
    print '(a,es12.4)', 'Max real-Schur reconstruction error: ', &
        maxval(abs(matmul(matmul(rz,rt),transpose(rz))-a))
end program linalg_schur_example
