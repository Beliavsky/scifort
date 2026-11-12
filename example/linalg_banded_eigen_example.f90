! SPDX-License-Identifier: MIT
program linalg_banded_eigen_example
    use scifort_kinds, only : dp
    use scifort_linalg_banded_eigen, only : eigh_tridiagonal, eig_banded
    implicit none
    integer :: status, i
    real(dp) :: diagonal(4), offdiag(3), eig(4), q(4,4), b(2,4)
    complex(dp) :: cb(2,4), cq(4,4)

    diagonal = 2.0_dp
    offdiag = -1.0_dp
    call eigh_tridiagonal(diagonal,offdiag,eig,status,q)
    if (status /= 0) error stop 'eigh_tridiagonal failed'
    print '(a)', 'eigh_tridiagonal eigenvalues:'
    print '(4f14.8)',eig

    ! Upper-band storage: row kd+1 stores the diagonal.
    b = 0.0_dp
    b(2,:) = diagonal
    b(1,2:4) = offdiag
    call eig_banded(b,eig,status,q)
    if (status /= 0) error stop 'real eig_banded failed'
    print '(a)', 'real eig_banded eigenvalues:'
    print '(4f14.8)',eig

    ! Lower-band storage: row 1 stores the real diagonal;
    ! the next row is the complex first subdiagonal.
    cb = cmplx(0.0_dp,0.0_dp,dp)
    cb(1,:) = cmplx(2.0_dp,0.0_dp,dp)
    cb(2,1:3) = cmplx(0.0_dp,1.0_dp,dp)
    call eig_banded(cb,eig,status,cq,lower=.true.)
    if (status /= 0) error stop 'complex eig_banded failed'
    print '(a)', 'complex Hermitian eig_banded eigenvalues:'
    print '(4f14.8)',eig
    print '(a)', 'first complex eigenvector (real, imag):'
    do i=1,4
        print '(2f15.8)',real(cq(i,1),dp),aimag(cq(i,1))
    end do
end program linalg_banded_eigen_example
