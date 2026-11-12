! SPDX-License-Identifier: MIT
program linalg_banded_bulge_example
    use scifort_kinds, only : dp
    use scifort_linalg_banded_tridiagonal, only : band_to_tridiagonal
    use scifort_linalg_banded_eigen, only : eig_banded, eigvals_banded
    implicit none
    integer, parameter :: n = 6, kd = 2
    real(dp) :: ab(kd+1,n), d(n), e(n-1), q(n,n), w(n), z(n,n), selected(2)
    integer :: j, status, found

    ab = 0.0_dp
    do j = 1,n
        ab(1,j) = 3.0_dp + 0.1_dp*real(j,dp)
        if (j < n) ab(2,j) = -0.4_dp
        if (j < n-1) ab(3,j) = 0.2_dp
    end do

    call band_to_tridiagonal(ab,d,e,status,q,lower=.true.)
    if (status /= 0) error stop 'banded reduction failed'
    print '(a)', 'Band-to-tridiagonal diagonal:'
    print '(6f12.7)', d
    print '(a)', 'Subdiagonal:'
    print '(5f12.7)', e

    call eig_banded(ab,w,status,z,lower=.true.,method='bulge')
    if (status /= 0) error stop 'full banded eigensystem failed'
    print '(a)', 'Eigenvalues (Givens bulge chasing + tridiagonal QL):'
    print '(6f12.7)', w

    call eigvals_banded(ab,selected,status,lower=.true.,method='sturm', &
        subset_by_index=[3,4],num_selected=found)
    if (status /= 0 .or. found /= 2) error stop 'selected spectrum failed'
    print '(a)', 'Two selected eigenvalues (banded Sturm):'
    print '(2f12.7)', selected
end program linalg_banded_bulge_example
