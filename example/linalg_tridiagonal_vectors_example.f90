! SPDX-License-Identifier: MIT
program linalg_tridiagonal_vectors_example
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success
    use scifort_linalg_banded_eigen, only : eigh_tridiagonal
    implicit none

    integer, parameter :: n = 2000, k = 3
    real(dp) :: d(n), e(n-1), eigenvalues(k), vectors(n,k)
    real(dp) :: residual(n), max_residual
    integer :: status, nfound, j

    d = 2.0_dp
    e = -1.0_dp
    call eigh_tridiagonal(d,e,eigenvalues,status,vectors, &
        subset_by_index=[998,1000], num_selected=nfound, &
        eigenvector_method='inverse')
    if (status /= linalg_status_success .or. nfound /= k) error stop 'eigh failed'

    max_residual = 0.0_dp
    do j = 1, k
        residual = (d-eigenvalues(j))*vectors(:,j)
        residual(1:n-1) = residual(1:n-1)+e*vectors(2:n,j)
        residual(2:n) = residual(2:n)+e*vectors(1:n-1,j)
        max_residual = max(max_residual,sqrt(sum(residual*residual)))
    end do
    print '(a,3(1x,f16.12))', 'selected eigenvalues:', eigenvalues
    print '(a,es12.4)', 'largest residual: ',max_residual
end program linalg_tridiagonal_vectors_example
