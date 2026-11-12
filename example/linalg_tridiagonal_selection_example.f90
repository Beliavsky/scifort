! SPDX-License-Identifier: MIT
program linalg_tridiagonal_selection_example
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success
    use scifort_linalg_banded_eigen, only : eigh_tridiagonal, eigvalsh_tridiagonal
    implicit none

    real(dp) :: d(6), e(5), selected(3), by_value(6), vectors(6,3)
    integer :: status, count, i

    d = 2.0_dp
    e = -1.0_dp
    call eigvalsh_tridiagonal(d, e, selected, status, &
        subset_by_index=[2,4], num_selected=count)
    if (status /= linalg_status_success) error stop 'index selection failed'
    print '(a,i0)', 'Number of indexed eigenvalues: ', count
    print '(3(f16.10))', selected

    call eigvalsh_tridiagonal(d, e, by_value, status, &
        subset_by_value=[1.0_dp,3.0_dp], num_selected=count)
    if (status /= linalg_status_success) error stop 'value selection failed'
    print '(a,i0)', 'Eigenvalues in (1,3]: ', count
    print '(3(f16.10))', by_value(1:count)

    call eigh_tridiagonal(d, e, selected, status, vectors, subset_by_index=[2,4])
    if (status /= linalg_status_success) error stop 'selected vectors failed'
    print '(a)', 'Selected eigenvector column norms:'
    do i = 1, 3
        print '(f16.10)', sqrt(sum(vectors(:,i)**2))
    end do
end program linalg_tridiagonal_selection_example
