! SPDX-License-Identifier: MIT
! Demonstrate (low, high] eigenvalue selection with fixed-capacity arrays.
program linalg_eigen_by_value_example
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success
    use scifort_linalg_dense_eigen, only : eigh, eigvalsh
    implicit none

    real(dp) :: a(4,4), b(4,4), w(4), z(4,4)
    complex(dp) :: ca(4,4), cb(4,4), cz(4,4)
    integer :: i, count, status

    a = 0.0_dp
    b = 0.0_dp
    do i = 1, 4
        b(i,i) = 1.0_dp
    end do
    a(1,1) = 1.0_dp
    a(2,2) = 2.0_dp
    a(3,3) = 2.0_dp
    a(4,4) = 4.0_dp

    ! Select both repeated eigenvalues equal to the inclusive upper bound.
    call eigh(a,w,status,z,subset_by_value=[1.0_dp,2.0_dp],num_selected=count)
    if (status /= linalg_status_success) error stop 'eigh by value failed'
    print '(a,i0)', 'standard: selected eigenvalues = ',count
    print '(a,*(f8.3,1x))', 'standard values: ',w(1:count)

    ! Type-1 generalized problem: A*x = lambda*B*x.
    b(2,2) = 2.0_dp
    call eigvalsh(a,w,status,b=b,subset_by_value=[1.0_dp,4.0_dp], &
        num_selected=count)
    if (status /= linalg_status_success) error stop 'generalized eigvalsh by value failed'
    print '(a,i0)', 'generalized: selected eigenvalues = ',count
    print '(a,*(f8.3,1x))', 'generalized values: ',w(1:count)

    ! The same API works for complex Hermitian-definite problems.
    ca = cmplx(a,0.0_dp,dp)
    cb = cmplx(b,0.0_dp,dp)
    ca(2,1) = cmplx(0.0_dp,0.1_dp,dp)
    ca(1,2) = conjg(ca(2,1))
    call eigh(ca,w,status,cz,b=cb,subset_by_value=[1.1_dp,4.5_dp], &
        num_selected=count)
    if (status /= linalg_status_success) error stop 'complex generalized eigh failed'
    print '(a,i0)', 'complex generalized: selected eigenvalues = ',count
    print '(a,*(f8.3,1x))', 'complex generalized values: ',w(1:count)
end program linalg_eigen_by_value_example
