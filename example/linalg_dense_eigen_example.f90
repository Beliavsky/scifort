! SPDX-License-Identifier: MIT
program linalg_dense_eigen_example
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success
    use scifort_linalg_dense_eigen, only : eigh, eigvalsh
    implicit none

    real(dp) :: a(3,3), w(3), vectors(3,3), largest(2)
    complex(dp) :: c(2,2), cv(2,2)
    real(dp) :: cw(2)
    integer :: status, i

    a = 0.0_dp
    do i=1,3
        a(i,i) = 2.0_dp
    end do
    a(2,1) = -1.0_dp
    a(3,2) = -1.0_dp
    ! The default is lower triangular storage; upper entries need not be set.
    call eigh(a, w, status, vectors)
    if (status /= linalg_status_success) error stop 'real eigh failed'
    print '(a,3f14.8)', 'Real symmetric eigenvalues:', w

    ! Fortran uses one-based inclusive indices. Return eigenvalues 2 and 3.
    call eigvalsh(a, largest, status, subset_by_index=[2,3])
    if (status /= linalg_status_success) error stop 'selected eigvalsh failed'
    print '(a,2f14.8)', 'Two largest eigenvalues:', largest

    c = cmplx(0.0_dp,0.0_dp,dp)
    c(1,1) = cmplx(2.0_dp,0.0_dp,dp)
    c(2,2) = cmplx(2.0_dp,0.0_dp,dp)
    c(2,1) = cmplx(0.0_dp,1.0_dp,dp)
    call eigh(c,cw,status,cv,lower=.true.)
    if (status /= linalg_status_success) error stop 'Hermitian eigh failed'
    print '(a,2f14.8)', 'Complex Hermitian eigenvalues:', cw
end program linalg_dense_eigen_example
