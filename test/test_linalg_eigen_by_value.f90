! SPDX-License-Identifier: MIT
! Native regression checks for half-open eigenvalue intervals and capacity outputs.
program test_linalg_eigen_by_value
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, &
        ieee_negative_inf, ieee_positive_inf
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape
    use scifort_linalg_dense_eigen, only : eigh, eigvalsh
    implicit none

    real(dp) :: a(4,4), b(4,4), values(4), vectors(4,4)
    real(dp) :: small(1), narrow(2), znarrow(4,1), nan, low_inf, high_inf
    real(dp) :: tiny_w(0), tiny_z(4,0), empty_matrix(0,0), empty_w(0)
    complex(dp) :: ca(4,4), cb(4,4), cv(4,4)
    integer :: status, nfound, i, typ

    a = 0.0_dp
    b = 0.0_dp
    do i = 1, 4
        b(i,i) = 1.0_dp
    end do
    a(1,1) = 1.0_dp
    a(2,2) = 2.0_dp
    a(3,3) = 2.0_dp
    a(4,4) = 3.0_dp
    low_inf = ieee_value(0.0_dp, ieee_negative_inf)
    high_inf = ieee_value(0.0_dp, ieee_positive_inf)
    nan = ieee_value(0.0_dp, ieee_quiet_nan)

    call eigh(a, values, status, vectors, subset_by_value=[1.0_dp,2.0_dp], &
        num_selected=nfound)
    call check(status == linalg_status_success, 'standard interval status')
    call check(nfound == 2, 'half-open interval and repeated eigenvalues')
    call check(all(abs(values(1:2) - 2.0_dp) < 1.0e-13_dp), 'selected eigenvalues')
    call check(all(values(3:4) == 0.0_dp), 'unused output initialized')
    call check(all(vectors(:,3:4) == 0.0_dp), 'unused eigenvectors initialized')
    call check(maxval(abs(matmul(a,vectors(:,1:2))-2.0_dp*vectors(:,1:2))) &
        < 1.0e-13_dp, 'selected eigenvector residual')

    call eigvalsh(a, values, status, subset_by_value=[low_inf,high_inf], &
        num_selected=nfound)
    call check(status == linalg_status_success .and. nfound == 4, 'unbounded interval')
    call check(maxval(abs(values - [1.0_dp,2.0_dp,2.0_dp,3.0_dp])) &
        < 1.0e-13_dp, 'unbounded spectrum')
    call eigvalsh(a, values, status, subset_by_value=[3.0_dp,4.0_dp], &
        num_selected=nfound)
    call check(status == linalg_status_success .and. nfound == 0, 'empty selection')
    call check(all(values == 0.0_dp), 'empty selection zeros')
    call eigh(a, tiny_w, status, tiny_z, subset_by_value=[3.0_dp,4.0_dp], &
        num_selected=nfound)
    call check(status == linalg_status_success .and. nfound == 0, 'empty buffers')

    call eigvalsh(a, small, status, subset_by_value=[1.0_dp,2.0_dp], &
        num_selected=nfound)
    call check(status == linalg_status_invalid_shape .and. nfound == 2, &
        'insufficient eigenvalue capacity reports count')
    call eigh(a, narrow, status, znarrow, subset_by_value=[1.0_dp,2.0_dp], &
        num_selected=nfound)
    call check(status == linalg_status_invalid_shape .and. nfound == 2, &
        'insufficient eigenvector capacity reports count')
    call eigvalsh(a, values, status, subset_by_value=[1.0_dp,1.0_dp])
    call check(status == linalg_status_invalid_shape, 'equal bounds rejected')
    call eigvalsh(a, values, status, subset_by_value=[3.0_dp,1.0_dp])
    call check(status == linalg_status_invalid_shape, 'reversed bounds rejected')
    call eigvalsh(a, values, status, subset_by_value=[nan,2.0_dp])
    call check(status == linalg_status_invalid_shape, 'nan bounds rejected')
    call eigvalsh(a, values, status, subset_by_value=[1.0_dp,2.0_dp], &
        subset_by_index=[1,2])
    call check(status == linalg_status_invalid_shape, 'conflicting selectors rejected')
    call eigvalsh(a, values, status, num_selected=nfound)
    call check(status == linalg_status_success .and. nfound == 4, 'legacy default count')
    call eigvalsh(a, narrow, status, subset_by_index=[2,3], num_selected=nfound)
    call check(status == linalg_status_success .and. nfound == 2, 'legacy index count')
    call eigvalsh(empty_matrix, empty_w, status, &
        subset_by_value=[low_inf,high_inf], num_selected=nfound)
    call check(status == linalg_status_success .and. nfound == 0, 'zero-order matrix')

    ! Generalized problem type 1 with exactly repeated generalized eigenvalues.
    b(2,2) = 2.0_dp
    call eigh(a, values, status, vectors, b=b, subset_by_value=[1.0_dp,2.0_dp], &
        num_selected=nfound)
    call check(status == linalg_status_success .and. nfound == 1, 'type 1 selection')
    call check(abs(values(1) - 2.0_dp) < 1.0e-13_dp, 'type 1 selected value')
    call check(maxval(abs(matmul(a,vectors(:,1:1)) - &
        matmul(b,vectors(:,1:1))*values(1))) < 1.0e-13_dp, 'type 1 residual')
    call eigh(a, values, status, vectors, b=b, subset_by_value=[high_inf,high_inf], &
        num_selected=nfound)
    call check(status == linalg_status_invalid_shape, 'generalized invalid interval')
    call eigh(a, tiny_w, status, tiny_z, b=b, &
        subset_by_value=[4.0_dp,5.0_dp], num_selected=nfound)
    call check(status == linalg_status_success .and. nfound == 0, &
        'generalized no matches with zero-sized output')
    call eigvalsh(a, small, status, b=b, &
        subset_by_value=[low_inf,high_inf], num_selected=nfound)
    call check(status == linalg_status_invalid_shape .and. nfound == 4, &
        'generalized insufficient capacity reports count')
    do typ = 2, 3
        call eigh(a, values, status, vectors, b=b, problem_type=typ, &
            subset_by_value=[2.1_dp,7.0_dp], num_selected=nfound)
        call check(status == linalg_status_success, 'generalized type 2 or 3 status')
        call check(nfound >= 1, 'generalized type 2 or 3 selected count')
    end do

    ca = cmplx(a,0.0_dp,dp)
    cb = cmplx(b,0.0_dp,dp)
    ca(2,1) = cmplx(0.0_dp,0.2_dp,dp)
    ca(1,2) = conjg(ca(2,1))
    call eigh(ca, values, status, cv, b=cb, subset_by_value=[-10.0_dp,3.5_dp], &
        num_selected=nfound)
    call check(status == linalg_status_success .and. nfound >= 1, &
        'complex generalized selection')
    call check(all(abs(aimag(cv(:,nfound+1:))) == 0.0_dp), &
        'complex unused columns zero')

    print '(a)', 'eigenvalue interval native tests passed'
contains
    subroutine check(condition, message)
        logical, intent(in) :: condition !! test outcome
        character(*), intent(in) :: message !! test identification
        if (.not. condition) then
            print '(a)', 'FAILED: '//message
            error stop 1
        end if
    end subroutine check
end program test_linalg_eigen_by_value
