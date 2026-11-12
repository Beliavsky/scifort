! SPDX-License-Identifier: MIT
program test_linalg_tridiagonal_selection
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, &
        ieee_positive_inf, ieee_negative_inf
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    use scifort_linalg_banded_eigen, only : eigh_tridiagonal, eigvalsh_tridiagonal
    implicit none

    real(dp) :: d(5), e(4), full(5), w(5), vec(5,5), a(5,5), huge_bounds(2)
    real(dp) :: two(2), selected(3), subspace(5,3), empty(0), infp, infn
    real(dp), allocatable :: dbig(:), ebig(:)
    real(dp) :: pi, target
    integer :: status, count, i, j, n

    d = 2.0_dp
    e = -1.0_dp
    call eigvalsh_tridiagonal(d, e, full, status)
    call check(status == linalg_status_success, 'full QL')
    call eigvalsh_tridiagonal(d, e, selected, status, &
        subset_by_index=[2,4], num_selected=count)
    call check(status == linalg_status_success .and. count == 3, 'index count')
    call check(maxval(abs(selected-full(2:4))) < 3.0e-14_dp, 'Sturm bisection index')
    call eigh_tridiagonal(d, e, selected, status, subspace, &
        subset_by_index=[2,4], num_selected=count)
    call check(status == linalg_status_success .and. count == 3, 'selected vectors')
    call check(maxval(abs(selected-full(2:4))) < 3.0e-14_dp, 'selected vector values')
    a = 0.0_dp
    do i = 1, 5
        a(i,i) = 2.0_dp
        if (i < 5) then
            a(i,i+1) = -1.0_dp
            a(i+1,i) = -1.0_dp
        end if
    end do
    do j = 1, 3
        call check(maxval(abs(matmul(a,subspace(:,j))-selected(j)*subspace(:,j))) &
            < 1.0e-13_dp, 'selected vector residual')
    end do
    call check(maxval(abs(matmul(transpose(subspace),subspace) - identity(3))) &
        < 1.0e-13_dp, 'selected vector orthogonality')

    d = [1.0_dp,1.0_dp,2.0_dp,2.0_dp,3.0_dp]
    e = 0.0_dp
    call eigvalsh_tridiagonal(d,e,w,status,subset_by_value=[1.0_dp,2.0_dp], &
        num_selected=count)
    call check(status == linalg_status_success .and. count == 2, 'half-open count')
    call check(maxval(abs(w(1:2)-2.0_dp)) < 3.0e-14_dp, 'half-open repeated values')
    call check(all(w(3:) == 0.0_dp), 'unused values zeroed')
    call eigh_tridiagonal(d,e,w,status,vec,subset_by_value=[1.0_dp,2.0_dp], &
        num_selected=count)
    call check(status == linalg_status_success .and. count == 2, 'vectors half-open count')
    call check(maxval(abs(w(1:2)-2.0_dp)) < 3.0e-14_dp, 'vectors half-open values')
    call check(all(w(3:) == 0.0_dp), 'unused values zeroed for vectors')
    infp = ieee_value(0.0_dp,ieee_positive_inf)
    infn = ieee_value(0.0_dp,ieee_negative_inf)
    call eigvalsh_tridiagonal(d,e,w,status,subset_by_value=[infn,infp], &
        num_selected=count)
    call check(status == linalg_status_success .and. count == 5, 'infinite endpoints')
    call check(maxval(abs(w-d)) < 3.0e-14_dp, 'infinite interval values')
    call eigvalsh_tridiagonal(d,e,empty,status,subset_by_value=[5.0_dp,6.0_dp], &
        num_selected=count)
    call check(status == linalg_status_success .and. count == 0, 'empty selection')
    call eigvalsh_tridiagonal(d,e,empty,status,subset_by_value=[-6.0_dp,-5.0_dp], &
        num_selected=count)
    call check(status == linalg_status_success .and. count == 0, 'empty lower selection')
    call eigvalsh_tridiagonal(d,e,two,status,subset_by_value=[1.0_dp,3.0_dp], &
        num_selected=count)
    call check(status == linalg_status_invalid_shape .and. count == 3, 'small buffer count')
    call check(all(two == 0.0_dp), 'small buffer remains zero')
    call eigh_tridiagonal(d,e,w,status,vec(:,1:1),subset_by_value=[1.0_dp,3.0_dp], &
        num_selected=count)
    call check(status == linalg_status_invalid_shape .and. count == 3, 'small vector buffer')

    ! The diagonal shortcut preserves tiny exact eigenvalues even when other
    ! diagonal entries differ by hundreds of orders of magnitude.
    d = [1.0e-200_dp,-1.0e200_dp,0.0_dp,1.0e200_dp,1.0e-200_dp]
    call eigvalsh_tridiagonal(d,e,selected,status,subset_by_index=[2,4])
    call check(status == linalg_status_success, 'heterogeneous diagonal status')
    call check(all(selected == [0.0_dp,1.0e-200_dp,1.0e-200_dp]), &
        'heterogeneous diagonal exact values')
    call eigvalsh_tridiagonal(d,e,two,status,subset_by_value=[0.0_dp,1.0e-190_dp], &
        num_selected=count)
    call check(status == linalg_status_success .and. count == 2, 'tiny exact diagonal count')
    call check(all(two == 1.0e-200_dp), 'tiny exact diagonal values')

    d = 0.0_dp
    call eigvalsh_tridiagonal(d,e,w,status,subset_by_value=[-1.0_dp,0.0_dp], &
        num_selected=count)
    call check(status == linalg_status_success .and. count == 5, 'zero matrix selection')
    call eigvalsh_tridiagonal(d,e,empty,status,subset_by_value=[0.0_dp,1.0_dp], &
        num_selected=count)
    call check(status == linalg_status_success .and. count == 0, 'zero matrix exclusion')
    d = 2.0e150_dp
    e = -1.0e150_dp
    call eigvalsh_tridiagonal(d,e,selected,status,subset_by_index=[2,4])
    call check(status == linalg_status_success, 'large-scale bisection')
    call check(maxval(abs(selected/1.0e150_dp - [1.0_dp,2.0_dp,3.0_dp])) &
        < 3.0e-14_dp, 'large-scale values')
    d = 2.0e-150_dp
    e = -1.0e-150_dp
    call eigvalsh_tridiagonal(d,e,selected,status,subset_by_index=[2,4])
    call check(status == linalg_status_success, 'small-scale bisection')
    call check(maxval(abs(selected/1.0e-150_dp - [1.0_dp,2.0_dp,3.0_dp])) &
        < 3.0e-14_dp, 'small-scale values')
    huge_bounds = [1.0e290_dp,1.0e300_dp]
    call eigvalsh_tridiagonal(d,e,empty,status,subset_by_value=huge_bounds, &
        num_selected=count)
    call check(status == linalg_status_success .and. count == 0, 'safe scale-bound clamp')
    call eigvalsh_tridiagonal(d,e,w,status,subset_by_value=[-1.0e300_dp,1.0e300_dp], &
        num_selected=count)
    call check(status == linalg_status_success .and. count == 5, 'safe wide bounds')

    call eigvalsh_tridiagonal(d,e,w,status,subset_by_index=[0,2])
    call check(status == linalg_status_invalid_shape, 'invalid index')
    call eigvalsh_tridiagonal(d,e,w,status,subset_by_index=[4,3])
    call check(status == linalg_status_invalid_shape, 'reversed index')
    call eigvalsh_tridiagonal(d,e,w,status,subset_by_value=[0.0_dp,0.0_dp])
    call check(status == linalg_status_invalid_shape, 'empty invalid interval')
    call eigvalsh_tridiagonal(d,e,w,status,subset_by_value=[infn,0.0_dp], &
        subset_by_index=[1,3])
    call check(status == linalg_status_invalid_shape, 'conflicting selection modes')
    call eigvalsh_tridiagonal(d,e,w,status,subset_by_index=[1,3],tol=-1.0_dp)
    call check(status == linalg_status_invalid_shape, 'invalid tolerance')
    call eigvalsh_tridiagonal(d,e,w,status,subset_by_value= &
        [ieee_value(0.0_dp,ieee_quiet_nan),1.0_dp])
    call check(status == linalg_status_invalid_shape, 'nan interval')
    d = 2.0_dp
    e = -1.0_dp
    call eigvalsh_tridiagonal(d,e,selected,status,subset_by_index=[1,3],max_iterations=1)
    call check(status == linalg_status_no_convergence, 'limited bisection iterations')
    call check(all(selected == 0.0_dp), 'nonconverged result zero')

    n = 10000
    allocate(dbig(n), ebig(n-1))
    dbig = 2.0_dp
    ebig = -1.0_dp
    call eigvalsh_tridiagonal(dbig,ebig,w,status,subset_by_index=[5000,5004], &
        num_selected=count)
    call check(status == linalg_status_success .and. count == 5, 'order 10000 slice')
    pi = acos(-1.0_dp)
    do i = 1, 5
        target = 2.0_dp-2.0_dp*cos(real(4999+i,dp)*pi/real(n+1,dp))
        call check(abs(w(i)-target) < 2.0e-13_dp, 'order 10000 Toeplitz eigenvalue')
    end do

    print '(a)', 'tridiagonal selection native tests passed'
contains
    subroutine check(ok, label)
        logical, intent(in) :: ok !! asserted condition
        character(len=*), intent(in) :: label !! failure description
        if (.not. ok) then
            print '(a)', 'FAILED: '//label
            error stop 1
        end if
    end subroutine check

    pure function identity(k) result(x)
        integer, intent(in) :: k !! identity matrix order
        real(dp) :: x(k,k)
        integer :: j
        x = 0.0_dp
        do j = 1, k
            x(j,j) = 1.0_dp
        end do
    end function identity
end program test_linalg_tridiagonal_selection
