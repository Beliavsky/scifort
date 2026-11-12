! SPDX-License-Identifier: MIT
program test_linalg_tridiagonal_vectors
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape
    use scifort_linalg_banded_eigen, only : eigh_tridiagonal
    implicit none

    real(dp), allocatable :: d(:), e(:), z(:,:), r(:), w(:)
    real(dp) :: pi, expected, residual, orthogonality
    real(dp) :: ds(6), es(5), ws(4), zs(6,4), wql(4), zql(6,4)
    integer :: n, i, j, k, status, count

    n = 10000
    allocate(d(n),e(n-1),z(n,5),r(n),w(5))
    d = 2.0_dp
    e = -1.0_dp
    call eigh_tridiagonal(d,e,w,status,z,subset_by_index=[5000,5004], &
        num_selected=count,eigenvector_method='inverse')
    call check(status == linalg_status_success .and. count == 5, '10k narrow vectors')
    pi = acos(-1.0_dp)
    do j = 1, 5
        expected = 2.0_dp - 2.0_dp*cos(real(4999+j,dp)*pi/real(n+1,dp))
        call check(abs(w(j)-expected) < 2.0e-13_dp, '10k eigenvalues')
        r = d*z(:,j)-w(j)*z(:,j)
        r(1:n-1) = r(1:n-1)+e*z(2:n,j)
        r(2:n) = r(2:n)+e*z(1:n-1,j)
        residual = sqrt(sum(r*r))
        call check(residual < 2.0e-10_dp, '10k eigenvector residual')
        do k = 1, j
            orthogonality = dot_product(z(:,j),z(:,k))
            if (j == k) then
                call check(abs(orthogonality-1.0_dp) < 5.0e-12_dp, '10k normalized')
            else
                call check(abs(orthogonality) < 5.0e-10_dp, '10k orthogonal')
            end if
        end do
    end do

    ds = [3.0_dp,2.0_dp,3.0_dp,2.0_dp,1.0_dp,1.0_dp]
    es = 0.0_dp
    call eigh_tridiagonal(ds,es,ws,status,zs,subset_by_index=[2,5], &
        eigenvector_method='inverse')
    call check(status == linalg_status_success, 'diagonal repeated vectors')
    call check(all(ws == [1.0_dp,2.0_dp,2.0_dp,3.0_dp]), &
        'diagonal repeated eigenvalue ordering')
    do j = 1, 4
        call check(abs(dot_product(zs(:,j),zs(:,j))-1.0_dp) < 1.0e-15_dp, &
            'diagonal repeated norm')
        do k = 1, j-1
            call check(abs(dot_product(zs(:,j),zs(:,k))) < 1.0e-15_dp, &
                'diagonal repeated orthogonality')
        end do
    end do

    ds = 2.0_dp
    es = -1.0_dp
    call eigh_tridiagonal(ds,es,ws,status,zs,subset_by_index=[2,5], &
        eigenvector_method='inverse')
    call check(status == linalg_status_success, 'explicit inverse method')
    call eigh_tridiagonal(ds,es,wql,status,zql,subset_by_index=[2,5], &
        eigenvector_method='ql')
    call check(status == linalg_status_success, 'explicit legacy QL method')
    call check(maxval(abs(ws-wql)) < 3.0e-14_dp, 'inverse vs QL values')
    call eigh_tridiagonal(ds,es,ws,status,zs,subset_by_index=[2,5], &
        eigenvector_method='invalid')
    call check(status == linalg_status_invalid_shape, 'invalid method rejected')
    call eigh_tridiagonal(ds,es,ws,status,zs(:,1:2),subset_by_index=[2,5], &
        num_selected=count, eigenvector_method='inverse')
    call check(status == linalg_status_invalid_shape .and. count == 4, &
        'small vector buffer and selected count')

    print '(a)', 'tridiagonal selected eigenvector tests passed'
contains
    subroutine check(ok,label)
        logical, intent(in) :: ok !! required assertion
        character(len=*), intent(in) :: label !! assertion description
        if (.not. ok) then
            print '(a)', 'FAILED: '//label
            error stop 1
        end if
    end subroutine check
end program test_linalg_tridiagonal_vectors
