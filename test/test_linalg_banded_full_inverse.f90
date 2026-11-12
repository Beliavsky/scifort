! SPDX-License-Identifier: MIT
program test_linalg_banded_full_inverse
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape
    use scifort_linalg_banded_eigen, only : eig_banded, eigvals_banded
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    implicit none

    integer, parameter :: n=12
    real(dp) :: ar(n,n), abr(3,n), wr(n), zr(n,n), wd(n), zd(n,n)
    real(dp) :: wrshort(n-1), zshort(n,n-1), nonfinite(3,n), lam, err
    complex(dp) :: ac(n,n), abc(3,n), zc(n,n)
    integer :: status, found, i, j, k
    real(dp), allocatable :: emptya(:, :), emptyw(:), emptyz(:, :)

    ar = 0.0_dp
    abr = 0.0_dp
    do j=1,n
        ar(j,j) = 3.0_dp+0.12_dp*real(j,dp)
        abr(1,j) = ar(j,j)
        do k=1,min(2,n-j)
            lam = 0.15_dp/real(k,dp)
            ar(j+k,j)=lam
            ar(j,j+k)=lam
            abr(k+1,j)=lam
        end do
    end do
    call eig_banded(abr,wr,status,zr,lower=.true.,method='inverse',num_selected=found)
    call check(status == linalg_status_success .and. found == n, 'real full inverse')
    call eig_banded(abr,wd,status,zd,lower=.true.,method='dense')
    call check(status == linalg_status_success, 'real dense')
    call check(maxval(abs(wr-wd)) < 1.0e-11_dp, 'real eigenvalues')
    call check(maxval(abs(matmul(ar,zr)-matmul(zr,diagonal(wr)))) < 2.0e-10_dp, &
        'real eigenpair residual')
    call check(maxval(abs(matmul(transpose(zr),zr)-identity())) < 1.0e-11_dp, &
        'real orthogonality')

    ac = cmplx(ar,0.0_dp,dp)
    abc = cmplx(abr,0.0_dp,dp)
    do j=1,n-1
        ac(j+1,j)=cmplx(0.0_dp,0.15_dp,dp)
        ac(j,j+1)=conjg(ac(j+1,j))
        abc(2,j)=ac(j+1,j)
    end do
    call eig_banded(abc,wr,status,zc,lower=.true.,method='inverse',num_selected=found)
    call check(status == linalg_status_success .and. found == n, 'complex full inverse')
    call check(maxval(abs(matmul(ac,zc)-matmul(zc,cdiagonal(wr)))) < 2.0e-10_dp, &
        'complex eigenpair residual')
    call check(maxval(abs(matmul(conjg(transpose(zc)),zc)-cidentity())) < 1.0e-11_dp, &
        'complex orthogonality')

    ! Upper packed complex input and values-only fallback to banded Sturm.
    abc = cmplx(0.0_dp,0.0_dp,dp)
    do j=1,n
        abc(3,j)=ac(j,j)
        do k=1,min(2,j-1)
            abc(3-k,j)=ac(j-k,j)
        end do
    end do
    call eig_banded(abc,wr,status,lower=.false.,method='inverse',num_selected=found)
    call check(status == linalg_status_success .and. found == n, 'values-only inverse')
    call eigvals_banded(abc,wd,status,lower=.false.,method='dense')
    call check(status == linalg_status_success, 'complex reference values')
    call check(maxval(abs(wr-wd)) < 1.0e-11_dp, 'complex eigenvalues')

    call eig_banded(abr,wr,status,zshort,lower=.true.,method='inverse', &
        num_selected=found)
    call check(status == linalg_status_invalid_shape .and. found == n, &
        'short full-eigenvector storage')
    call eig_banded(abr,wrshort,status,lower=.true.,method='inverse', &
        num_selected=found)
    call check(status == linalg_status_invalid_shape .and. found == n, &
        'short full-spectrum eigenvalue storage')
    call eig_banded(abr,wr,status,zr,lower=.true.,method='unknown')
    call check(status == linalg_status_invalid_shape, 'invalid strategy')

    ! Unused cells of upper-packed input must not be inspected.
    nonfinite = 0.0_dp
    nonfinite(3,:) = 2.0_dp
    nonfinite(2,2:n) = -0.25_dp
    nonfinite(1,:) = ieee_value(0.0_dp,ieee_quiet_nan)
    nonfinite(2,1) = ieee_value(0.0_dp,ieee_quiet_nan)
    nonfinite(1,3:n) = 0.0_dp
    call eig_banded(nonfinite,wr,status,zr,method='inverse')
    call check(status == linalg_status_success, 'ignore unused upper-packed cells')

    allocate(emptya(1,0),emptyw(0),emptyz(0,0))
    call eig_banded(emptya,emptyw,status,emptyz,method='inverse',num_selected=found)
    call check(status == linalg_status_success .and. found == 0, 'zero-sized full spectrum')
    print '(a)', 'banded full-spectrum inverse native tests passed'

contains
    subroutine check(ok, label)
        logical, intent(in) :: ok !! required condition
        character(len=*), intent(in) :: label !! test description
        if (.not. ok) then
            print '(a)', 'FAILED: '//label
            error stop 1
        end if
    end subroutine check

    pure function identity() result(q)
        real(dp) :: q(n,n)
        integer :: t
        q=0.0_dp
        do t=1,n
            q(t,t)=1.0_dp
        end do
    end function identity

    pure function cidentity() result(q)
        complex(dp) :: q(n,n)
        q=cmplx(identity(),0.0_dp,dp)
    end function cidentity

    pure function diagonal(w) result(q)
        real(dp), intent(in) :: w(:) !! eigenvalues
        real(dp) :: q(size(w),size(w))
        integer :: t
        q=0.0_dp
        do t=1,size(w)
            q(t,t)=w(t)
        end do
    end function diagonal

    pure function cdiagonal(w) result(q)
        real(dp), intent(in) :: w(:) !! eigenvalues
        complex(dp) :: q(size(w),size(w))
        q=cmplx(diagonal(w),0.0_dp,dp)
    end function cdiagonal
end program test_linalg_banded_full_inverse
