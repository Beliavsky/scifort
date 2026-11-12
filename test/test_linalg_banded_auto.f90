! SPDX-License-Identifier: MIT
program test_linalg_banded_auto
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape
    use scifort_linalg_banded_eigen, only : eigvals_banded, eig_banded
    implicit none
    integer, parameter :: n=208, kd=3
    real(dp) :: ab(kd+1,n), wauto(n), wbulge(n), wdense(n), weigen(n)
    real(dp) :: z(n,n), selected(4)
    complex(dp) :: cab(kd+1,n), cz(n,n)
    real(dp) :: residual
    integer :: i, j, k, st, count

    ab = 0.0_dp
    cab = cmplx(0.0_dp,0.0_dp,dp)
    do j=1,n
        ab(1,j) = 1.0_dp+0.6_dp*sin(0.37_dp*real(j,dp))
        cab(kd+1,j) = cmplx(ab(1,j),0.0_dp,dp)
        do k=1,min(kd,n-j)
            i=j+k
            ab(k+1,j) = 0.08_dp*cos(real(5*i+j,dp))
            cab(kd+1-k,i) = cmplx(ab(k+1,j),0.07_dp*sin(real(i+3*j,dp)),dp)
        end do
    end do

    call eigvals_banded(ab,wauto,st,lower=.true.)
    call require(st == linalg_status_success,'real auto eigvals')
    call eigvals_banded(ab,wbulge,st,lower=.true.,method='bulge')
    call require(st == linalg_status_success,'real bulge eigvals')
    call require(all(wauto == wbulge),'real auto selects bulge for narrow bands')
    call eigvals_banded(ab,wdense,st,lower=.true.,method='dense')
    call require(st == linalg_status_success,'real dense eigvals')
    call require(maxval(abs(wauto-wdense)) < 2.0e-11_dp,'real spectra agree')
    call eig_banded(ab,weigen,st,lower=.true.)
    call require(st == linalg_status_success,'real eig_banded values-only auto')
    call require(all(wauto == weigen),'real eig_banded auto selects bulge')

    call eig_banded(ab,weigen,st,z,lower=.true.)
    call require(st == linalg_status_success,'real auto with eigenvectors')
    call require(maxval(abs(weigen-wdense)) < 2.0e-11_dp,'real dense vector default')
    call eigvals_banded(ab,selected,st,lower=.true.,subset_by_index=[100,103], &
        num_selected=count)
    call require(st == linalg_status_success .and. count == 4,'selected auto unchanged')
    call require(maxval(abs(selected-wdense(100:103))) < 2.0e-8_dp, &
        'selected Sturm agrees with full values')
    call eigvals_banded(ab,wauto(:n-1),st,lower=.true.)
    call require(st == linalg_status_invalid_shape,'short full eigvals rejected')

    call eigvals_banded(cab,wauto,st)
    call require(st == linalg_status_success,'complex auto eigvals')
    call eigvals_banded(cab,wbulge,st,method='bulge')
    call require(st == linalg_status_success,'complex bulge eigvals')
    call require(all(wauto == wbulge),'complex auto selects bulge')
    call eigvals_banded(cab,wdense,st,method='dense')
    call require(st == linalg_status_success,'complex dense eigvals')
    call require(maxval(abs(wauto-wdense)) < 2.0e-11_dp,'complex spectra agree')
    call eig_banded(cab,weigen,st)
    call require(st == linalg_status_success,'complex eig_banded values-only auto')
    call require(all(wauto == weigen),'complex eig_banded auto selects bulge')
    call eig_banded(cab,weigen,st,cz)
    call require(st == linalg_status_success,'complex auto with eigenvectors')
    call require(maxval(abs(weigen-wdense)) < 2.0e-11_dp,'complex dense vector default')
    residual = 0.0_dp
    do j=1,n
        residual = max(residual,abs(sum(conjg(cz(:,j))*cz(:,j))-1.0_dp))
    end do
    call require(residual < 2.0e-11_dp,'complex eigenvector norms')
    print '(a)', 'banded auto eigenvalue native tests passed'

contains
    subroutine require(condition,label)
        logical, intent(in) :: condition !! required numerical assertion
        character(len=*), intent(in) :: label !! informative failed assertion label
        if (.not. condition) then
            print '(a)', 'FAILED: '//label
            error stop 1
        end if
    end subroutine require
end program test_linalg_banded_auto
