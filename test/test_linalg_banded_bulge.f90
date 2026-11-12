! SPDX-License-Identifier: MIT
program test_linalg_banded_bulge
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape
    use scifort_linalg_banded_eigen, only : eig_banded, eigvals_banded
    implicit none
    integer, parameter :: n=12, kd=3
    real(dp) :: ab(kd+1,n), a(n,n), w(n), selected(4), vectors(n,n), sv(n,4), vals(4)
    real(dp) :: wscale(n), wref(n), scl, residual
    complex(dp) :: cab(kd+1,n), ca(n,n), cq(n,n), cv(n,4)
    integer :: st, found, i, j, k

    ab = 0.0_dp
    a = 0.0_dp
    do j=1,n
        ab(1,j) = 2.0_dp+cos(real(j,dp))
        a(j,j) = ab(1,j)
        do k=1,min(kd,n-j)
            ab(k+1,j) = 0.2_dp*cos(real(3*j+5*k,dp))
            a(j+k,j) = ab(k+1,j)
            a(j,j+k) = ab(k+1,j)
        end do
    end do
    call eig_banded(ab,w,st,vectors,lower=.true.,method='bulge',num_selected=found)
    call check(st == linalg_status_success .and. found == n,'real full eigenpairs')
    wref = w
    residual = maxval(abs(matmul(a,vectors)-matmul(vectors,diag(w))))
    call check(residual < 2.0e-12_dp,'real full residual')
    call eig_banded(ab,selected,st,sv,lower=.true.,method='bulge', &
        subset_by_index=[3,6],num_selected=found)
    call check(st == linalg_status_success .and. found == 4,'real selected eigenpairs')
    call check(maxval(abs(selected-w(3:6))) < 1.0e-12_dp,'selected eigenvalues')
    call check(maxval(abs(matmul(a,sv)-matmul(sv,diag(selected)))) < 1.0e-12_dp, &
        'selected real residual')
    call eigvals_banded(ab,vals,st,lower=.true.,method='bulge', &
        subset_by_index=[3,6],num_selected=found)
    call check(st == linalg_status_success .and. found == 4,'real selected values')
    call check(maxval(abs(selected-vals)) < 1.0e-13_dp,'selected eigvalsh results')
    call eig_banded(ab,selected(1:3),st,lower=.true.,method='bulge', &
        subset_by_index=[3,6],num_selected=found)
    call check(st == linalg_status_invalid_shape .and. found == 4, &
        'undersized selected values')
    call eig_banded(ab,selected,st,lower=.true.,method='bulge', &
        subset_by_value=[w(4),w(7)],num_selected=found)
    call check(st == linalg_status_success .and. found == 3,'half-open interval')
    call check(maxval(abs(selected(1:3)-w(5:7))) < 1.0e-12_dp,'interval values')
    call eig_banded(ab,w,st,lower=.true.,method='unrecognized')
    call check(st == linalg_status_invalid_shape,'reject unknown method')

    cab = cmplx(0.0_dp,0.0_dp,dp)
    ca = cmplx(0.0_dp,0.0_dp,dp)
    do j=1,n
        ca(j,j) = cmplx(1.0_dp+cos(real(j,dp)),0.0_dp,dp)
        cab(kd+1,j) = ca(j,j)
        do k=1,min(kd,j-1)
            i=j-k
            ca(i,j) = cmplx(0.15_dp*sin(real(5*i+j,dp)), &
                0.18_dp*cos(real(7*j+i,dp)),dp)
            ca(j,i) = conjg(ca(i,j))
            cab(kd+1-k,j) = ca(i,j)
        end do
    end do
    call eig_banded(cab,w,st,cq,method='bulge')
    call check(st == linalg_status_success,'complex full eigenpairs')
    residual = maxval(abs(matmul(ca,cq)-matmul(cq,diag_complex(w))))
    call check(residual < 3.0e-12_dp,'complex full residual')
    call eig_banded(cab,selected,st,cv,method='bulge', &
        subset_by_index=[2,5],num_selected=found)
    call check(st == linalg_status_success .and. found == 4,'complex selected eigenpairs')
    call check(maxval(abs(selected-w(2:5))) < 1.0e-12_dp,'complex selection')
    call check(maxval(abs(matmul(ca,cv)-matmul(cv,diag_complex(selected)))) < 3.0e-12_dp, &
        'complex selection residual')

    scl = 1.0e150_dp
    ab = ab*scl
    call eigvals_banded(ab,wscale,st,lower=.true.,method='bulge')
    call check(st == linalg_status_success,'large-scale eigenvalues')
    call check(maxval(abs(wscale/scl-wref)) < 2.0e-11_dp,'large-scale spectrum')
    ab = ab/scl ! Restore the original data before applying the small scale.
    ab = ab*1.0e-150_dp
    call eigvals_banded(ab,wscale,st,lower=.true.,method='bulge')
    call check(st == linalg_status_success,'small-scale eigenvalues')
    call check(maxval(abs(wscale/1.0e-150_dp-wref)) < 2.0e-11_dp, &
        'small-scale spectrum')
    print '(a)', 'banded bulge-chasing eigen native tests passed'

contains
    function diag(x) result(y)
        real(dp), intent(in) :: x(:) !! eigenvalues
        real(dp) :: y(size(x),size(x))
        integer :: l
        y = 0.0_dp
        do l=1,size(x)
            y(l,l) = x(l)
        end do
    end function diag
    function diag_complex(x) result(y)
        real(dp), intent(in) :: x(:) !! eigenvalues
        complex(dp) :: y(size(x),size(x))
        integer :: l
        y = cmplx(0.0_dp,0.0_dp,dp)
        do l=1,size(x)
            y(l,l) = cmplx(x(l),0.0_dp,dp)
        end do
    end function diag_complex
    subroutine check(condition,label)
        logical, intent(in) :: condition !! numerical assertion
        character(len=*), intent(in) :: label !! assertion name
        if (.not. condition) then
            print '(a)', 'FAILED: '//label
            error stop 1
        end if
    end subroutine check
end program test_linalg_banded_bulge
