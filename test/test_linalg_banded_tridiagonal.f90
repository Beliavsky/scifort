! SPDX-License-Identifier: MIT
program test_linalg_banded_tridiagonal
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    use scifort_linalg_banded_tridiagonal, only : band_to_tridiagonal
    implicit none
    integer :: n, kd, mode, st, i, j, k
    real(dp), allocatable :: a(:,:), ab(:,:), q(:,:), d(:), e(:), t(:,:), res(:,:)
    complex(dp), allocatable :: ca(:,:), cab(:,:), cq(:,:), ct(:,:), cres(:,:)
    real(dp) :: scale, error

    do mode = 1,2
        do n = 0,65,13
            do kd = 0,9,3
                allocate(a(n,n),ab(kd+1,n),q(n,n),d(n),e(max(0,n-1)),t(n,n),res(n,n))
                allocate(ca(n,n),cab(kd+1,n),cq(n,n),ct(n,n),cres(n,n))
                a = 0.0_dp
                ca = cmplx(0.0_dp,0.0_dp,dp)
                ab = 0.0_dp
                cab = cmplx(0.0_dp,0.0_dp,dp)
                do j = 1,n
                    a(j,j) = 2.0_dp+sin(real(5*j,dp))
                    ca(j,j) = cmplx(a(j,j),0.0_dp,dp)
                    do k = 1,min(kd,n-j)
                        i = j+k
                        a(i,j) = 0.37_dp*cos(real(7*i+13*j,dp))
                        a(j,i) = a(i,j)
                        ca(i,j) = cmplx(a(i,j),0.31_dp*sin(real(3*i+11*j,dp)),dp)
                        ca(j,i) = conjg(ca(i,j))
                    end do
                end do
                do j = 1,n
                    do k = 0,min(kd,n-j)
                        if (mode == 1) then
                            ab(k+1,j) = a(j+k,j)
                            cab(k+1,j) = ca(j+k,j)
                        else
                            ab(kd+1-k,j+k) = a(j,j+k)
                            cab(kd+1-k,j+k) = ca(j,j+k)
                        end if
                    end do
                end do
                call band_to_tridiagonal(ab,d,e,st,q,lower=mode==1)
                call require(st == linalg_status_success,'real band-to-tridiagonal status')
                call require(all(e >= 0.0_dp),'real positive offdiagonals')
                t = 0.0_dp
                do j = 1,n
                    t(j,j) = d(j)
                end do
                do j = 1,n-1
                    t(j,j+1) = e(j)
                    t(j+1,j) = e(j)
                end do
                scale = max(1.0_dp,maxval(abs(a)))
                res = matmul(q,t)
                res = matmul(res,transpose(q))
                error = maxval(abs(res-a))/scale
                call require(error < 1.0e-11_dp,'real similarity reconstruction')
                call band_to_tridiagonal(cab,d,e,st,cq,lower=mode==1)
                call require(st == linalg_status_success,'complex band-to-tridiagonal status')
                call require(all(e >= 0.0_dp),'complex positive offdiagonals')
                ct = cmplx(0.0_dp,0.0_dp,dp)
                do j = 1,n
                    ct(j,j) = cmplx(d(j),0.0_dp,dp)
                end do
                do j = 1,n-1
                    ct(j,j+1) = cmplx(e(j),0.0_dp,dp)
                    ct(j+1,j) = cmplx(e(j),0.0_dp,dp)
                end do
                scale = max(1.0_dp,maxval(abs(ca)))
                cres = matmul(cq,ct)
                cres = matmul(cres,conjg(transpose(cq)))
                error = maxval(abs(cres-ca))/scale
                call require(error < 1.0e-11_dp,'complex similarity reconstruction')
                deallocate(a,ab,q,d,e,t,res,ca,cab,cq,ct,cres)
            end do
        end do
    end do
    allocate(ab(2,3),d(3),e(2),q(3,3))
    ab = 0.0_dp
    ab(2,:) = 1.0_dp
    ab(1,1) = ieee_value(0.0_dp,ieee_quiet_nan) ! unused upper entry
    call band_to_tridiagonal(ab,d,e,st,q)
    call require(st == linalg_status_success,'unused NaN accepted')
    ab(1,2) = ieee_value(0.0_dp,ieee_quiet_nan)
    call band_to_tridiagonal(ab,d,e,st,q)
    call require(st == linalg_status_no_convergence,'valid NaN rejected')
    call band_to_tridiagonal(ab,d,e(1:1),st,q)
    call require(st == linalg_status_invalid_shape,'invalid shape')
    print '(a)', 'banded bulge-chasing tridiagonal native tests passed'

contains
    subroutine require(condition,label)
        logical, intent(in) :: condition !! required invariant
        character(len=*), intent(in) :: label !! failure diagnostic
        if (.not. condition) then
            print '(a)', 'FAILED: '//label
            error stop 1
        end if
    end subroutine require
end program test_linalg_banded_tridiagonal
