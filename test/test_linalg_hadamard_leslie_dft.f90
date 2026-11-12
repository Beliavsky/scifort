! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_linalg_hadamard_leslie_dft
    use, intrinsic :: iso_fortran_env, only : int32, int64
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape
    use scifort_linalg_hadamard, only : hadamard
    use scifort_linalg_leslie, only : leslie
    use scifort_linalg_dft, only : dft
    implicit none

    integer :: i, j, n, power, status
    integer(int32) :: h_bad(3,3), h32(4,4)
    integer(int64) :: h64(4,4), f64(3), s64(2), l64(3,3)
    real(dp) :: hreal(4,4), f(4), s(3), l(4,4)
    real(dp) :: fb(3,4), sb(1,3), lb(3,4,4), lb_bad(2,4,4)
    complex(dp) :: hc(4,4), z(4), zs(3), lc(4,4), lmixed(4,4)
    complex(dp) :: zbatch(1,4), lbatchc(3,4,4), z_survival(2,3), zs_batch(3,3)
    complex(dp) :: bad_dft(3,3), empty(0,0), one(1,1)
    complex(dp), allocatable :: fourier(:,:), gram(:,:), v(:), transformed(:)
    integer(int64), allocatable :: h(:,:), gram_integer(:,:)
    real(dp), parameter :: tol = 1.0e-12_dp

    call hadamard(4,h32,status)
    call check(status == linalg_status_success,'Hadamard int32 status')
    call check(all(h32(1,:) == [1,1,1,1]),'Hadamard first row')
    call check(all(h32(2,:) == [1,-1,1,-1]),'Hadamard second row')
    call check(all(h32(4,:) == [1,-1,-1,1]),'Hadamard last row')
    call hadamard(4,h64,status)
    call check(status == 0 .and. all(h64 == int(h32,int64)),'Hadamard int64')
    call hadamard(4,hreal,status)
    call check(status == 0 .and. all(hreal == real(h32,dp)),'Hadamard real64')
    call hadamard(4,hc,status)
    call check(status == 0 .and. all(hc == cmplx(hreal,0.0_dp,dp)),'Hadamard complex64')
    do power=0,8
        n=2**power
        allocate(h(n,n),gram_integer(n,n))
        call hadamard(n,h,status)
        call check(status == 0,'Hadamard power-of-two status')
        call check(all(abs(h) == 1_int64),'Hadamard entries +-1')
        gram_integer = matmul(h,transpose(h))
        do i=1,n
            gram_integer(i,i) = gram_integer(i,i)-int(n,int64)
        end do
        call check(all(gram_integer == 0_int64),'Hadamard orthogonality')
        deallocate(h,gram_integer)
    end do
    h_bad = 5
    call hadamard(3,h_bad,status)
    call check(status == linalg_status_invalid_shape .and. all(h_bad == 0),'invalid order')
    h_bad = 5
    call hadamard(4,h_bad,status)
    call check(status == linalg_status_invalid_shape .and. all(h_bad == 0),'invalid Hadamard shape')
    call hadamard(0,h_bad,status)
    call check(status == linalg_status_invalid_shape,'zero Hadamard order')

    f = [0.1_dp,2.0_dp,1.0_dp,0.1_dp]
    s = [0.2_dp,0.8_dp,0.7_dp]
    call leslie(f,s,l,status)
    call check(status == linalg_status_success,'Leslie example status')
    call check(all(l(1,:) == f),'Leslie first row')
    do i=2,4
        call check(l(i,i-1) == s(i-1),'Leslie subdiagonal')
        do j=1,4
            if (j /= i-1) call check(l(i,j) == 0.0_dp,'Leslie zero entries')
        end do
    end do
    f64 = [2_int64,4_int64,3_int64]
    s64 = [5_int64,7_int64]
    call leslie(f64,s64,l64,status)
    call check(status == 0 .and. all(l64(1,:) == f64) .and. l64(3,2) == 7_int64, &
        'integer Leslie example')
    z = cmplx(f,[1.0_dp,-0.5_dp,0.0_dp,2.0_dp],dp)
    zs = cmplx(s,[2.0_dp,-1.0_dp,0.25_dp],dp)
    call leslie(z,zs,lc,status)
    call check(status == 0 .and. all(lc(1,:) == z) .and. lc(3,2) == zs(2), &
        'complex Leslie')
    call leslie(f,zs,lmixed,status)
    call check(status == 0 .and. all(lmixed(1,:) == cmplx(f,0.0_dp,dp)) .and. &
        lmixed(3,2) == zs(2),'real/complex Leslie')
    call leslie(z,s,lmixed,status)
    call check(status == 0 .and. all(lmixed(1,:) == z) .and. &
        lmixed(3,2) == cmplx(s(2),0.0_dp,dp),'complex/real Leslie')
    do i=1,3
        fb(i,:) = f*real(i,dp)
    end do
    sb(1,:) = s
    call leslie(fb,sb,lb,status)
    call check(status == 0 .and. all(lb(3,1,:) == fb(3,:)), &
        'broadcast Leslie fecundities')
    call check(all(lb(:,3,2) == s(2)),'broadcast Leslie survival')
    zbatch(1,:) = z
    zs_batch(1,:) = zs
    zs_batch(2,:) = zs*2.0_dp
    zs_batch(3,:) = zs*3.0_dp
    call leslie(zbatch,zs_batch,lbatchc,status)
    call check(status == 0 .and. all(lbatchc(3,1,:) == z) .and. &
        lbatchc(3,4,3) == zs_batch(3,3),'complex Leslie batch broadcasting')
    call leslie(fb,zs_batch,lbatchc,status)
    call check(status == 0 .and. all(lbatchc(2,1,:) == cmplx(fb(2,:),0.0_dp,dp)) .and. &
        lbatchc(2,3,2) == zs_batch(2,2),'mixed Leslie batch real/complex')
    call leslie(zbatch,sb,lbatchc,status)
    call check(status == linalg_status_invalid_shape, 'wrong Leslie batch output rejected')
    z_survival(1,:) = zs
    z_survival(2,:) = zs
    call leslie(fb,z_survival,lbatchc,status)
    call check(status == linalg_status_invalid_shape,'incompatible batch sizes')
    lb_bad = -1.0_dp
    call leslie(fb,sb,lb_bad,status)
    call check(status == linalg_status_invalid_shape .and. all(lb_bad == 0.0_dp), &
        'bad Leslie batch output shape')
    l = 4.0_dp
    call leslie(f,s(:2),l,status)
    call check(status == linalg_status_invalid_shape .and. all(l == 0.0_dp), &
        'bad Leslie vector lengths')

    call dft(0,empty,status)
    call check(status == 0,'empty DFT')
    call dft(1,one,status,scale='sqrtn')
    call check(status == 0 .and. one(1,1) == cmplx(1.0_dp,0.0_dp,dp),'DFT order 1')
    do n=2,65
        allocate(fourier(n,n),gram(n,n),v(n),transformed(n))
        call dft(n,fourier,status,scale='sqrtn')
        call check(status == 0,'unitary DFT status')
        gram = matmul(fourier,conjg(transpose(fourier)))
        do i=1,n
            gram(i,i) = gram(i,i)-cmplx(1.0_dp,0.0_dp,dp)
        end do
        call check(maxval(abs(gram)) < tol,'unitary DFT orthogonality')
        do i=1,n
            v(i) = cmplx(cos(real(i,dp)),sin(real(i*i,dp)),dp)
        end do
        call dft(n,fourier,status)
        transformed = matmul(fourier,v)
        call check(abs(transformed(1)-sum(v)) < 3.0e-13_dp,'DFT zero frequency')
        call dft(n,fourier,status,scale='n')
        call check(maxval(abs(matmul(fourier,v)*real(n,dp)-transformed)) < 2.0e-11_dp, &
            'DFT normalized mode')
        deallocate(fourier,gram,v,transformed)
    end do
    bad_dft = 3.0_dp
    call dft(4,bad_dft,status)
    call check(status == linalg_status_invalid_shape .and. all(bad_dft == cmplx(0.0_dp,0.0_dp,dp)), &
        'DFT wrong shape')
    bad_dft = 3.0_dp
    call dft(3,bad_dft,status,scale='sqrt')
    call check(status == linalg_status_invalid_shape .and. all(bad_dft == cmplx(0.0_dp,0.0_dp,dp)), &
        'DFT wrong scale')
    print '(a)', 'test_linalg_hadamard_leslie_dft: PASS'

contains
    subroutine check(condition,label)
        logical, intent(in) :: condition !! required assertion
        character(*), intent(in) :: label !! failure diagnostic
        if (.not. condition) then
            print '(a,1x,a)', 'FAILED:',trim(label)
            error stop 1
        end if
    end subroutine check
end program test_linalg_hadamard_leslie_dft
