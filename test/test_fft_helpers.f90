! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_fft_helpers
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_positive_inf, ieee_is_nan
    use scifort_fft
    implicit none
    real(dp) :: x(3,4), y(3,4), z(3,4), expected(3,4), gaps(6,4)
    real(dp) :: empty(0,3), empty_out(0,3), half_real(3,4)
    complex(dp) :: cx(3,4), cy(3,4), cz(3,4), half(3,3)
    complex(dp) :: batch(2,3,4), bout(2,3,4), bref(2,3,4), scalar
    real(dp) :: high(2,1,1,1,1,1,1,1,1,1,1,1,1,1,3)
    real(dp) :: hout(2,1,1,1,1,1,1,1,1,1,1,1,1,1,3)
    real(dp) :: hback(2,1,1,1,1,1,1,1,1,1,1,1,1,1,3)
    integer :: n, previous, following, a, b, code, i, j, checks, failures
    checks = 0
    failures = 0
    previous = 1
    do n = 1, 2048
        if (smooth(n)) previous = n
        following = n
        do while (.not. smooth(following))
            following = following+1
        end do
        a = next_fast_len(n,status=code)
        call check(code == 0 .and. a == following,'next length agrees with brute search')
        b = prev_fast_len(n,status=code)
        call check(code == 0 .and. b == previous,'previous length agrees with brute search')
        call check(next_fast_len(n,real=.true.) == a,'real and complex backend next lengths')
        call check(prev_fast_len(n,real=.false.) == b,'real and complex backend previous lengths')
    end do
    a = next_fast_len(0,status=code)
    call check(code == 1 .and. a == 0,'zero target rejected')
    a = prev_fast_len(-1,status=code)
    call check(code == 1 .and. a == 0,'negative target rejected')
    a = next_fast_len(huge(0),status=code)
    call check(code == 1 .and. a == 0,'next length integer overflow rejected')
    a = prev_fast_len(huge(0),status=code)
    call check(code == 0 .and. a > 0 .and. smooth(a),'previous maximum integer is smooth')
    do j = 1, 4
        do i = 1, 3
            x(i,j) = real(i+10*j,dp)
        end do
    end do
    cx = cmplx(x,-x,dp)
    expected = cshift(cshift(x,-1,dim=1),-2,dim=2)
    call fftshift_nd(x,y,status=code)
    call check(code == 0 .and. all(y == expected),'odd/even two-dimensional shift')
    call ifftshift_nd(y,z,status=code)
    call check(code == 0 .and. all(z == x),'real shift inverse')
    call fftshift_nd(cx,cy,axes=[2,1],status=code)
    call check(code == 0 .and. all(cy == cmplx(expected,-expected,dp)),'complex reordered axes')
    call ifftshift_nd(cy,cz,axes=[1,2])
    call check(all(cz == cx),'complex shift inverse')
    call fftshift_nd(x,y,axes=[2])
    call check(all(y == cshift(x,-2,dim=2)),'single-axis shift')
    call fftshift_nd(x,y,axes=[integer ::])
    call check(all(y == x),'empty axes identity')
    gaps(1:6:2,:) = x
    gaps(2:6:2,:) = -99
    call fftshift_nd(gaps(1:6:2,:),gaps(2:6:2,:),status=code)
    call check(code == 0 .and. all(gaps(2:6:2,:) == expected),'strided input and output')
    call check(all(gaps(1:6:2,:) == x),'strided input preserved')
    high = reshape([(real(i,dp),i=1,6)],shape(high))
    call fftshift_nd(high,hout,axes=[15,1],status=code)
    call check(code == 0,'rank fifteen shift')
    call ifftshift_nd(hout,hback,axes=[15,1])
    call check(all(hback == high),'rank fifteen inverse')
    call fftshift_nd(empty,empty_out,status=code)
    call check(code == 0,'empty array shift')
    call fftshift_nd(x,y,axes=[0],status=code)
    call check(code == 1 .and. all(ieee_is_nan(y)),'invalid axis')
    call fftshift_nd(x,y,axes=[1,1],status=code)
    call check(code == 1,'duplicate axes')
    call fftshift_nd(x,y,axes=[3],status=code)
    call check(code == 1,'out of range axis')
    call fftshift_nd(x,empty_out,status=code)
    call check(code == 1,'wrong shift output shape')
    x(1,1) = ieee_value(0.0_dp,ieee_quiet_nan)
    x(2,1) = ieee_value(0.0_dp,ieee_positive_inf)
    call fftshift_nd(x,y,status=code)
    call ifftshift_nd(y,z)
    call check(code == 0 .and. ieee_is_nan(z(1,1)),'NaN preserved by shifts')
    call check(z(2,1) == x(2,1),'infinity preserved by shifts')
    x = real(cx,dp)
    call fft2(cx,cy,status=code)
    call fftn(cx,cz)
    call check(code == 0 .and. maxval(abs(cy-cz)) == 0,'fft2 complex default axes')
    call fft2(x,cy)
    call fftn(x,cz)
    call check(maxval(abs(cy-cz)) == 0,'fft2 real')
    call ifft2(cx,cy)
    call ifftn(cx,cz)
    call check(maxval(abs(cy-cz)) == 0,'ifft2 complex')
    call ifft2(x,cy)
    call ifftn(x,cz)
    call check(maxval(abs(cy-cz)) == 0,'ifft2 real')
    call rfft2(x,half)
    call irfft2(half,y,s=[3,4])
    call check(maxval(abs(y-x)) < 1e-12_dp,'real 2D pair')
    call ihfft2(x,half)
    call hfft2(half,y,s=[3,4])
    call check(maxval(abs(y-x)) < 1e-12_dp,'Hermitian 2D pair')
    half_real = real(cx,dp)
    call irfft2(half_real,y,s=[3,4],status=code)
    call irfftn(half_real,z,s=[3,4])
    call check(code == 0 .and. maxval(abs(y-z)) == 0,'irfft2 real half input')
    call hfft2(half_real,y,s=[3,4],status=code)
    call hfftn(half_real,z,s=[3,4])
    call check(code == 0 .and. maxval(abs(y-z)) == 0,'hfft2 real half input')
    batch = reshape([(cmplx(real(i,dp),-real(i,dp),dp),i=1,24)],shape(batch))
    call fft2(batch,bout,status=code)
    call fftn(batch,bref,axes=[2,3])
    call check(code == 0 .and. maxval(abs(bout-bref)) == 0,'2D batch default last two axes')
    call fft2(batch,bout,axes=[3,1])
    call fftn(batch,bref,axes=[3,1])
    call check(maxval(abs(bout-bref)) == 0,'2D selected axes')
    call fft2(cx,cy,axes=[1],status=code)
    call check(code == 1,'2D requires two axes')
    call fft2(cx,cy,s=[3],status=code)
    call check(code == 1,'2D requires two lengths')
    call fft2(cx(:,1),cy(:,1),status=code)
    call check(code == 1,'2D requires rank at least two')
    scalar = (1.0_dp,0.0_dp)
    call fft2(scalar,cy,status=code)
    call check(code == 1,'scalar input rejected')
    if (failures > 0) error stop 'FFT helper tests failed'
    print '(a,i0,a)', 'FFT helpers: ',checks,' checks passed'
contains
    subroutine check(ok,label)
        implicit none
        logical, intent(in) :: ok !! condition required for success
        character(*), intent(in) :: label !! failure description
        checks = checks+1
        if (ok) return
        failures = failures+1
        print '(a)', 'FAIL: '//label
    end subroutine check
    pure logical function smooth(value) result(ok)
        implicit none
        integer, intent(in) :: value !! positive integer to factor
        integer :: remaining, i
        integer, parameter :: primes(4) = [2,3,5,7]
        remaining = value
        do i = 1, 4
            do while (mod(remaining,primes(i)) == 0)
                remaining = remaining/primes(i)
            end do
        end do
        ok = remaining == 1
    end function smooth
end program test_fft_helpers
