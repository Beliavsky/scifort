! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_fft
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_positive_inf, ieee_is_nan
    use scifort_fft
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    character(8), parameter :: norms(3) = [character(8) :: 'backward','forward','ortho']
    complex(dp), allocatable :: x(:), y(:), z(:), expected(:), gap(:), saved(:)
    real(dp), allocatable :: r(:), a(:), b(:), c(:), padded(:)
    real(dp) :: angle, weight, factor, logical_size, energy, nan, inf
    integer :: n, j, k, mode, kind, orth, family, code, checks, failures
    logical :: orthogonal
    checks = 0
    failures = 0
    nan = ieee_value(0.0_dp,ieee_quiet_nan)
    inf = ieee_value(0.0_dp,ieee_positive_inf)
    do n = 1, 33
        x = [(cmplx(sin(real(j,dp)),cos(0.7_dp*real(j,dp)),dp),j=1,n)]
        r = real(x,dp)
        expected = x
        do k = 0, n-1
            expected(k+1) = (0.0_dp,0.0_dp)
            do j = 0, n-1
                angle = -2.0_dp*pi*real(j*k,dp)/real(n,dp)
                expected(k+1) = expected(k+1)+x(j+1)*cmplx(cos(angle),sin(angle),dp)
            end do
        end do
        do mode = 1, 3
            factor = 1.0_dp
            if (mode == 2) factor = 1.0_dp/real(n,dp)
            if (mode == 3) factor = 1.0_dp/sqrt(real(n,dp))
            y = fft(x,norm=trim(norms(mode)),status=code)
            call check(code == 0,'complex FFT status')
            call close_complex(y,expected*factor,3.0e-13_dp,'independent DFT definition')
            z = ifft(y,norm=trim(norms(mode)),status=code)
            call check(code == 0,'complex inverse status')
            call close_complex(z,x,3.0e-13_dp,'complex FFT round trip')
            y = fft(r,norm=trim(norms(mode)))
            z = fft(cmplx(r,0.0_dp,dp),norm=trim(norms(mode)))
            call close_complex(y,z,0.0_dp,'real FFT overload')
            y = ifft(r,norm=trim(norms(mode)))
            z = ifft(cmplx(r,0.0_dp,dp),norm=trim(norms(mode)))
            call close_complex(y,z,0.0_dp,'real inverse overload')
            y = rfft(r,norm=trim(norms(mode)))
            call check(size(y) == n/2+1,'real half-spectrum size')
            a = irfft(y,n=n,norm=trim(norms(mode)))
            call close_real(a,r,3.0e-13_dp,'real FFT round trip including odd lengths')
            y = ihfft(r,norm=trim(norms(mode)))
            a = hfft(y,n=n,norm=trim(norms(mode)))
            call close_real(a,r,3.0e-13_dp,'Hermitian FFT round trip')
        end do
        gap = [(cmplx(-111.0_dp,7.0_dp,dp),j=1,2*n)]
        gap(1:2*n:2) = x
        saved = gap
        y = fft(gap(1:2*n:2))
        call close_complex(y,expected,3.0e-13_dp,'strided input')
        call close_complex(gap,saved,0.0_dp,'input and gaps preserved')
    end do
    do n = 1, 12
        r = [(sin(0.73_dp*real(j,dp)),j=1,n)]
        do family = 1, 2
            do kind = 1, 4
                if (family == 1 .and. kind == 1 .and. n == 1) cycle
                do mode = 1, 3
                    do orth = 0, 1
                        orthogonal = orth == 1
                        logical_size = 2.0_dp*real(n,dp)
                        if (kind == 1) then
                            if (family == 1) logical_size = 2.0_dp*real(n-1,dp)
                            if (family == 2) logical_size = 2.0_dp*real(n+1,dp)
                        end if
                        factor = 1.0_dp
                        if (mode == 2) factor = 1.0_dp/logical_size
                        if (mode == 3) factor = 1.0_dp/sqrt(logical_size)
                        b = r*0.0_dp
                        ! Direct defining matrices are independent of the FFT extension.
                        do k = 0, n-1
                            do j = 0, n-1
                                weight = 2.0_dp
                                select case (kind)
                                case (1)
                                    if (family == 1) then
                                        angle = pi*real(j*k,dp)/real(n-1,dp)
                                        if (j == 0 .or. j == n-1) weight = 1.0_dp
                                    else
                                        angle = pi*real((j+1)*(k+1),dp)/real(n+1,dp)
                                    end if
                                case (2)
                                    angle = pi*real(2*j+1,dp)*real(k,dp)/real(2*n,dp)
                                    if (family == 2) angle = &
                                        pi*real(2*j+1,dp)*real(k+1,dp)/real(2*n,dp)
                                case (3)
                                    angle = pi*real(2*k+1,dp)*real(j,dp)/real(2*n,dp)
                                    if (family == 1 .and. j == 0) weight = 1.0_dp
                                    if (family == 2) then
                                        angle = pi*real(2*k+1,dp)*real(j+1,dp)/real(2*n,dp)
                                        if (j == n-1) weight = 1.0_dp
                                    end if
                                case (4)
                                    angle = pi*real(2*j+1,dp)*real(2*k+1,dp)/real(4*n,dp)
                                end select
                                if (orthogonal) then
                                    if (family == 1) then
                                        if (kind == 1 .and. (j == 0 .or. j == n-1)) weight = weight*sqrt(2.0_dp)
                                        if (kind == 3 .and. j == 0) weight = weight*sqrt(2.0_dp)
                                        if (kind == 1 .and. (k == 0 .or. k == n-1)) weight = weight/sqrt(2.0_dp)
                                        if (kind == 2 .and. k == 0) weight = weight/sqrt(2.0_dp)
                                    else
                                        if (kind == 3 .and. j == n-1) weight = weight*sqrt(2.0_dp)
                                        if (kind == 2 .and. k == n-1) weight = weight/sqrt(2.0_dp)
                                    end if
                                end if
                                if (family == 1) then
                                    b(k+1) = b(k+1)+factor*weight*r(j+1)*cos(angle)
                                else
                                    b(k+1) = b(k+1)+factor*weight*r(j+1)*sin(angle)
                                end if
                            end do
                        end do
                        if (family == 1) then
                            a = dct(r,type=kind,norm=trim(norms(mode)),orthogonalize=orthogonal,status=code)
                            c = idct(a,type=kind,norm=trim(norms(mode)),orthogonalize=orthogonal)
                            y = dct(cmplx(r,r*0.3_dp,dp),type=kind,norm=trim(norms(mode)), &
                                orthogonalize=orthogonal)
                        else
                            a = dst(r,type=kind,norm=trim(norms(mode)),orthogonalize=orthogonal,status=code)
                            c = idst(a,type=kind,norm=trim(norms(mode)),orthogonalize=orthogonal)
                            y = dst(cmplx(r,r*0.3_dp,dp),type=kind,norm=trim(norms(mode)), &
                                orthogonalize=orthogonal)
                        end if
                        call check(code == 0,'real transform status')
                        call close_real(a,b,2.0e-13_dp,'independent cosine/sine matrix')
                        call close_real(c,r,2.0e-13_dp,'cosine/sine inverse round trip')
                        call close_complex(y,cmplx(a,a*0.3_dp,dp),2.0e-13_dp,'complex linearity')
                        if (mode == 3 .and. orthogonal) then
                            energy = sum(r*r)
                            call check(abs(sum(a*a)-energy) <= 2.0e-13_dp*max(1.0_dp,energy), &
                                'orthonormal energy preservation')
                        end if
                    end do
                end do
            end do
        end do
    end do
    r = [1.0_dp,2.0_dp,3.0_dp]
    y = fft(r,n=5)
    padded = [r,0.0_dp,0.0_dp]
    z = fft(padded)
    call close_complex(y,z,0.0_dp,'FFT zero padding')
    y = fft(r,n=2)
    z = fft(r(:2))
    call close_complex(y,z,0.0_dp,'FFT cropping')
    a = dct(r,n=5)
    b = dct(padded)
    call close_real(a,b,0.0_dp,'DCT zero padding')
    a = dst(r,n=2)
    b = dst(r(:2))
    call close_real(a,b,0.0_dp,'DST cropping')
    a = fftfreq(4)
    call close_real(a,[0.0_dp,0.25_dp,-0.5_dp,-0.25_dp],0.0_dp,'even signed frequency bins')
    a = rfftfreq(5,d=-0.25_dp)
    call close_real(a,[0.0_dp,-0.8_dp,-1.6_dp],1.0e-15_dp,'odd real bins, negative spacing')
    a = fftfreq(4,d=huge(1.0_dp),status=code)
    call check(code == 0 .and. a(2) > 0.0_dp,'large spacing retains subnormal frequency')
    a = fftfreq(1,d=nearest(0.0_dp,1.0_dp),status=code)
    call check(code == 0 .and. a(1) == 0.0_dp,'singleton bins with smallest spacing')
    r = [1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp]
    a = fftshift(r)
    call close_real(a,[4.0_dp,5.0_dp,1.0_dp,2.0_dp,3.0_dp],0.0_dp,'odd FFT shift')
    a = ifftshift(r)
    call close_real(a,[3.0_dp,4.0_dp,5.0_dp,1.0_dp,2.0_dp],0.0_dp,'odd inverse shift')
    x = cmplx(r,r,dp)
    y = ifftshift(fftshift(x))
    call close_complex(y,x,0.0_dp,'complex shift round trip')
    a = fftshift([nan,1.0_dp,inf])
    call check(ieee_is_nan(a(2)) .and. a(1) == inf,'shifts preserve nonfinite values')
    x = [cmplx(1.0_dp,99.0_dp,dp),cmplx(2.0_dp,3.0_dp,dp),cmplx(4.0_dp,77.0_dp,dp)]
    a = irfft(x,n=4)
    x(1) = cmplx(1.0_dp,0.0_dp,dp)
    x(3) = cmplx(4.0_dp,0.0_dp,dp)
    b = irfft(x,n=4)
    call close_real(a,b,0.0_dp,'self-conjugate endpoint imaginary parts ignored')
    y = fft([real(dp) ::],n=3,status=code)
    call check(code == 0 .and. size(y) == 3 .and. all(abs(y) == 0.0_dp),'empty FFT with explicit padding')
    y = rfft([real(dp) ::],status=code)
    call check(code == 1 .and. size(y) == 0,'empty transform without length rejected')
    a = hfft([cmplx(2.0_dp,7.0_dp,dp)],status=code)
    call check(code == 1 .and. size(a) == 0,'Hermitian default requires a positive inferred length')
    a = hfft([cmplx(2.0_dp,7.0_dp,dp)],n=1,status=code)
    call check(code == 0 .and. size(a) == 1 .and. a(1) == 2.0_dp,'explicit singleton Hermitian transform')
    y = fft([1.0_dp],n=0,status=code)
    call check(code == 1 .and. size(y) == 0,'zero FFT length rejected')
    y = ifft([1.0_dp],norm='invalid',status=code)
    call check(code == 1 .and. size(y) == 0,'invalid normalization rejected')
    y = fft([nan],status=code)
    call check(code == 1 .and. size(y) == 0,'NaN real input rejected without traps')
    y = ifft([cmplx(0.0_dp,inf,dp)],status=code)
    call check(code == 1 .and. size(y) == 0,'infinite imaginary input rejected')
    a = dct([1.0_dp],type=1,status=code)
    call check(code == 1 .and. size(a) == 0,'singleton DCT I rejected')
    a = dst([3.0_dp],type=1,status=code)
    call check(code == 0 .and. size(a) == 1 .and. abs(a(1)-6.0_dp) < 1.0e-14_dp,'singleton DST I')
    a = dct([1.0_dp],type=0,status=code)
    call check(code == 1 .and. size(a) == 0,'invalid cosine type rejected')
    y = idst([cmplx(nan,0.0_dp,dp)],status=code)
    call check(code == 1 .and. size(y) == 0,'invalid complex sine input rejected')
    a = dst([1.0_dp],n=huge(1),status=code)
    call check(code == 1 .and. size(a) == 0,'FFT extension integer overflow rejected before allocation')
    a = fftfreq(0,status=code)
    call check(code == 1 .and. size(a) == 0,'invalid bin count')
    a = fftfreq(3,d=0.0_dp,status=code)
    call check(code == 1 .and. size(a) == 0,'zero sample spacing')
    a = fftfreq(3,d=nan,status=code)
    call check(code == 1 .and. size(a) == 0,'NaN sample spacing')
    a = rfftfreq(3,d=inf,status=code)
    call check(code == 1 .and. size(a) == 0,'infinite sample spacing')
    a = rfftfreq(4,d=nearest(0.0_dp,1.0_dp),status=code)
    call check(code == 1 .and. size(a) == 0,'unrepresentable bins rejected without overflow')
    if (failures /= 0) error stop 1
    print '(a,i0,a)', 'Standalone FFT: ',checks,' checks passed'
contains
    subroutine check(ok,label)
        implicit none
        logical, intent(in) :: ok !! expected invariant
        character(*), intent(in) :: label !! failure description
        checks = checks+1
        if (ok) return
        failures = failures+1
        print '(2a)', 'FAIL: ',label
    end subroutine check
    subroutine close_real(actual,reference,tolerance,label)
        implicit none
        real(dp), intent(in) :: actual(:) !! computed vector
        real(dp), intent(in) :: reference(:) !! independent expected vector
        real(dp), intent(in) :: tolerance !! scaled maximum-error bound
        character(*), intent(in) :: label !! failure description
        call check(size(actual) == size(reference),label//' shape')
        if (size(actual) /= size(reference)) return
        call check(all(abs(actual-reference) <= tolerance*max(1.0_dp,maxval(abs(reference)))),label)
    end subroutine close_real
    subroutine close_complex(actual,reference,tolerance,label)
        implicit none
        complex(dp), intent(in) :: actual(:) !! computed vector
        complex(dp), intent(in) :: reference(:) !! independent expected vector
        real(dp), intent(in) :: tolerance !! scaled maximum-error bound
        character(*), intent(in) :: label !! failure description
        call check(size(actual) == size(reference),label//' shape')
        if (size(actual) /= size(reference)) return
        call check(all(abs(actual-reference) <= tolerance*max(1.0_dp,maxval(abs(reference)))),label)
    end subroutine close_complex
end program test_fft
