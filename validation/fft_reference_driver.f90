! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
! Input-driven validation program; deliberately outside the unattended tests.
program fft_reference_driver
    use scifort_fft
    implicit none
    complex(dp), allocatable :: x(:), y(:)
    real(dp), allocatable :: r(:)
    integer :: op, real_input, nx, n, kind, mode, orth, code, ios, j
    real(dp) :: d, dreal, dimag
    character(8), parameter :: norms(0:2) = [character(8) :: 'backward','forward','ortho']
    do
        read(*,*,iostat=ios) op,real_input,nx,n,kind,mode,orth,d
        if (ios < 0) exit
        if (ios /= 0) error stop 'invalid driver header'
        allocate(x(nx))
        do j = 1, nx
            read(*,*) dreal, dimag
            x(j) = cmplx(dreal,dimag,dp)
        end do
        if (orth < 0) then
            call evaluate()
        else
            call evaluate(orth == 1)
        end if
        write(*,'(i0,1x,i0)') code,size(y)
        do j = 1, size(y)
            write(*,'(2(es26.17e3,1x))') real(y(j),dp),aimag(y(j))
        end do
        deallocate(x,y)
    end do
contains
    subroutine evaluate(orthogonal)
        implicit none
        logical, intent(in), optional :: orthogonal !! optional endpoint weighting
        code = 0
        select case (op)
        case (1)
            if (real_input == 1) then
                y = fft(real(x,dp),n,trim(norms(mode)),code)
            else
                y = fft(x,n,trim(norms(mode)),code)
            end if
        case (2)
            if (real_input == 1) then
                y = ifft(real(x,dp),n,trim(norms(mode)),code)
            else
                y = ifft(x,n,trim(norms(mode)),code)
            end if
        case (3)
            y = rfft(real(x,dp),n,trim(norms(mode)),code)
        case (4,5)
            if (op == 4) then
                if (real_input == 1) then
                    r = irfft(real(x,dp),n,trim(norms(mode)),code)
                else
                    r = irfft(x,n,trim(norms(mode)),code)
                end if
            else
                if (real_input == 1) then
                    r = hfft(real(x,dp),n,trim(norms(mode)),code)
                else
                    r = hfft(x,n,trim(norms(mode)),code)
                end if
            end if
            y = cmplx(r,0.0_dp,dp)
        case (6)
            y = ihfft(real(x,dp),n,trim(norms(mode)),code)

        case (7)
            if (real_input == 1) then
                r = dct(real(x,dp),kind,n,trim(norms(mode)),orthogonal,code)
                y = cmplx(r,0.0_dp,dp)
            else
                y = dct(x,kind,n,trim(norms(mode)),orthogonal,code)
            end if
        case (8)
            if (real_input == 1) then
                r = idct(real(x,dp),kind,n,trim(norms(mode)),orthogonal,code)
                y = cmplx(r,0.0_dp,dp)
            else
                y = idct(x,kind,n,trim(norms(mode)),orthogonal,code)
            end if
        case (9)
            if (real_input == 1) then
                r = dst(real(x,dp),kind,n,trim(norms(mode)),orthogonal,code)
                y = cmplx(r,0.0_dp,dp)
            else
                y = dst(x,kind,n,trim(norms(mode)),orthogonal,code)
            end if
        case (10)
            if (real_input == 1) then
                r = idst(real(x,dp),kind,n,trim(norms(mode)),orthogonal,code)
                y = cmplx(r,0.0_dp,dp)
            else
                y = idst(x,kind,n,trim(norms(mode)),orthogonal,code)
            end if
        case (11)
            r = fftfreq(n,d,code)
            y = cmplx(r,0.0_dp,dp)
        case (12)
            r = rfftfreq(n,d,code)
            y = cmplx(r,0.0_dp,dp)
        case (13)
            if (real_input == 1) then
                y = cmplx(fftshift(real(x,dp)),0.0_dp,dp)
            else
                y = fftshift(x)
            end if
        case (14)
            if (real_input == 1) then
                y = cmplx(ifftshift(real(x,dp)),0.0_dp,dp)
            else
                y = ifftshift(x)
            end if
        case default
            error stop 'unknown driver operation'
        end select
    end subroutine evaluate
end program fft_reference_driver
