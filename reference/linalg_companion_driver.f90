! SPDX-License-Identifier: MIT
! Simple text protocol used by tools/compare_linalg_companion.py.
! Input:  number of cases, then for each: mode type n, coefficient vector.
! Modes C/P/F are companion/Fiedler companion/Fiedler; types R/Z/I.
! Real and integer cases use one coefficient line; complex use separate real
! and imaginary lines. Output: status and each matrix row (real/imag separate).
program linalg_companion_driver
    use scifort_kinds, only : dp
    use scifort_linalg_companion, only : companion, fiedler_companion
    use scifort_linalg_fiedler, only : fiedler
    implicit none
    integer :: cases, t, n, d, i, s
    character(1) :: mode, dtype
    real(dp), allocatable :: a(:), b(:), m(:,:)
    complex(dp), allocatable :: z(:), zm(:,:)
    integer, allocatable :: ia(:), im(:,:)

    read(*,*) cases
    do t=1,cases
        read(*,*) mode,dtype,n
        d=n
        if (mode /= 'F') d=max(0,n-1)
        allocate(a(n), b(n), z(n), ia(n), m(d,d), zm(d,d), im(d,d))
        select case(dtype)
        case('R')
            if (n > 0) read(*,*) a
            select case(mode)
            case('C')
                call companion(a,m,s)
            case('P')
                call fiedler_companion(a,m,s)
            case('F')
                call fiedler(a,m,s)
            end select
        case('Z')
            if (n > 0) then
                read(*,*) a
                read(*,*) b
            end if
            z=cmplx(a,b,kind=dp)
            select case(mode)
            case('C')
                call companion(z,zm,s)
            case('P')
                call fiedler_companion(z,zm,s)
            case('F')
                call fiedler(z,m,s)
            end select
        case('I')
            if (n > 0) read(*,*) ia
            select case(mode)
            case('C')
                call companion(ia,m,s)
            case('P')
                call fiedler_companion(ia,m,s)
            case('F')
                call fiedler(ia,im,s)
                m=real(im,dp)
            end select
        end select
        write(*,'(I0)') s
        do i=1,d
            if (dtype == 'Z' .and. mode /= 'F') then
                write(*,'(*(ES26.17E3,1x))') real(zm(i,:),dp)
                write(*,'(*(ES26.17E3,1x))') aimag(zm(i,:))
            else
                write(*,'(*(ES26.17E3,1x))') m(i,:)
            end if
        end do
        deallocate(a,b,z,ia,m,zm,im)
    end do
end program linalg_companion_driver
