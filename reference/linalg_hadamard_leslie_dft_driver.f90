! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
! Streaming driver for tools/compare_linalg_hadamard_leslie_dft.py.
! Input: case count; then mode, dtype, n, bf, bs, scale and any Leslie data.
! Each success emits status and one (real,imag) pair per matrix entry.
program linalg_hadamard_leslie_dft_driver
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_kinds, only : dp
    use scifort_linalg_hadamard, only : hadamard
    use scifort_linalg_leslie, only : leslie
    use scifort_linalg_dft, only : dft
    implicit none
    integer :: count, case_id, n, bf, bs, b, i, j, k, status, nalloc
    character(1) :: mode, dtype
    character(10) :: scale
    real(dp), allocatable :: hr(:,:)
    integer(int64), allocatable :: hi(:,:)
    complex(dp), allocatable :: hc(:,:), fourier(:,:)

    read(*,*) count
    do case_id=1,count
        read(*,*) mode,dtype,n,bf,bs,scale
        nalloc = max(0,n)
        select case (mode)
        case ('H')
            select case (dtype)
            case ('I')
                allocate(hi(nalloc,nalloc))
                call hadamard(n,hi,status)
                write(*,'(i0)') status
                if (status == 0) then
                    do i=1,n
                        do j=1,n
                            call emit(cmplx(real(hi(i,j),dp),0.0_dp,dp))
                        end do
                    end do
                end if
                deallocate(hi)
            case ('R')
                allocate(hr(nalloc,nalloc))
                call hadamard(n,hr,status)
                write(*,'(i0)') status
                if (status == 0) then
                    do i=1,n
                        do j=1,n
                            call emit(cmplx(hr(i,j),0.0_dp,dp))
                        end do
                    end do
                end if
                deallocate(hr)
            case ('C')
                allocate(hc(nalloc,nalloc))
                call hadamard(n,hc,status)
                write(*,'(i0)') status
                if (status == 0) then
                    do i=1,n
                        do j=1,n
                            call emit(hc(i,j))
                        end do
                    end do
                end if
                deallocate(hc)
            case default
                error stop 'unknown Hadamard dtype'
            end select
        case ('D')
            allocate(fourier(nalloc,nalloc))
            call dft(n,fourier,status,scale=trim(scale))
            write(*,'(i0)') status
            if (status == 0) then
                do i=1,n
                    do j=1,n
                        call emit(fourier(i,j))
                    end do
                end do
            end if
            deallocate(fourier)
        case ('L','B')
            b = max(bf,bs)
            select case (dtype)
            case ('R')
                if (mode == 'L') then
                    block
                        real(dp), allocatable :: f(:)
                        real(dp), allocatable :: survival(:)
                        real(dp), allocatable :: a(:,:)
                        allocate(f(nalloc),survival(max(0,n-1)),a(nalloc,nalloc))
                        read(*,*) f
                        read(*,*) survival
                        call leslie(f,survival,a,status)
                        write(*,'(i0)') status
                        if (status == 0) then
                            do i=1,n
                                do j=1,n
                                    call emit(cmplx(a(i,j),kind=dp))
                                end do
                            end do
                        end if
                    end block
                else
                    block
                        real(dp), allocatable :: f(:,:)
                        real(dp), allocatable :: survival(:,:)
                        real(dp), allocatable :: a(:,:,:)
                        allocate(f(bf,nalloc),survival(bs,max(0,n-1)),a(b,nalloc,nalloc))
                        do k=1,bf
                            read(*,*) f(k,:)
                        end do
                        do k=1,bs
                            read(*,*) survival(k,:)
                        end do
                        call leslie(f,survival,a,status)
                        write(*,'(i0)') status
                        if (status == 0) then
                            do k=1,b
                                do i=1,n
                                    do j=1,n
                                        call emit(cmplx(a(k,i,j),kind=dp))
                                    end do
                                end do
                            end do
                        end if
                    end block
                end if
            case ('C')
                if (mode == 'L') then
                    block
                        complex(dp), allocatable :: f(:)
                        complex(dp), allocatable :: survival(:)
                        complex(dp), allocatable :: a(:,:)
                        allocate(f(nalloc),survival(max(0,n-1)),a(nalloc,nalloc))
                        read(*,*) f
                        read(*,*) survival
                        call leslie(f,survival,a,status)
                        write(*,'(i0)') status
                        if (status == 0) then
                            do i=1,n
                                do j=1,n
                                    call emit(cmplx(a(i,j),kind=dp))
                                end do
                            end do
                        end if
                    end block
                else
                    block
                        complex(dp), allocatable :: f(:,:)
                        complex(dp), allocatable :: survival(:,:)
                        complex(dp), allocatable :: a(:,:,:)
                        allocate(f(bf,nalloc),survival(bs,max(0,n-1)),a(b,nalloc,nalloc))
                        do k=1,bf
                            read(*,*) f(k,:)
                        end do
                        do k=1,bs
                            read(*,*) survival(k,:)
                        end do
                        call leslie(f,survival,a,status)
                        write(*,'(i0)') status
                        if (status == 0) then
                            do k=1,b
                                do i=1,n
                                    do j=1,n
                                        call emit(cmplx(a(k,i,j),kind=dp))
                                    end do
                                end do
                            end do
                        end if
                    end block
                end if
            case ('X')
                if (mode == 'L') then
                    block
                        real(dp), allocatable :: f(:)
                        complex(dp), allocatable :: survival(:)
                        complex(dp), allocatable :: a(:,:)
                        allocate(f(nalloc),survival(max(0,n-1)),a(nalloc,nalloc))
                        read(*,*) f
                        read(*,*) survival
                        call leslie(f,survival,a,status)
                        write(*,'(i0)') status
                        if (status == 0) then
                            do i=1,n
                                do j=1,n
                                    call emit(cmplx(a(i,j),kind=dp))
                                end do
                            end do
                        end if
                    end block
                else
                    block
                        real(dp), allocatable :: f(:,:)
                        complex(dp), allocatable :: survival(:,:)
                        complex(dp), allocatable :: a(:,:,:)
                        allocate(f(bf,nalloc),survival(bs,max(0,n-1)),a(b,nalloc,nalloc))
                        do k=1,bf
                            read(*,*) f(k,:)
                        end do
                        do k=1,bs
                            read(*,*) survival(k,:)
                        end do
                        call leslie(f,survival,a,status)
                        write(*,'(i0)') status
                        if (status == 0) then
                            do k=1,b
                                do i=1,n
                                    do j=1,n
                                        call emit(cmplx(a(k,i,j),kind=dp))
                                    end do
                                end do
                            end do
                        end if
                    end block
                end if
            case ('Y')
                if (mode == 'L') then
                    block
                        complex(dp), allocatable :: f(:)
                        real(dp), allocatable :: survival(:)
                        complex(dp), allocatable :: a(:,:)
                        allocate(f(nalloc),survival(max(0,n-1)),a(nalloc,nalloc))
                        read(*,*) f
                        read(*,*) survival
                        call leslie(f,survival,a,status)
                        write(*,'(i0)') status
                        if (status == 0) then
                            do i=1,n
                                do j=1,n
                                    call emit(cmplx(a(i,j),kind=dp))
                                end do
                            end do
                        end if
                    end block
                else
                    block
                        complex(dp), allocatable :: f(:,:)
                        real(dp), allocatable :: survival(:,:)
                        complex(dp), allocatable :: a(:,:,:)
                        allocate(f(bf,nalloc),survival(bs,max(0,n-1)),a(b,nalloc,nalloc))
                        do k=1,bf
                            read(*,*) f(k,:)
                        end do
                        do k=1,bs
                            read(*,*) survival(k,:)
                        end do
                        call leslie(f,survival,a,status)
                        write(*,'(i0)') status
                        if (status == 0) then
                            do k=1,b
                                do i=1,n
                                    do j=1,n
                                        call emit(cmplx(a(k,i,j),kind=dp))
                                    end do
                                end do
                            end do
                        end if
                    end block
                end if
            case ('I')
                block
                    integer(int64), allocatable :: f(:), survival(:), a(:,:)
                    allocate(f(nalloc),survival(max(0,n-1)),a(nalloc,nalloc))
                    read(*,*) f
                    read(*,*) survival
                    call leslie(f,survival,a,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do i=1,n
                            do j=1,n
                                call emit(cmplx(real(a(i,j),dp),0.0_dp,dp))
                            end do
                        end do
                    end if
                end block
            case default
                error stop 'unknown Leslie dtype'
            end select
        case default
            error stop 'unknown mode'
        end select
    end do

contains
    subroutine emit(value)
        complex(dp), intent(in) :: value !! printed as real and imaginary double-precision values
        write(*,'(2(es26.17e3,1x))') real(value,dp),aimag(value)
    end subroutine emit
end program linalg_hadamard_leslie_dft_driver
