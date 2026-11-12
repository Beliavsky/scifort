! SPDX-License-Identifier: MIT
! Bulk reference driver for the assumed-rank overlap-add interface.
program oaconvolve_assumed_rank_driver
    use, intrinsic :: iso_fortran_env, only : input_unit,output_unit
    use scifort_kinds, only : dp
    use scifort_signal_fftconvolve, only : oaconvolve_nd
    implicit none
    integer :: n(15),m(15),dims(15), axes(15), r, kindflag, imode, naxes
    integer :: status, i, ios
    real(dp) :: re,im
    character(5) :: selected_mode
    complex(dp), allocatable :: z1(:), z2(:), zr(:)
    real(dp), allocatable :: a(:),b(:)

    do
        read(input_unit,*,iostat=ios) r,kindflag,imode,naxes
        if (ios /= 0) exit
        if (r < 1 .or. r > 15) error stop 'driver rank'
        read(input_unit,*) n(:r)
        read(input_unit,*) m(:r)
        read(input_unit,*) dims(:r)
        read(input_unit,*) axes(:naxes)
        select case (imode)
        case (1)
            selected_mode='full'
        case (2)
            selected_mode='same'
        case (3)
            selected_mode='valid'
        case default
            error stop 'mode'
        end select
        allocate(z1(product(n(:r))),z2(product(m(:r))),a(product(n(:r))),b(product(m(:r))))
        do i=1,size(z1)
            read(input_unit,*) re,im
            z1(i) = cmplx(re,im,dp)
        end do
        do i=1,size(z2)
            read(input_unit,*) re,im
            z2(i) = cmplx(re,im,dp)
        end do
        a = real(z1,dp)
        b = real(z2,dp)
        select case(r)
        case (1)
            block
                real(dp), allocatable :: yr(:)
                complex(dp), allocatable :: yc(:)
                allocate(yr(dims(1)))
                allocate(yc(dims(1)))
                select case (kindflag)
                case (0)
                    call oaconvolve_nd(reshape(a,[n(1)]),reshape(b,[m(1)]), &
                        yr,selected_mode,axes(:naxes),status)
                    zr = cmplx(reshape(yr,[size(yr)]),0.0_dp,dp)
                case (1)
                    call oaconvolve_nd(reshape(z1,[n(1)]),reshape(z2,[m(1)]), &
                        yc,selected_mode,axes(:naxes),status)
                    zr = reshape(yc,[size(yc)])
                case (2)
                    call oaconvolve_nd(reshape(a,[n(1)]),reshape(z2,[m(1)]), &
                        yc,selected_mode,axes(:naxes),status)
                    zr = reshape(yc,[size(yc)])
                case (3)
                    call oaconvolve_nd(reshape(z1,[n(1)]),reshape(b,[m(1)]), &
                        yc,selected_mode,axes(:naxes),status)
                    zr = reshape(yc,[size(yc)])
                end select
            end block
        case (2)
            block
                real(dp), allocatable :: yr(:,:)
                complex(dp), allocatable :: yc(:,:)
                allocate(yr(dims(1),dims(2)))
                allocate(yc(dims(1),dims(2)))
                select case (kindflag)
                case (0)
                    call oaconvolve_nd(reshape(a,[n(1),n(2)]),reshape(b,[m(1),m(2)]), &
                        yr,selected_mode,axes(:naxes),status)
                    zr = cmplx(reshape(yr,[size(yr)]),0.0_dp,dp)
                case (1)
                    call oaconvolve_nd(reshape(z1,[n(1),n(2)]),reshape(z2,[m(1),m(2)]), &
                        yc,selected_mode,axes(:naxes),status)
                    zr = reshape(yc,[size(yc)])
                case (2)
                    call oaconvolve_nd(reshape(a,[n(1),n(2)]),reshape(z2,[m(1),m(2)]), &
                        yc,selected_mode,axes(:naxes),status)
                    zr = reshape(yc,[size(yc)])
                case (3)
                    call oaconvolve_nd(reshape(z1,[n(1),n(2)]),reshape(b,[m(1),m(2)]), &
                        yc,selected_mode,axes(:naxes),status)
                    zr = reshape(yc,[size(yc)])
                end select
            end block
        case (3)
            block
                real(dp), allocatable :: yr(:,:,:)
                complex(dp), allocatable :: yc(:,:,:)
                allocate(yr(dims(1),dims(2),dims(3)))
                allocate(yc(dims(1),dims(2),dims(3)))
                select case (kindflag)
                case (0)
                    call oaconvolve_nd(reshape(a,[n(1),n(2),n(3)]),reshape(b,[m(1),m(2),m(3)]), &
                        yr,selected_mode,axes(:naxes),status)
                    zr = cmplx(reshape(yr,[size(yr)]),0.0_dp,dp)
                case (1)
                    call oaconvolve_nd(reshape(z1,[n(1),n(2),n(3)]),reshape(z2,[m(1),m(2),m(3)]), &
                        yc,selected_mode,axes(:naxes),status)
                    zr = reshape(yc,[size(yc)])
                case (2)
                    call oaconvolve_nd(reshape(a,[n(1),n(2),n(3)]),reshape(z2,[m(1),m(2),m(3)]), &
                        yc,selected_mode,axes(:naxes),status)
                    zr = reshape(yc,[size(yc)])
                case (3)
                    call oaconvolve_nd(reshape(z1,[n(1),n(2),n(3)]),reshape(b,[m(1),m(2),m(3)]), &
                        yc,selected_mode,axes(:naxes),status)
                    zr = reshape(yc,[size(yc)])
                end select
            end block
        case (5)
            block
                real(dp), allocatable :: yr(:,:,:,:,:)
                complex(dp), allocatable :: yc(:,:,:,:,:)
                allocate(yr(dims(1),dims(2),dims(3),dims(4),dims(5)))
                allocate(yc(dims(1),dims(2),dims(3),dims(4),dims(5)))
                select case (kindflag)
                case (0)
                    call oaconvolve_nd(reshape(a,[n(1),n(2),n(3),n(4),n(5)]),reshape(b,[m(1),m(2),m(3),m(4),m(5)]), &
                        yr,selected_mode,axes(:naxes),status)
                    zr = cmplx(reshape(yr,[size(yr)]),0.0_dp,dp)
                case (1)
                    call oaconvolve_nd(reshape(z1,[n(1),n(2),n(3),n(4),n(5)]),reshape(z2,[m(1),m(2),m(3),m(4),m(5)]), &
                        yc,selected_mode,axes(:naxes),status)
                    zr = reshape(yc,[size(yc)])
                case (2)
                    call oaconvolve_nd(reshape(a,[n(1),n(2),n(3),n(4),n(5)]),reshape(z2,[m(1),m(2),m(3),m(4),m(5)]), &
                        yc,selected_mode,axes(:naxes),status)
                    zr = reshape(yc,[size(yc)])
                case (3)
                    call oaconvolve_nd(reshape(z1,[n(1),n(2),n(3),n(4),n(5)]),reshape(b,[m(1),m(2),m(3),m(4),m(5)]), &
                        yc,selected_mode,axes(:naxes),status)
                    zr = reshape(yc,[size(yc)])
                end select
            end block
        case (15)
            block
                real(dp), allocatable :: yr(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
                complex(dp), allocatable :: yc(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
                allocate(yr(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6),dims(7),dims(8), &
                  dims(9),dims(10),dims(11),dims(12),dims(13),dims(14),dims(15)))
                allocate(yc(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6),dims(7),dims(8), &
                  dims(9),dims(10),dims(11),dims(12),dims(13),dims(14),dims(15)))
                select case (kindflag)
                case (0)
                    call oaconvolve_nd(reshape(a,[n]),reshape(b,[m]), &
                        yr,selected_mode,axes(:naxes),status)
                    zr = cmplx(reshape(yr,[size(yr)]),0.0_dp,dp)
                case (1)
                    call oaconvolve_nd(reshape(z1,[n]),reshape(z2,[m]), &
                        yc,selected_mode,axes(:naxes),status)
                    zr = reshape(yc,[size(yc)])
                case (2)
                    call oaconvolve_nd(reshape(a,[n]),reshape(z2,[m]), &
                        yc,selected_mode,axes(:naxes),status)
                    zr = reshape(yc,[size(yc)])
                case (3)
                    call oaconvolve_nd(reshape(z1,[n]),reshape(b,[m]), &
                        yc,selected_mode,axes(:naxes),status)
                    zr = reshape(yc,[size(yc)])
                end select
            end block
        end select
        write(output_unit,'(i0,1x,i0)') status, size(zr)
        do i=1,size(zr)
            write(output_unit,'(2es26.17)') real(zr(i),dp),aimag(zr(i))
        end do
        deallocate(z1,z2,a,b,zr)
    end do
end program oaconvolve_assumed_rank_driver

