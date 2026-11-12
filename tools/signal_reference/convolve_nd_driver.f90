! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
! stdin/stdout driver for independent SciPy convolution/correlation validation.
program convolve_nd_driver
    use scifort_kinds, only : dp
    use scifort_signal_convolve_nd, only : convolve_nd, correlate_nd, convolution_shape_nd
    implicit none
    integer :: kind, operation, algorithm, mode_id, naxes, code, i, ios, n(15), m(15)
    integer, allocatable :: dims(:), axes(:)
    character(len=6) :: modestr, methodstr
    complex(dp), allocatable :: flat_a(:), flat_b(:), flat_y(:)
    real(dp), allocatable :: ra(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
    complex(dp), allocatable :: ca(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
    real(dp), allocatable :: rb(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
    complex(dp), allocatable :: cb(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
    real(dp), allocatable :: ry(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
    complex(dp), allocatable :: cy(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
    do
        read(*,*,iostat=ios) kind, operation, algorithm, mode_id, naxes
        if (ios/=0) exit
        read(*,*) n
        read(*,*) m
        if (allocated(axes)) deallocate(axes)
        if (naxes>=0) then
            allocate(axes(naxes))
            if (naxes>0) read(*,*) axes
        end if
        allocate(flat_a(product(n)),flat_b(product(m)))
        do i=1,size(flat_a)
            read(*,*) flat_a(i)
        end do
        do i=1,size(flat_b)
            read(*,*) flat_b(i)
        end do
        select case(mode_id)
        case (1)
            modestr='full'
        case (2)
            modestr='same'
        case (3)
            modestr='valid'
        end select
        select case(algorithm)
        case(1)
            methodstr='direct'
        case(2)
            methodstr='fft'
        case(3)
            methodstr='auto'
        end select
        dims=convolution_shape_nd(n,m,trim(modestr),axes)
        if (size(dims)/=15) error stop 'invalid driver layout'
        if (kind==0) then
            allocate(ra(n(1),n(2),n(3),n(4),n(5),n(6),n(7),n(8),n(9),n(10),n(11),n(12),n(13),n(14),n(15)), &
                rb(m(1),m(2),m(3),m(4),m(5),m(6),m(7),m(8),m(9),m(10),m(11),m(12),m(13),m(14),m(15)), &
                ry(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6),dims(7),dims(8),dims(9),dims(10),dims(11), &
                dims(12),dims(13),dims(14),dims(15)))
            ra=reshape(real(flat_a,dp),[n(1),n(2),n(3),n(4),n(5),n(6),n(7),n(8),n(9),n(10),n(11),n(12),n(13), &
                n(14),n(15)])
            rb=reshape(real(flat_b,dp),[m(1),m(2),m(3),m(4),m(5),m(6),m(7),m(8),m(9),m(10),m(11),m(12),m(13), &
                m(14),m(15)])
            ry=0.0_dp
            if (operation==0) then
                call convolve_nd(ra,rb,ry,trim(modestr),trim(methodstr),axes,code)
            else
                call correlate_nd(ra,rb,ry,trim(modestr),trim(methodstr),axes,code)
            end if
            flat_y=cmplx(reshape(ry,[size(ry)]),0.0_dp,dp)
            deallocate(ra,rb,ry)
        else
            allocate(cy(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6),dims(7),dims(8),dims(9),dims(10), &
                dims(11),dims(12),dims(13),dims(14),dims(15)))
            cy=(0.0_dp,0.0_dp)
            select case(kind)
            case(1)
                allocate(ca(n(1),n(2),n(3),n(4),n(5),n(6),n(7),n(8),n(9),n(10),n(11),n(12),n(13),n(14),n(15)), &
                cb(m(1),m(2),m(3),m(4),m(5),m(6),m(7),m(8),m(9),m(10),m(11),m(12),m(13),m(14),m(15)))
                ca=reshape(flat_a,[n(1),n(2),n(3),n(4),n(5),n(6),n(7),n(8),n(9),n(10),n(11),n(12),n(13),n(14),n(15)])
                cb=reshape(flat_b,[m(1),m(2),m(3),m(4),m(5),m(6),m(7),m(8),m(9),m(10),m(11),m(12),m(13),m(14),m(15)])
                if (operation==0) then
                    call convolve_nd(ca,cb,cy,trim(modestr),trim(methodstr),axes,code)
                else
                    call correlate_nd(ca,cb,cy,trim(modestr),trim(methodstr),axes,code)
                end if
                deallocate(ca,cb)
            case(2)
                allocate(ra(n(1),n(2),n(3),n(4),n(5),n(6),n(7),n(8),n(9),n(10),n(11),n(12),n(13),n(14),n(15)), &
                cb(m(1),m(2),m(3),m(4),m(5),m(6),m(7),m(8),m(9),m(10),m(11),m(12),m(13),m(14),m(15)))
                ra=reshape(real(flat_a,dp),[n(1),n(2),n(3),n(4),n(5),n(6),n(7),n(8),n(9),n(10),n(11),n(12), &
                n(13),n(14),n(15)])
                cb=reshape(flat_b,[m(1),m(2),m(3),m(4),m(5),m(6),m(7),m(8),m(9),m(10),m(11),m(12),m(13),m(14),m(15)])
                if (operation==0) then
                    call convolve_nd(ra,cb,cy,trim(modestr),trim(methodstr),axes,code)
                else
                    call correlate_nd(ra,cb,cy,trim(modestr),trim(methodstr),axes,code)
                end if
                deallocate(ra,cb)
            case(3)
                allocate(ca(n(1),n(2),n(3),n(4),n(5),n(6),n(7),n(8),n(9),n(10),n(11),n(12),n(13),n(14),n(15)), &
                rb(m(1),m(2),m(3),m(4),m(5),m(6),m(7),m(8),m(9),m(10),m(11),m(12),m(13),m(14),m(15)))
                ca=reshape(flat_a,[n(1),n(2),n(3),n(4),n(5),n(6),n(7),n(8),n(9),n(10),n(11),n(12),n(13),n(14),n(15)])
                rb=reshape(real(flat_b,dp),[m(1),m(2),m(3),m(4),m(5),m(6),m(7),m(8),m(9),m(10),m(11),m(12), &
                m(13),m(14),m(15)])
                if (operation==0) then
                    call convolve_nd(ca,rb,cy,trim(modestr),trim(methodstr),axes,code)
                else
                    call correlate_nd(ca,rb,cy,trim(modestr),trim(methodstr),axes,code)
                end if
                deallocate(ca,rb)
            end select
            flat_y=reshape(cy,[size(cy)])
            deallocate(cy)
        end if
        print '(i0,*(1x,i0))',code,dims
        do i=1,size(flat_y)
            print '(es26.17,1x,es26.17)',real(flat_y(i),dp),aimag(flat_y(i))
        end do
        deallocate(flat_a,flat_b,flat_y)
    end do
end program convolve_nd_driver
