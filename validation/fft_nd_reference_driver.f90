! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program fft_nd_reference_driver
    use scifort_fft
    implicit none
    real(dp), allocatable :: rx(:,:,:,:), ry(:,:,:,:)
    complex(dp), allocatable :: cx(:,:,:,:), cy(:,:,:,:), flat(:)
    integer, allocatable :: lengths(:), axes(:)
    integer :: op, real_input, kind, mode, orth, ns, na, xd(4), yd(4), code, ios, j
    real(dp) :: re, im
    character(8), parameter :: norms(3) = [character(8) :: 'backward','forward','ortho']
    logical, allocatable :: orthogonal
    do
        read(*,*,iostat=ios) op,real_input,kind,mode,orth,ns,na,xd,yd
        if (ios /= 0) exit
        if (allocated(lengths)) deallocate(lengths)
        if (allocated(axes)) deallocate(axes)
        if (allocated(orthogonal)) deallocate(orthogonal)
        if (ns >= 0) then
            allocate(lengths(ns))
            if (ns > 0) read(*,*) lengths
        end if
        if (na >= 0) then
            allocate(axes(na))
            if (na > 0) read(*,*) axes
        end if
        if (orth >= 0) then
            allocate(orthogonal)
            orthogonal = orth == 1
        end if
        allocate(flat(product(xd)))
        do j = 1, size(flat)
            read(*,*) re,im
            flat(j) = cmplx(re,im,dp)
        end do
        cx = reshape(flat,xd)
        rx = real(cx,dp)
        allocate(cy(yd(1),yd(2),yd(3),yd(4)),ry(yd(1),yd(2),yd(3),yd(4)))
        select case(op)
        case(1)
            if (real_input == 1) then
                call fftn(rx,cy,s=lengths,axes=axes,norm=trim(norms(mode+1)),status=code)
            else
                call fftn(cx,cy,s=lengths,axes=axes,norm=trim(norms(mode+1)),status=code)
            end if
        case(2)
            if (real_input == 1) then
                call ifftn(rx,cy,s=lengths,axes=axes,norm=trim(norms(mode+1)),status=code)
            else
                call ifftn(cx,cy,s=lengths,axes=axes,norm=trim(norms(mode+1)),status=code)
            end if
        case(3)
            call rfftn(rx,cy,s=lengths,axes=axes,norm=trim(norms(mode+1)),status=code)
        case(4)
            if (real_input == 1) then
                call irfftn(rx,ry,s=lengths,axes=axes,norm=trim(norms(mode+1)),status=code)
            else
                call irfftn(cx,ry,s=lengths,axes=axes,norm=trim(norms(mode+1)),status=code)
            end if
            cy = cmplx(ry,0.0_dp,dp)
        case(5)
            if (real_input == 1) then
                call hfftn(rx,ry,s=lengths,axes=axes,norm=trim(norms(mode+1)),status=code)
            else
                call hfftn(cx,ry,s=lengths,axes=axes,norm=trim(norms(mode+1)),status=code)
            end if
            cy = cmplx(ry,0.0_dp,dp)
        case(6)
            call ihfftn(rx,cy,s=lengths,axes=axes,norm=trim(norms(mode+1)),status=code)
        case(7)
            if (real_input == 1) then
                call dctn(rx,ry,s=lengths,axes=axes,norm=trim(norms(mode+1)),status=code,type=kind,orthogonalize=orthogonal)
            else
                call dctn(cx,cy,s=lengths,axes=axes,norm=trim(norms(mode+1)),status=code,type=kind,orthogonalize=orthogonal)
            end if
            if (real_input == 1) cy = cmplx(ry,0.0_dp,dp)
        case(8)
            if (real_input == 1) then
                call idctn(rx,ry,s=lengths,axes=axes,norm=trim(norms(mode+1)),status=code,type=kind,orthogonalize=orthogonal)
            else
                call idctn(cx,cy,s=lengths,axes=axes,norm=trim(norms(mode+1)),status=code,type=kind,orthogonalize=orthogonal)
            end if
            if (real_input == 1) cy = cmplx(ry,0.0_dp,dp)
        case(9)
            if (real_input == 1) then
                call dstn(rx,ry,s=lengths,axes=axes,norm=trim(norms(mode+1)),status=code,type=kind,orthogonalize=orthogonal)
            else
                call dstn(cx,cy,s=lengths,axes=axes,norm=trim(norms(mode+1)),status=code,type=kind,orthogonalize=orthogonal)
            end if
            if (real_input == 1) cy = cmplx(ry,0.0_dp,dp)
        case(10)
            if (real_input == 1) then
                call idstn(rx,ry,s=lengths,axes=axes,norm=trim(norms(mode+1)),status=code,type=kind,orthogonalize=orthogonal)
            else
                call idstn(cx,cy,s=lengths,axes=axes,norm=trim(norms(mode+1)),status=code,type=kind,orthogonalize=orthogonal)
            end if
            if (real_input == 1) cy = cmplx(ry,0.0_dp,dp)
        end select
        flat = reshape(cy,[size(cy)])
        write(*,'(2(i0,1x))') code,size(flat)
        do j = 1, size(flat)
            write(*,'(2(es26.17e3,1x))') real(flat(j),dp),aimag(flat(j))
        end do
        deallocate(flat,cx,rx,cy,ry)
    end do
end program fft_nd_reference_driver
