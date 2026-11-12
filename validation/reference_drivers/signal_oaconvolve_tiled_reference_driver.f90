! SPDX-License-Identifier: MIT
! Stream I/O driver for cross-checking tiled overlap-add against SciPy.
program signal_oaconvolve_tiled_reference_driver
    use, intrinsic :: iso_fortran_env, only : int32, iostat_end
    use scifort_kinds, only : dp
    use scifort_signal_fftconvolve, only : oaconvolve_nd
    implicit none
    integer :: ui, uo, stat, ios, ctype
    integer(int32) :: header(9)
    integer :: n1, n2, m1, m2, q1, q2, axis_code, mode_code
    integer(int32) :: code
    character(5) :: mode
    character(1024) :: file_in, file_out
    real(dp), allocatable :: ar(:,:), br(:,:), yr(:,:)
    complex(dp), allocatable :: ac(:,:), bc(:,:), yc(:,:)

    if (command_argument_count() /= 2) error stop 'usage: driver input.bin output.bin'
    call get_command_argument(1,file_in)
    call get_command_argument(2,file_out)
    open(newunit=ui,file=trim(file_in),form='unformatted',access='stream',status='old',action='read')
    open(newunit=uo,file=trim(file_out),form='unformatted',access='stream',status='replace',action='write')
    do
        read(ui,iostat=ios) header
        if (ios == iostat_end) exit
        if (ios /= 0) error stop 'input header read error'
        n1=header(1)
        n2=header(2)
        m1=header(3)
        m2=header(4)
        q1=header(5)
        q2=header(6)
        axis_code=header(7)
        mode_code=header(8)
        ctype=header(9)
        select case(mode_code)
        case(1)
            mode='full'
        case(2)
            mode='same'
        case(3)
            mode='valid'
        case default
            error stop 'invalid mode_code'
        end select
        select case(ctype)
        case(1)
            allocate(ar(n1,n2),br(m1,m2),yr(q1,q2))
            read(ui) ar,br
            if (axis_code == 0) then
                call oaconvolve_nd(ar,br,yr,mode=trim(mode),status=stat)
            else
                call oaconvolve_nd(ar,br,yr,mode=trim(mode),axes=[axis_code],status=stat)
            end if
            code=int(stat,int32)
            write(uo) code,yr
            deallocate(ar,br,yr)
        case(2)
            allocate(ac(n1,n2),bc(m1,m2),yc(q1,q2))
            read(ui) ac,bc
            if (axis_code == 0) then
                call oaconvolve_nd(ac,bc,yc,mode=trim(mode),status=stat)
            else
                call oaconvolve_nd(ac,bc,yc,mode=trim(mode),axes=[axis_code],status=stat)
            end if
            code=int(stat,int32)
            write(uo) code,yc
            deallocate(ac,bc,yc)
        case(3)
            allocate(ar(n1,n2),bc(m1,m2),yc(q1,q2))
            read(ui) ar,bc
            if (axis_code == 0) then
                call oaconvolve_nd(ar,bc,yc,mode=trim(mode),status=stat)
            else
                call oaconvolve_nd(ar,bc,yc,mode=trim(mode),axes=[axis_code],status=stat)
            end if
            code=int(stat,int32)
            write(uo) code,yc
            deallocate(ar,bc,yc)
        case(4)
            allocate(ac(n1,n2),br(m1,m2),yc(q1,q2))
            read(ui) ac,br
            if (axis_code == 0) then
                call oaconvolve_nd(ac,br,yc,mode=trim(mode),status=stat)
            else
                call oaconvolve_nd(ac,br,yc,mode=trim(mode),axes=[axis_code],status=stat)
            end if
            code=int(stat,int32)
            write(uo) code,yc
            deallocate(ac,br,yc)
        case default
            error stop 'invalid ctype'
        end select
    end do
    close(ui)
    close(uo)
end program signal_oaconvolve_tiled_reference_driver
