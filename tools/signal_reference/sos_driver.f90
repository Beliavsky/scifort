! SPDX-License-Identifier: MIT
! Test-only driver for optional comparisons against scipy.signal.
program compare_signal_sos
    use scifort_kinds, only : dp
    use scifort_signal_sos, only : sosfilt, sosfilt_zi
    use scifort_signal_filtfilt, only : filtfilt, sosfiltfilt
    implicit none
    character(len=1) :: tag
    character(len=20) :: mode
    real(dp), allocatable :: sos(:,:), b(:), a(:), x(:), zi(:,:), zf(:,:), y(:)
    complex(dp), allocatable :: csos(:,:), cb(:), ca(:), cx(:), czi(:,:), czf(:,:), cy(:)
    integer :: ios, i, ns, n, na, nb, pad, stat
    do
        read(*,*,iostat=ios) tag
        if (ios /= 0) exit
        select case(tag)
        case('R','C')
            read(*,*) ns,n
            if (tag == 'R') then
                allocate(sos(ns,6), x(n), zi(ns,2))
                do i=1,ns
                    read(*,*) sos(i,:)
                end do
                if (n > 0) read(*,*) x
                do i=1,ns
                    read(*,*) zi(i,:)
                end do
                y = sosfilt(sos,x,zi=zi,zf=zf,status=stat)
                call emit_real(stat,[y,reshape(zf,[size(zf)])])
                deallocate(sos,x,zi)
            else
                allocate(csos(ns,6), cx(n), czi(ns,2))
                do i=1,ns
                    read(*,*) csos(i,:)
                end do
                if (n > 0) read(*,*) cx
                do i=1,ns
                    read(*,*) czi(i,:)
                end do
                cy = sosfilt(csos,cx,zi=czi,zf=czf,status=stat)
                call emit_complex(stat,[cy,reshape(czf,[size(czf)])])
                deallocate(csos,cx,czi)
            end if
        case('Z','Y')
            read(*,*) ns
            if (tag=='Z') then
                allocate(sos(ns,6))
                do i=1,ns
                    read(*,*) sos(i,:)
                end do
                zi = sosfilt_zi(sos,status=stat)
                call emit_real(stat,reshape(zi,[size(zi)]))
                deallocate(sos)
            else
                allocate(csos(ns,6))
                do i=1,ns
                    read(*,*) csos(i,:)
                end do
                czi = sosfilt_zi(csos,status=stat)
                call emit_complex(stat,reshape(czi,[size(czi)]))
                deallocate(csos)
            end if
        case('F','G')
            read(*,*) nb,na,n,mode,pad
            if (tag=='F') then
                allocate(b(nb),a(na),x(n))
                read(*,*) b
                read(*,*) a
                read(*,*) x
                if (pad == -999) then
                    y = filtfilt(b,a,x,padtype=trim(mode),status=stat)
                else
                    y = filtfilt(b,a,x,padtype=trim(mode),padlen=pad,status=stat)
                end if
                call emit_real(stat,y)
                deallocate(b,a,x)
            else
                allocate(cb(nb),ca(na),cx(n))
                read(*,*) cb
                read(*,*) ca
                read(*,*) cx
                if (pad == -999) then
                    cy = filtfilt(cb,ca,cx,padtype=trim(mode),status=stat)
                else
                    cy = filtfilt(cb,ca,cx,padtype=trim(mode),padlen=pad,status=stat)
                end if
                call emit_complex(stat,cy)
                deallocate(cb,ca,cx)
            end if
        case('S','T')
            read(*,*) ns,n,mode,pad
            if (tag=='S') then
                allocate(sos(ns,6),x(n))
                do i=1,ns
                    read(*,*) sos(i,:)
                end do
                read(*,*) x
                if (pad == -999) then
                    y = sosfiltfilt(sos,x,padtype=trim(mode),status=stat)
                else
                    y = sosfiltfilt(sos,x,padtype=trim(mode),padlen=pad,status=stat)
                end if
                call emit_real(stat,y)
                deallocate(sos,x)
            else
                allocate(csos(ns,6),cx(n))
                do i=1,ns
                    read(*,*) csos(i,:)
                end do
                read(*,*) cx
                if (pad == -999) then
                    cy = sosfiltfilt(csos,cx,padtype=trim(mode),status=stat)
                else
                    cy = sosfiltfilt(csos,cx,padtype=trim(mode),padlen=pad,status=stat)
                end if
                call emit_complex(stat,cy)
                deallocate(csos,cx)
            end if
        case default
            error stop 'unknown SOS command'
        end select
    end do
contains
    subroutine emit_real(code, values)
        integer, intent(in) :: code !! return status
        real(dp), intent(in) :: values(:) !! values
        integer :: j
        write(*,'(i0,1x,i0)',advance='no') code,size(values)
        do j=1,size(values)
            write(*,'(1x,es26.17e3)',advance='no') values(j)
        end do
        write(*,*)
    end subroutine
    subroutine emit_complex(code, values)
        integer, intent(in) :: code !! return status
        complex(dp), intent(in) :: values(:) !! values
        integer :: j
        write(*,'(i0,1x,i0)',advance='no') code,size(values)
        do j=1,size(values)
            write(*,'(1x,es26.17e3,1x,es26.17e3)',advance='no') &
                real(values(j),dp),aimag(values(j))
        end do
        write(*,*)
    end subroutine
end program
