! SPDX-License-Identifier: MIT
! Compare independently against the installed SciPy version (driver only).
program periodic_sequences_driver
    use scifort_kinds, only : dp
    use scifort_signal_periodic, only : sawtooth, square
    use scifort_signal_sequences, only : max_len_seq, max_len_seq_result
    implicit none

    character(len=128) :: header
    integer :: ios, op, n, count, ntaps, status
    real(dp) :: control
    real(dp), allocatable :: t(:), width(:), y(:)
    integer, allocatable :: initial(:), taps(:)
    type(max_len_seq_result) :: mls

    do
        read(*,'(a)',iostat=ios) header
        if (ios /= 0) exit
        read(header,*) op, n, control
        if (op <= 4) then
            allocate(t(n))
            if (n > 0) read(*,*) t
            if (op == 2 .or. op == 4) then
                allocate(width(n))
                if (n > 0) read(*,*) width
            end if
            select case (op)
            case (1)
                y = sawtooth(t,width=control,status=status)
            case (2)
                y = sawtooth(t,width=width,status=status)
            case (3)
                y = square(t,duty=control,status=status)
            case (4)
                y = square(t,duty=width,status=status)
            end select
            print '(i0,1x,i0)', status,size(y)
            print '(*(es25.16,1x))', y
            deallocate(t)
            if (allocated(width)) deallocate(width)
            deallocate(y)
        else if (op == 5) then
            count = int(control)
            allocate(initial(n))
            read(*,*) initial
            read(*,*) ntaps
            if (ntaps > 0) then
                allocate(taps(ntaps))
                read(*,*) taps
                mls = max_len_seq(n,state=initial,length=count,taps=taps)
                deallocate(taps)
            else
                mls = max_len_seq(n,state=initial,length=count)
            end if
            print '(i0,1x,i0,1x,i0)', mls%status,size(mls%seq),size(mls%state)
            print '(*(i0,1x))', mls%seq
            print '(*(i0,1x))', mls%state
            deallocate(initial)
        else
            error stop 'unknown driver operation'
        end if
    end do
end program periodic_sequences_driver
