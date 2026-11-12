! SPDX-License-Identifier: MIT
! Text-driven SciPy comparison for logarithmic-time MLS skip and rewind.
program sequence_skip_driver
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_signal_sequences, only : max_len_seq, max_len_seq_result
    implicit none

    integer :: nc, nbits, n, ntaps, i, ios
    integer(int64) :: skip
    integer, allocatable :: state(:), taps(:)
    type(max_len_seq_result) :: result

    read(*,*,iostat=ios) nc
    if (ios /= 0) error stop 'missing case count'
    do i = 1, nc
        read(*,*,iostat=ios) nbits, n, skip, ntaps
        if (ios /= 0) error stop 'invalid case header'
        if (nbits < 2 .or. nbits > 32 .or. n < 0 .or. ntaps < 0) error stop 'invalid driver case'
        allocate(state(nbits), taps(ntaps))
        read(*,*,iostat=ios) state
        if (ios /= 0) error stop 'invalid initial state'
        if (ntaps > 0) then
            read(*,*,iostat=ios) taps
            if (ios /= 0) error stop 'invalid taps'
            result = max_len_seq(nbits, state=state, length=n, taps=taps, skip=skip)
        else
            result = max_len_seq(nbits, state=state, length=n, skip=skip)
        end if
        write(*,'(*(i0,1x))') result%status, size(result%seq), size(result%state)
        if (size(result%seq) > 0) then
            write(*,'(*(i0,1x))') result%seq
        else
            print '(a)', ''
        end if
        if (size(result%state) > 0) then
            write(*,'(*(i0,1x))') result%state
        else
            print '(a)', ''
        end if
        deallocate(state, taps)
    end do
end program sequence_skip_driver
