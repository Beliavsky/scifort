! SPDX-License-Identifier: MIT
! Skip-ahead/rewind of maximum-length and arbitrary Fibonacci sequences.
program test_signal_sequences_skip
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_signal_sequences, only : max_len_seq, max_len_seq_result, &
        sequence_status_success, sequence_status_invalid_input
    implicit none

    type(max_len_seq_result) :: reference, jumped, initial, back, second, streaming
    integer, parameter :: distances(*) = [0, 1, 2, 7, 13, 67, 511, 1025]
    integer :: n, k, skip_small, cycle
    integer(int64) :: huge_skip
    integer, allocatable :: initial_state(:)

    do n = 2, 32
        do k = 1, size(distances)
            skip_small = distances(k)
            reference = max_len_seq(n, length=skip_small+35)
            jumped = max_len_seq(n, skip=int(skip_small,int64), length=35)
            call check_success(reference)
            call check_success(jumped)
            if (any(reference%seq(skip_small+1:) /= jumped%seq)) &
                error stop 'skip-ahead samples differ from direct generation'
            if (any(reference%state /= jumped%state)) error stop 'skip-ahead final state'

            initial = max_len_seq(n, skip=int(skip_small,int64), length=0)
            back = max_len_seq(n, state=initial%state, &
                skip=-int(skip_small,int64), length=0)
            call check_success(back)
            if (any(back%state /= 1)) error stop 'reverse skip returned wrong state'
        end do
    end do

    ! All default primitive polynomials return to the initial state after a
    ! full cycle. A cycle is small enough to check exactly for these widths.
    do n = 2, 16
        cycle = 2**n-1
        reference = max_len_seq(n,length=31)
        jumped = max_len_seq(n,skip=int(cycle,int64),length=31)
        if (any(reference%seq /= jumped%seq)) error stop 'primitive cycle samples'
        if (any(reference%state /= jumped%state)) error stop 'primitive cycle state'
    end do

    ! Nonprimitive and repeated/custom feedback taps, and nonbinary seeds,
    ! must obey exactly the same binary-clock convention as the old API.
    initial_state = [1, 0, 2, 0, 1, 0, 0, -1]
    reference = max_len_seq(8,state=initial_state,taps=[1,3,3,6],length=175)
    jumped = max_len_seq(8,state=initial_state,taps=[6,3,1],skip=97_int64,length=78)
    call check_success(reference)
    call check_success(jumped)
    if (any(reference%seq(98:) /= jumped%seq)) error stop 'custom taps skip'
    if (any(reference%state /= jumped%state)) error stop 'custom taps state'
    initial = max_len_seq(8,state=initial_state,taps=[6,3,1],skip=97_int64,length=0)
    back = max_len_seq(8,state=initial%state,taps=[6,3,1],skip=-97_int64,length=0)
    if (any(back%state /= merge(1,0,initial_state/=0))) error stop 'custom taps rewind'

    ! 63-bit and most-negative int64 distances must neither overflow nor
    ! perform a linear scan. Check inverse composition and state continuation.
    huge_skip = huge(0_int64)
    do n = 2, 32, 5
        initial = max_len_seq(n,length=0,skip=huge_skip)
        back = max_len_seq(n,state=initial%state,length=0,skip=-huge_skip)
        if (any(back%state /= 1)) error stop 'large positive skip inverse'
        initial = max_len_seq(n,length=0,skip=-huge_skip-1_int64)
        second = max_len_seq(n,state=initial%state,length=0,skip=huge_skip)
        back = max_len_seq(n,state=second%state,length=0,skip=1_int64)
        if (any(back%state /= 1)) error stop 'minimum int64 rewind'
        jumped = max_len_seq(n,length=35,skip=huge_skip)
        streaming = max_len_seq(n,state=initial%state,length=17,skip=huge_skip)
        ! Above 'initial' is at -2**63; so adding 2**63-1 gives -1.
        reference = max_len_seq(n,length=17,skip=-1_int64)
        if (any(streaming%seq /= reference%seq)) error stop 'huge mixed skip'
        if (any(streaming%state /= reference%state)) error stop 'huge mixed state'
    end do

    jumped = max_len_seq(4,length=0,skip=-1_int64)
    call check_success(jumped)
    if (size(jumped%seq) /= 0) error stop 'zero length after reverse skip'
    back = max_len_seq(4,state=jumped%state,length=0,skip=1_int64)
    if (any(back%state /= 1)) error stop 'zero length reverse state'

    jumped = max_len_seq(1,length=5,skip=10_int64)
    if (jumped%status /= sequence_status_invalid_input) error stop 'invalid width'
    jumped = max_len_seq(7,state=[0,0,0,0,0,0,0],length=1,skip=15_int64)
    if (jumped%status /= sequence_status_invalid_input) error stop 'zero state must fail'
    jumped = max_len_seq(7,taps=[0],length=1,skip=-15_int64)
    if (jumped%status /= sequence_status_invalid_input) error stop 'invalid taps must fail'

    print '(a)', 'test_signal_sequences_skip passed'
contains
    subroutine check_success(r)
        type(max_len_seq_result), intent(in) :: r !! expected successful sequence
        if (r%status /= sequence_status_success) error stop 'sequence status failure'
    end subroutine check_success
end program test_signal_sequences_skip
