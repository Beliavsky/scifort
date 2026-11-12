! SPDX-License-Identifier: MIT
! Self-contained tests for Fibonacci maximum-length binary sequences.
program test_signal_sequences
    use scifort_signal_sequences, only : max_len_seq, max_len_seq_result, &
        sequence_status_success, sequence_status_invalid_input
    implicit none

    type(max_len_seq_result) :: a, b, c
    integer :: n, period, lag
    integer, allocatable :: pm(:)

    a = max_len_seq(4)
    call successful(a)
    call assert_bits(a%seq, [1,1,1,1,0,1,0,1,1,0,0,1,0,0,0])
    call assert_bits(a%state, [1,1,1,1])

    a = max_len_seq(8, length=20)
    call successful(a)
    call assert_bits(a%seq, [1,1,1,1,1,1,1,1,0,1,1,0,1,1,0,0,1,1,1,1])
    a = max_len_seq(8, state=[0,1,0,1,0,1,0,1],length=25,taps=[7,6,1])
    call successful(a)
    call assert_bits(a%seq, [0,1,0,1,0,1,0,1,0,0,1,0,0,1,0,0,0,1,0,1,0,0,0,0,1])
    call assert_bits(a%state, [0,0,1,1,0,1,0,0])
    b = max_len_seq(8, state=[0,5,0,3,0,1,0,8],length=25,taps=[1,6,7])
    call successful(b)
    call assert_bits(a%seq,b%seq)
    call assert_bits(a%state,b%state)

    ! Streaming a sequence in blocks must preserve both bits and final state.
    a = max_len_seq(8,length=255)
    b = max_len_seq(8,length=63)
    c = max_len_seq(8,state=b%state,length=192)
    call successful(a)
    call successful(b)
    call successful(c)
    call assert_bits(a%seq(1:63),b%seq)
    call assert_bits(a%seq(64:255),c%seq)
    call assert_bits(a%state,c%state)

    ! A full primitive cycle visits every nonzero state and has two-level
    ! periodic autocorrelation after converting {0,1} to {-1,+1}.
    do n = 2, 10
        period = 2**n-1
        a = max_len_seq(n)
        call successful(a)
        if (size(a%seq) /= period) error stop 'sequence full cycle length'
        if (count(a%seq == 1) /= 2**(n-1)) error stop 'sequence one count'
        if (count(a%seq == 0) /= 2**(n-1)-1) error stop 'sequence zero count'
        if (any(a%state /= 1)) error stop 'sequence return to initial state'
        pm = 2*a%seq-1
        if (sum(pm*pm) /= period) error stop 'MLS autocorr zero lag'
        do lag = 1, period-1
            if (sum(pm * cshift(pm,lag)) /= -1) error stop 'MLS autocorrelation'
        end do
    end do

    ! All supported default feedback polynomials reproduce an MLS prefix
    ! and keep their feedback shift register binary.
    do n = 2, 32
        a = max_len_seq(n,length=3*n+11)
        call successful(a)
        if (any(a%seq(1:n) /= 1)) error stop 'sequence prefix'
        if (any(a%state < 0) .or. any(a%state > 1)) error stop 'sequence final bits'
        if (size(a%state) /= n) error stop 'sequence final state size'
    end do
    a = max_len_seq(32)
    if (a%status /= sequence_status_invalid_input) error stop 'unsafe default length'
    a = max_len_seq(20,length=0)
    call successful(a)
    if (size(a%seq) /= 0 .or. any(a%state /= 1)) error stop 'zero length sequence'
    a = max_len_seq(4,length=0,state=[1,0,0,0])
    call successful(a)
    call assert_bits(a%state,[1,0,0,0])

    a = max_len_seq(1)
    call invalid(a)
    a = max_len_seq(33,length=4)
    call invalid(a)
    a = max_len_seq(5,state=[1,1,1])
    call invalid(a)
    a = max_len_seq(5,state=[0,0,0,0,0])
    call invalid(a)
    a = max_len_seq(5,length=-1)
    call invalid(a)
    a = max_len_seq(5,taps=[integer ::])
    call invalid(a)
    a = max_len_seq(5,taps=[5])
    call invalid(a)
    a = max_len_seq(5,taps=[-1])
    call invalid(a)

    print '(a)', 'test_signal_sequences passed'
contains
    subroutine successful(x)
        type(max_len_seq_result), intent(in) :: x !! result expected to be valid
        if (x%status /= sequence_status_success) error stop 'MLS expected success'
    end subroutine successful

    subroutine invalid(x)
        type(max_len_seq_result), intent(in) :: x !! rejected request
        if (x%status /= sequence_status_invalid_input) error stop 'MLS expected failure'
    end subroutine invalid

    subroutine assert_bits(actual, expected)
        integer, intent(in) :: actual(:) !! actual binary sequence
        integer, intent(in) :: expected(:) !! expected binary sequence
        if (size(actual) /= size(expected)) error stop 'MLS bits size mismatch'
        if (any(actual /= expected)) then
            print '(a,*(i0,1x))', 'actual:   ',actual
            print '(a,*(i0,1x))', 'expected: ',expected
            error stop 'MLS bits mismatch'
        end if
    end subroutine assert_bits
end program test_signal_sequences
