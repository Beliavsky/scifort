! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Maximum-length binary sequences from Fibonacci linear-feedback registers.
! Uses the published feedback-tap convention described by SciPy max_len_seq.
! This is an independent iterative formulation, not a port of SciPy source.
module scifort_signal_sequences
    use, intrinsic :: iso_fortran_env, only : int64
    implicit none
    private

    integer, parameter, public :: sequence_status_success = 0
    integer, parameter, public :: sequence_status_invalid_input = 1
    integer, parameter, public :: sequence_status_allocation_failed = 2

    type, public :: max_len_seq_result
        integer, allocatable :: seq(:) !! generated bits, values 0 or 1
        integer, allocatable :: state(:) !! register state after generating bits
        integer :: status = sequence_status_invalid_input !! zero on success
    end type max_len_seq_result

    public :: max_len_seq

contains

    ! One Fibonacci-register clock, represented by low-order bits of an
    ! int64. The high register bit contains the XOR of the selected feedback
    ! bits. This works for every supported width, including 32 bits, without
    ! relying on signed default-integer overflow.
    pure function clock_mask(bits, feedback, nbits, reverse) result(next_bits)
        integer(int64), intent(in) :: bits !! packed current register state
        integer(int64), intent(in) :: feedback !! parity mask including bit zero
        integer, intent(in) :: nbits !! register width (2..32)
        logical, intent(in) :: reverse !! invert rather than advance a clock
        integer(int64) :: next_bits
        integer(int64) :: width_mask
        logical :: next_first

        width_mask = shiftl(1_int64, nbits) - 1_int64
        if (reverse) then
            ! After a forward step next(1:n-1)=previous(2:n).
            ! Recover previous(1) from the feedback bit and the known taps.
            next_bits = iand(shiftl(bits, 1), width_mask)
            next_first = btest(bits, nbits-1) .neqv. &
                (mod(popcnt(iand(bits, shiftr(feedback, 1))),2) == 1)
            if (next_first) next_bits = ibset(next_bits, 0)
        else
            next_bits = shiftr(bits, 1)
            if (mod(popcnt(iand(bits, feedback)),2) == 1) &
                next_bits = ibset(next_bits, nbits-1)
        end if
    end function clock_mask

    ! Apply a GF(2) linear map represented by images of the basis vectors.
    pure function apply_columns(columns, bits, nbits) result(mapped)
        integer(int64), intent(in) :: columns(32) !! binary transition columns
        integer(int64), intent(in) :: bits !! packed input register
        integer, intent(in) :: nbits !! active register width
        integer(int64) :: mapped
        integer :: j

        mapped = 0_int64
        do j = 1, nbits
            if (btest(bits, j-1)) mapped = ieor(mapped, columns(j))
        end do
    end function apply_columns

    ! Repeated squaring of a binary transition matrix. O(nbits**2 * log(skip))
    ! avoids materializing skipped bits, supports positive and negative skips,
    ! and handles the most-negative int64 without overflowing ABS(skip).
    pure subroutine jump_register(state, feedback, skip)
        integer, intent(inout) :: state(:) !! binary register, updated in place
        integer(int64), intent(in) :: feedback !! one-clock XOR mask
        integer(int64), intent(in) :: skip !! signed number of clock steps
        integer(int64) :: columns(32), squared(32), bits, remaining
        integer :: j, nbits
        logical :: reverse

        nbits = size(state)
        bits = 0_int64
        do j = 1, nbits
            if (state(j) /= 0) bits = ibset(bits,j-1)
        end do
        if (skip == 0_int64) return
        reverse = skip < 0_int64
        do j = 1, nbits
            columns(j) = clock_mask(ibset(0_int64,j-1),feedback,nbits,reverse)
        end do
        if (reverse) then
            ! One reverse clock consumes the extra step that cannot be
            ! represented as a positive int64 for skip = -2**63.
            bits = apply_columns(columns,bits,nbits)
            remaining = -(skip+1_int64)
        else
            remaining = skip
        end if
        do while (remaining > 0_int64)
            if (btest(remaining,0)) bits = apply_columns(columns,bits,nbits)
            remaining = shiftr(remaining,1)
            if (remaining == 0_int64) exit
            do j = 1, nbits
                squared(j) = apply_columns(columns,columns(j),nbits)
            end do
            columns = squared
        end do
        do j = 1, nbits
            state(j) = merge(1,0,btest(bits,j-1))
        end do
    end subroutine jump_register

    subroutine primitive_taps(nbits, taps)
        integer, intent(in) :: nbits !! shift-register width, in [2,32]
        integer, allocatable, intent(out) :: taps(:) !! default feedback locations

        ! Primitive-polynomial tap positions as documented by the public
        ! newwaveinstruments.com M-sequence table linked from SciPy's manual.
        select case (nbits)
        case (2)
            taps = [1]
        case (3)
            taps = [2]
        case (4)
            taps = [3]
        case (5)
            taps = [3]
        case (6)
            taps = [5]
        case (7)
            taps = [6]
        case (8)
            taps = [7,6,1]
        case (9)
            taps = [5]
        case (10)
            taps = [7]
        case (11)
            taps = [9]
        case (12)
            taps = [11,10,4]
        case (13)
            taps = [12,11,8]
        case (14)
            taps = [13,12,2]
        case (15)
            taps = [14]
        case (16)
            taps = [15,13,4]
        case (17)
            taps = [14]
        case (18)
            taps = [11]
        case (19)
            taps = [18,17,14]
        case (20)
            taps = [17]
        case (21)
            taps = [19]
        case (22)
            taps = [21]
        case (23)
            taps = [18]
        case (24)
            taps = [23,22,17]
        case (25)
            taps = [22]
        case (26)
            taps = [25,24,20]
        case (27)
            taps = [26,25,22]
        case (28)
            taps = [25]
        case (29)
            taps = [27]
        case (30)
            taps = [29,28,7]
        case (31)
            taps = [28]
        case (32)
            taps = [31,30,10]
        case default
            taps = [integer ::]
        end select
    end subroutine primitive_taps

    function max_len_seq(nbits, state, length, taps, skip) result(out)
        integer, intent(in) :: nbits !! number of shift-register bits, from 2 to 32
        integer, intent(in), optional :: state(:) !! initial register, nonzero entries become 1
        integer, intent(in), optional :: length !! requested sample count; full cycle for nbits<=20
        integer, intent(in), optional :: taps(:) !! SciPy-position feedback taps, 1..nbits-1
        integer(int64), intent(in), optional :: skip !! signed steps before output, including negative jumps
        type(max_len_seq_result) :: out

        integer, allocatable :: positions(:)
        logical :: active(0:32)
        integer :: bit, i, j, n, allocate_status
        integer(int64) :: feedback

        allocate(out%seq(0), out%state(0))
        if (nbits < 2 .or. nbits > 32) return
        if (present(state)) then
            if (size(state) /= nbits) return
        end if
        if (present(taps)) then
            if (size(taps) == 0) return
            if (any(taps < 1) .or. any(taps >= nbits)) return
            positions = taps
        else
            call primitive_taps(nbits, positions)
        end if
        if (present(length)) then
            if (length < 0) return
            n = length
        else
            ! Never attempt to allocate an implicit multi-gigabyte cycle.
            ! For longer registers, require the caller to specify length.
            if (nbits > 20) return
            n = 2**nbits-1
        end if

        deallocate(out%state)
        allocate(out%state(nbits))
        out%state = 1
        if (present(state)) out%state = merge(1, 0, state /= 0)
        if (all(out%state == 0)) return

        deallocate(out%seq)
        allocate(out%seq(n), stat=allocate_status)
        if (allocate_status /= 0) then
            allocate(out%seq(0))
            out%status = sequence_status_allocation_failed
            return
        end if
        ! Treat taps as a set, like the binary feedback mask. A repeated
        ! position is counted only once, rather than XORing it twice.
        active = .false.
        do j = 1, size(positions)
            active(positions(j)) = .true.
        end do
        ! Only the initial register is advanced: skip-ahead never allocates
        ! or generates the intermediate sequence. The same feedback tap set
        ! is used for the normal sample loop and the jump matrix.
        if (present(skip)) then
            if (skip /= 0_int64) then
                feedback = ibset(0_int64,0)
                do j = 1, nbits-1
                    if (active(j)) feedback = ibset(feedback,j)
                end do
                call jump_register(out%state,feedback,skip)
            end if
        end if

        do i = 1, n
            out%seq(i) = out%state(1)
            bit = out%state(1)
            do j = 1, nbits-1
                if (active(j)) bit = ieor(bit, out%state(j+1))
            end do
            out%state(1:nbits-1) = out%state(2:nbits)
            out%state(nbits) = bit
        end do
        out%status = sequence_status_success
    end function max_len_seq

end module scifort_signal_sequences
