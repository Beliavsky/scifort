! SPDX-License-Identifier: MIT
! Benchmark copied versus indexed typed elementwise callback payloads.
! Build with: gfortran -O2 -std=f2018 -ffree-line-length-none ...
program benchmark_elementwise_typed_args
    use scifort_kinds, only : dp
    use scifort_optimize_elementwise_batch, only : elementwise_typed_argument, &
        elementwise_typed_argument_from
    use scifort_optimize_root_scalar, only : find_root_vectorized_args_into, &
        find_root_vectorized_indexed_args_into, bracket_root_vectorized_args_into, &
        bracket_root_vectorized_indexed_args_into, find_root_result, bracket_root_result
    use scifort_optimize_minimize_scalar, only : find_minimum_vectorized_args_into, &
        find_minimum_vectorized_indexed_args_into, bracket_minimum_vectorized_args_into, &
        bracket_minimum_vectorized_indexed_args_into, find_minimum_result, bracket_minimum_result
    implicit none

    type :: context
        real(dp) :: target
        real(dp), allocatable :: payload(:)
    end type context

    type(context), allocatable :: source(:)
    type(elementwise_typed_argument) :: args(1)
    type(find_root_result), allocatable :: roots(:)
    type(bracket_root_result), allocatable :: root_brackets(:)
    type(find_minimum_result), allocatable :: minima(:)
    type(bracket_minimum_result), allocatable :: minimum_brackets(:)
    real(dp), allocatable :: root_input(:,:), minimum_input(:,:), xl(:), xm(:)
    real(dp) :: t0, t1, copied_time, indexed_time, copied_checksum, indexed_checksum
    integer :: n, payload_len, repeats, i, j, run, method, copied_work, indexed_work
    logical :: copied_ok, indexed_ok
    character(len=20) :: name(4)

    n = 256
    payload_len = 512
    repeats = 3
    call integer_option(1, n)
    call integer_option(2, payload_len)
    call integer_option(3, repeats)
    if (n < 1 .or. payload_len < 1 .or. repeats < 1) error stop 'Positive benchmark dimensions required'

    allocate(source(n), roots(n), root_brackets(n), minima(n), minimum_brackets(n))
    allocate(root_input(2,n), minimum_input(3,n), xl(n), xm(n))
    do i = 1, n
        source(i)%target = 0.05_dp + 0.50_dp * real(i,dp)/real(n,dp)
        allocate(source(i)%payload(payload_len))
        source(i)%payload = real(i,dp) / real(n,dp)
    end do
    args(1) = elementwise_typed_argument_from(source)
    root_input(1,:) = -2.0_dp
    root_input(2,:) = 2.0_dp
    minimum_input(1,:) = -2.0_dp
    minimum_input(2,:) = 0.0_dp
    minimum_input(3,:) = 2.0_dp
    xl = -2.0_dp
    xm = -1.5_dp
    name = [character(len=20) :: 'find_root', 'find_minimum', 'bracket_root', 'bracket_minimum']

    print '(a,i0,a,i0,a,i0)', 'lanes=', n, ' payload_reals=', payload_len, ' repeats=', repeats
    print '(a)', 'Times are process CPU seconds, excluding one-time argument construction.'
    print '(a)', 'operation            copied_s   indexed_s   copied/indexed    nfev'
    do method = 1, 4
        call run_case(method, .false., copied_checksum, copied_work, copied_ok)
        call run_case(method, .true., indexed_checksum, indexed_work, indexed_ok)
        call check_parity()
        copied_time = 0.0_dp
        indexed_time = 0.0_dp
        do run = 1, repeats
            if (mod(run,2) == 0) then
                call timed_case(method, .true., indexed_time, indexed_checksum, indexed_work, indexed_ok)
                call timed_case(method, .false., copied_time, copied_checksum, copied_work, copied_ok)
            else
                call timed_case(method, .false., copied_time, copied_checksum, copied_work, copied_ok)
                call timed_case(method, .true., indexed_time, indexed_checksum, indexed_work, indexed_ok)
            end if
            call check_parity()
        end do
        copied_time = copied_time / real(repeats,dp)
        indexed_time = indexed_time / real(repeats,dp)
        print '(a20,2f11.5,f16.2,i10)', name(method), copied_time, indexed_time, &
            copied_time/max(indexed_time,tiny(1.0_dp)), copied_work
    end do
contains
    subroutine integer_option(which, value)
        integer, intent(in) :: which !! one-based command line position
        integer, intent(inout) :: value !! default or requested positive integer
        integer :: ierr
        character(len=64) :: raw
        if (command_argument_count() < which) return
        call get_command_argument(which, raw)
        read(raw,*,iostat=ierr) value
        if (ierr /= 0) error stop 'Invalid integer benchmark option'
    end subroutine integer_option

    subroutine check_parity()
        if (.not. copied_ok .or. .not. indexed_ok) error stop 'Unconverged benchmark solve'
        if (copied_work /= indexed_work) error stop 'Evaluation count mismatch'
        if (abs(copied_checksum-indexed_checksum) > &
                1.0e-10_dp * max(1.0_dp,abs(copied_checksum))) error stop 'Numerical mismatch'
    end subroutine check_parity

    subroutine timed_case(which, indexed, seconds, checksum, work, ok)
        integer, intent(in) :: which !! solver selector
        logical, intent(in) :: indexed !! indexed callback rather than copying
        real(dp), intent(inout) :: seconds !! accumulated CPU time
        real(dp), intent(out) :: checksum !! aggregate numerical result
        integer, intent(out) :: work !! total objective evaluation count
        logical, intent(out) :: ok !! all solves converged
        call cpu_time(t0)
        call run_case(which, indexed, checksum, work, ok)
        call cpu_time(t1)
        seconds = seconds + (t1 - t0)
    end subroutine timed_case

    subroutine run_case(which, indexed, checksum, work, ok)
        integer, intent(in) :: which !! solver selector
        logical, intent(in) :: indexed !! indexed callback rather than copying
        real(dp), intent(out) :: checksum !! aggregate numerical result
        integer, intent(out) :: work !! total objective evaluation count
        logical, intent(out) :: ok !! all solves converged
        select case (which)
        case (1)
            if (indexed) then
                call find_root_vectorized_indexed_args_into(root_indexed, root_input, args, roots)
            else
                call find_root_vectorized_args_into(root_copied, root_input, args, roots)
            end if
            checksum = sum(roots%x)
            work = sum(roots%nfev)
            ok = all(roots%success)
        case (2)
            if (indexed) then
                call find_minimum_vectorized_indexed_args_into(minimum_indexed, minimum_input, args, minima)
            else
                call find_minimum_vectorized_args_into(minimum_copied, minimum_input, args, minima)
            end if
            checksum = sum(minima%x)
            work = sum(minima%nfev)
            ok = all(minima%success)
        case (3)
            if (indexed) then
                call bracket_root_vectorized_indexed_args_into(root_indexed, xl, args, root_brackets, xr0=-1.5_dp)
            else
                call bracket_root_vectorized_args_into(root_copied, xl, args, root_brackets, xr0=-1.5_dp)
            end if
            checksum = sum(root_brackets%xl) + sum(root_brackets%xr)
            work = sum(root_brackets%nfev)
            ok = all(root_brackets%success)
        case (4)
            if (indexed) then
                call bracket_minimum_vectorized_indexed_args_into(minimum_indexed, xm, args, minimum_brackets)
            else
                call bracket_minimum_vectorized_args_into(minimum_copied, xm, args, minimum_brackets)
            end if
            checksum = sum(minimum_brackets%xl) + sum(minimum_brackets%xm) + sum(minimum_brackets%xr)
            work = sum(minimum_brackets%nfev)
            ok = all(minimum_brackets%success)
        end select
    end subroutine run_case

    subroutine root_indexed(x, inputs, indices, fx)
        real(dp), intent(in) :: x(:) !! packed points
        type(elementwise_typed_argument), intent(in) :: inputs(:) !! original arguments
        integer, intent(in) :: indices(:,:) !! original argument offsets
        real(dp), intent(out) :: fx(:) !! objective values
        integer :: k
        do k = 1, size(x)
            select type(v => inputs(1)%values(indices(1,k))%value)
            type is (context)
                fx(k) = x(k) - v%target
            class default
                error stop 'Incorrect indexed type'
            end select
        end do
    end subroutine root_indexed

    subroutine root_copied(x, inputs, fx)
        real(dp), intent(in) :: x(:) !! packed points
        type(elementwise_typed_argument), intent(in) :: inputs(:) !! copied parameters
        real(dp), intent(out) :: fx(:) !! objective values
        integer :: k
        do k = 1, size(x)
            select type(v => inputs(1)%values(k)%value)
            type is (context)
                fx(k) = x(k) - v%target
            class default
                error stop 'Incorrect copied type'
            end select
        end do
    end subroutine root_copied

    subroutine minimum_indexed(x, inputs, indices, fx)
        real(dp), intent(in) :: x(:) !! packed points
        type(elementwise_typed_argument), intent(in) :: inputs(:) !! original arguments
        integer, intent(in) :: indices(:,:) !! original argument offsets
        real(dp), intent(out) :: fx(:) !! objective values
        integer :: k
        real(dp) :: delta
        do k = 1, size(x)
            select type(v => inputs(1)%values(indices(1,k))%value)
            type is (context)
                delta = x(k) - v%target
                fx(k) = delta*delta
            class default
                error stop 'Incorrect indexed type'
            end select
        end do
    end subroutine minimum_indexed

    subroutine minimum_copied(x, inputs, fx)
        real(dp), intent(in) :: x(:) !! packed points
        type(elementwise_typed_argument), intent(in) :: inputs(:) !! copied parameters
        real(dp), intent(out) :: fx(:) !! objective values
        integer :: k
        real(dp) :: delta
        do k = 1, size(x)
            select type(v => inputs(1)%values(k)%value)
            type is (context)
                delta = x(k) - v%target
                fx(k) = delta*delta
            class default
                error stop 'Incorrect copied type'
            end select
        end do
    end subroutine minimum_copied
end program benchmark_elementwise_typed_args
