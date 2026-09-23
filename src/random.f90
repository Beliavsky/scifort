! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_random
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_kinds, only : dp
    implicit none
    private

    integer(int64), parameter :: m1 = 4294967087_int64
    integer(int64), parameter :: m2 = 4294944443_int64
    integer(int64), parameter :: a12 = 1403580_int64
    integer(int64), parameter :: a13n = 810728_int64
    integer(int64), parameter :: a21 = 527612_int64
    integer(int64), parameter :: a23n = 1370589_int64
    integer(int64), parameter :: bit26 = 67108864_int64
    integer(int64), parameter :: accept26 = (m1 / bit26) * bit26
    integer(int64), parameter :: bucket_width26 = accept26 / bit26

    integer, parameter, public :: rng_status_ok = 0
    integer, parameter, public :: rng_status_invalid_state = 1

    type, public :: rng_state
        private
        integer(int64) :: s1(3) = [12345_int64, 12345_int64, 12345_int64]
        integer(int64) :: s2(3) = [12345_int64, 12345_int64, 12345_int64]
    end type rng_state

    interface rng_seed
        module procedure rng_seed_default
        module procedure rng_seed_int64
    end interface rng_seed

    public :: rng_fill_uniform
    public :: rng_get_state
    public :: rng_seed
    public :: rng_set_state
    public :: rng_uniform

contains

    subroutine rng_seed_default(state, seed)
        type(rng_state), intent(out) :: state !! generator state initialized from the seed
        integer, intent(in) :: seed !! deterministic scalar seed; any integer value is accepted

        call rng_seed_int64(state, int(seed, int64))
    end subroutine rng_seed_default

    subroutine rng_seed_int64(state, seed)
        type(rng_state), intent(out) :: state !! generator state initialized from the seed
        integer(int64), intent(in) :: seed !! deterministic scalar seed; any int64 value is accepted

        integer(int64) :: value

        value = modulo(seed, m2 - 1_int64)
        if (value == 0_int64) value = m2 - 1_int64
        state%s1 = value
        state%s2 = value
    end subroutine rng_seed_int64

    subroutine rng_get_state(state, values)
        type(rng_state), intent(in) :: state !! generator state to serialize
        integer(int64), intent(out) :: values(6) !! six recurrence components in stream order

        values(1:3) = state%s1
        values(4:6) = state%s2
    end subroutine rng_get_state

    subroutine rng_set_state(state, values, status)
        type(rng_state), intent(inout) :: state !! generator state, unchanged for invalid values
        integer(int64), intent(in) :: values(6) !! six recurrence components to restore
        integer, intent(out) :: status !! zero on success, one if the state vector is invalid

        if (.not. valid_state(values)) then
            status = rng_status_invalid_state
            return
        end if

        state%s1 = values(1:3)
        state%s2 = values(4:6)
        status = rng_status_ok
    end subroutine rng_set_state

    function rng_uniform(state) result(u)
        type(rng_state), intent(inout) :: state !! explicit generator state advanced by this draw
        real(dp) :: u

        integer(int64) :: hi
        integer(int64) :: j
        integer(int64) :: lo

        hi = next_26_bits(state)
        lo = next_26_bits(state)
        j = hi * bit26 + lo

        ! j is uniform on 0 .. 2**52 - 1. The half-step offset gives an
        ! open interval and an exactly representable binary64 grid.
        u = scale(real(j, dp), -52) + scale(1.0_dp, -53)
    end function rng_uniform

    subroutine rng_fill_uniform(state, values)
        type(rng_state), intent(inout) :: state !! explicit generator state advanced by the draws
        real(dp), intent(out) :: values(:) !! independent uniforms on the open interval (0, 1)

        integer :: i

        do i = 1, size(values)
            values(i) = rng_uniform(state)
        end do
    end subroutine rng_fill_uniform

    function next_26_bits(state) result(bits)
        type(rng_state), intent(inout) :: state !! generator state advanced until a word is accepted
        integer(int64) :: bits

        integer(int64) :: raw

        do
            raw = next_raw(state) - 1_int64
            if (raw < accept26) exit
        end do
        bits = raw / bucket_width26
    end function next_26_bits

    function next_raw(state) result(z)
        type(rng_state), intent(inout) :: state !! MRG32k3a state advanced by one recurrence step
        integer(int64) :: z

        integer(int64) :: p1
        integer(int64) :: p2

        p1 = modulo(a12 * state%s1(2) - a13n * state%s1(1), m1)
        state%s1(1) = state%s1(2)
        state%s1(2) = state%s1(3)
        state%s1(3) = p1

        p2 = modulo(a21 * state%s2(3) - a23n * state%s2(1), m2)
        state%s2(1) = state%s2(2)
        state%s2(2) = state%s2(3)
        state%s2(3) = p2

        z = p1 - p2
        if (z <= 0_int64) z = z + m1
    end function next_raw

    pure function valid_state(values) result(valid)
        integer(int64), intent(in) :: values(6) !! candidate MRG32k3a recurrence state
        logical :: valid

        valid = all(values(1:3) >= 0_int64) .and. all(values(1:3) < m1) .and. &
            all(values(4:6) >= 0_int64) .and. all(values(4:6) < m2) .and. &
            any(values(1:3) /= 0_int64) .and. any(values(4:6) /= 0_int64)
    end function valid_state

end module scifort_random
