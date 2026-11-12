! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Independently implemented periodic waves and discrete unit impulses.
! The periodic waves use their piecewise definitions over one 2*pi cycle.
module scifort_signal_periodic
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    implicit none
    private

    integer, parameter, public :: periodic_status_success = 0
    integer, parameter, public :: periodic_status_invalid_input = 1
    real(dp), parameter :: twopi = 2.0_dp*acos(-1.0_dp)

    interface sawtooth
        module procedure sawtooth_scalar
        module procedure sawtooth_array
    end interface
    interface square
        module procedure square_scalar
        module procedure square_array
    end interface
    interface unit_impulse
        module procedure unit_impulse_1d
        module procedure unit_impulse_2d
    end interface

    public :: sawtooth, square, unit_impulse

contains

    pure real(dp) function sawtooth_value(t, width) result(y)
        real(dp), intent(in) :: t !! phase in radians
        real(dp), intent(in) :: width !! fraction of period spent rising, in [0,1]
        real(dp) :: phase, cutoff

        phase = modulo(t, twopi)/twopi
        cutoff = width
        if (phase < cutoff) then
            y = 2.0_dp*phase/width - 1.0_dp
        else
            y = 1.0_dp - 2.0_dp*(phase-cutoff)/(1.0_dp-cutoff)
        end if
    end function sawtooth_value

    pure real(dp) function square_value(t, duty) result(y)
        real(dp), intent(in) :: t !! phase in radians
        real(dp), intent(in) :: duty !! fraction of cycle at +1, in [0,1]
        real(dp) :: phase

        phase = modulo(t, twopi)/twopi
        y = -1.0_dp
        if (phase < duty) y = 1.0_dp
    end function square_value

    function sawtooth_scalar(t, width, status) result(y)
        real(dp), intent(in) :: t(:) !! sample phases in radians
        real(dp), intent(in), optional :: width !! fraction rising, default 1
        integer, intent(out), optional :: status !! zero on success
        real(dp), allocatable :: y(:)
        real(dp) :: w
        integer :: i

        allocate(y(size(t)))
        y = ieee_value(0.0_dp, ieee_quiet_nan)
        if (present(status)) status = periodic_status_invalid_input
        w = 1.0_dp
        if (present(width)) w = width
        if (.not. ieee_is_finite(w)) return
        if (w < 0.0_dp .or. w > 1.0_dp) return
        if (.not. all(ieee_is_finite(t))) return
        do i = 1, size(t)
            y(i) = sawtooth_value(t(i), w)
        end do
        if (present(status)) status = periodic_status_success
    end function sawtooth_scalar

    function sawtooth_array(t, width, status) result(y)
        real(dp), intent(in) :: t(:) !! sample phases in radians
        real(dp), intent(in) :: width(:) !! rising fractions, one per sample, all in [0,1]
        integer, intent(out), optional :: status !! zero on success
        real(dp), allocatable :: y(:)
        integer :: i

        allocate(y(size(t)))
        y = ieee_value(0.0_dp, ieee_quiet_nan)
        if (present(status)) status = periodic_status_invalid_input
        if (size(t) /= size(width)) return
        if (.not. all(ieee_is_finite(t))) return
        if (.not. all(ieee_is_finite(width))) return
        if (any(width < 0.0_dp .or. width > 1.0_dp)) return
        do i = 1, size(t)
            y(i) = sawtooth_value(t(i), width(i))
        end do
        if (present(status)) status = periodic_status_success
    end function sawtooth_array

    function square_scalar(t, duty, status) result(y)
        real(dp), intent(in) :: t(:) !! sample phases in radians
        real(dp), intent(in), optional :: duty !! fraction high, default 0.5, in [0,1]
        integer, intent(out), optional :: status !! zero on success
        real(dp), allocatable :: y(:)
        real(dp) :: d
        integer :: i

        allocate(y(size(t)))
        y = ieee_value(0.0_dp, ieee_quiet_nan)
        if (present(status)) status = periodic_status_invalid_input
        d = 0.5_dp
        if (present(duty)) d = duty
        if (.not. ieee_is_finite(d)) return
        if (d < 0.0_dp .or. d > 1.0_dp) return
        if (.not. all(ieee_is_finite(t))) return
        do i = 1, size(t)
            y(i) = square_value(t(i), d)
        end do
        if (present(status)) status = periodic_status_success
    end function square_scalar

    function square_array(t, duty, status) result(y)
        real(dp), intent(in) :: t(:) !! sample phases in radians
        real(dp), intent(in) :: duty(:) !! high fractions, one per sample, all in [0,1]
        integer, intent(out), optional :: status !! zero on success
        real(dp), allocatable :: y(:)
        integer :: i

        allocate(y(size(t)))
        y = ieee_value(0.0_dp, ieee_quiet_nan)
        if (present(status)) status = periodic_status_invalid_input
        if (size(t) /= size(duty)) return
        if (.not. all(ieee_is_finite(t))) return
        if (.not. all(ieee_is_finite(duty))) return
        if (any(duty < 0.0_dp .or. duty > 1.0_dp)) return
        do i = 1, size(t)
            y(i) = square_value(t(i), duty(i))
        end do
        if (present(status)) status = periodic_status_success
    end function square_array

    function unit_impulse_1d(n, idx, centered, status) result(y)
        integer, intent(in) :: n !! number of samples, nonnegative
        integer, intent(in), optional :: idx !! impulse location, one-based, default 1
        logical, intent(in), optional :: centered !! choose floor(n/2)+1 when true
        integer, intent(out), optional :: status !! zero on success
        real(dp), allocatable :: y(:)
        integer :: k

        if (present(status)) status = periodic_status_invalid_input
        allocate(y(max(0,n)))
        y = 0.0_dp
        if (n < 0) return
        k = 1
        if (present(centered)) then
            if (centered) k = n/2+1
        end if
        if (present(idx)) k = idx
        if (n > 0) then
            if (k < 1 .or. k > n) then
                y = ieee_value(0.0_dp, ieee_quiet_nan)
                return
            end if
            y(k) = 1.0_dp
        else
            if (present(idx)) return
        end if
        if (present(status)) status = periodic_status_success
    end function unit_impulse_1d

    function unit_impulse_2d(shape, idx, centered, status) result(y)
        integer, intent(in) :: shape(2) !! row and column counts, both nonnegative
        integer, intent(in), optional :: idx(2) !! one-based row and column impulse indices
        logical, intent(in), optional :: centered !! select floor(shape/2)+1 when true
        integer, intent(out), optional :: status !! zero on success
        real(dp), allocatable :: y(:,:)
        integer :: k(2)

        if (present(status)) status = periodic_status_invalid_input
        allocate(y(max(0,shape(1)),max(0,shape(2))))
        y = 0.0_dp
        if (any(shape < 0)) return
        k = 1
        if (present(centered)) then
            if (centered) k = shape/2+1
        end if
        if (present(idx)) k = idx
        if (all(shape > 0)) then
            if (any(k < 1) .or. any(k > shape)) then
                y = ieee_value(0.0_dp, ieee_quiet_nan)
                return
            end if
            y(k(1),k(2)) = 1.0_dp
        else
            if (present(idx)) return
        end if
        if (present(status)) status = periodic_status_success
    end function unit_impulse_2d

end module scifort_signal_periodic
