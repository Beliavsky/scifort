! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Excursion-based peak/valley detection, inspired by the algorithm described
! by Eli Billauer and the MIT-licensed PEAKS package by Jason Christopherson.
! Independent SciFort implementation: no upstream source code copied.
! https://github.com/jchristopherson/peaks
module scifort_signal_peak_excursion
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    implicit none
    private

    integer, parameter, public :: excursion_status_success = 0
    integer, parameter, public :: excursion_status_invalid_input = 1
    public :: peak_detect

contains

    pure subroutine peak_detect(x, delta, imx, mx, imn, mn, status)
        real(dp), intent(in) :: x(:) !! sampled one-dimensional signal
        real(dp), intent(in) :: delta !! strictly positive minimum peak-to-valley excursion
        integer, allocatable, intent(out) :: imx(:) !! one-based maxima locations
        real(dp), allocatable, intent(out) :: mx(:) !! maxima amplitudes
        integer, allocatable, intent(out) :: imn(:) !! one-based minima locations
        real(dp), allocatable, intent(out) :: mn(:) !! minima amplitudes
        integer, intent(out), optional :: status !! zero on success; one on invalid input

        integer, allocatable :: tmxi(:), tmni(:)
        real(dp), allocatable :: tmx(:), tmn(:)
        integer :: i, n, nmxi, nmni, candidate_max_i, candidate_min_i
        real(dp) :: candidate_max, candidate_min
        logical :: looking_for_maximum

        allocate(imx(0), mx(0), imn(0), mn(0))
        if (present(status)) status = excursion_status_invalid_input
        if (.not. ieee_is_finite(delta)) return
        if (delta <= 0.0_dp) return
        if (any(.not. ieee_is_finite(x))) return

        n = size(x)
        if (n == 0) then
            if (present(status)) status = excursion_status_success
            return
        end if
        allocate(tmxi(n), tmni(n), tmx(n), tmn(n))
        nmxi = 0
        nmni = 0
        candidate_max = x(1)
        candidate_min = x(1)
        candidate_max_i = 1
        candidate_min_i = 1
        looking_for_maximum = .true.

        do i = 1, n
            if (x(i) > candidate_max) then
                candidate_max = x(i)
                candidate_max_i = i
            end if
            if (x(i) < candidate_min) then
                candidate_min = x(i)
                candidate_min_i = i
            end if
            if (looking_for_maximum) then
                if (x(i) < candidate_max - delta) then
                    nmxi = nmxi + 1
                    tmxi(nmxi) = candidate_max_i
                    tmx(nmxi) = candidate_max
                    candidate_min = x(i)
                    candidate_min_i = i
                    looking_for_maximum = .false.
                end if
            else
                if (x(i) > candidate_min + delta) then
                    nmni = nmni + 1
                    tmni(nmni) = candidate_min_i
                    tmn(nmni) = candidate_min
                    candidate_max = x(i)
                    candidate_max_i = i
                    looking_for_maximum = .true.
                end if
            end if
        end do

        imx = tmxi(:nmxi)
        mx = tmx(:nmxi)
        imn = tmni(:nmni)
        mn = tmn(:nmni)
        if (present(status)) status = excursion_status_success
    end subroutine peak_detect

end module scifort_signal_peak_excursion
