! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Minimum-distance STFT synthesis windows for a prescribed analysis window
! and integer hop. The dual constraints separate into disjoint residue classes
! modulo the hop, so no dense linear algebra is needed.
module scifort_signal_stft_dual
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    implicit none
    private

    integer, parameter, public :: stft_dual_status_success = 0
    integer, parameter, public :: stft_dual_status_invalid_input = 1
    integer, parameter, public :: stft_dual_status_noninvertible = 2

    interface closest_stft_dual_window
        module procedure closest_stft_dual_real
        module procedure closest_stft_dual_complex
    end interface
    public :: closest_stft_dual_window

contains

    subroutine closest_stft_dual_complex(win, hop, dual, alpha, desired_dual, scaled, status)
        complex(dp), intent(in) :: win(:) !! STFT analysis window
        integer, intent(in) :: hop !! analysis hop in samples
        complex(dp), allocatable, intent(out) :: dual(:) !! closest (possibly scaled) dual window
        complex(dp), intent(out), optional :: alpha !! scaling of analysis window
        complex(dp), intent(in), optional :: desired_dual(:) !! target synthesis window; default ones
        logical, intent(in), optional :: scaled !! optimize scaling factor, default true
        integer, intent(out), optional :: status !! success, invalid arguments, or singular overlap
        complex(dp), allocatable :: wanted(:), canonical(:), projected(:)
        complex(dp) :: scaling, inner
        real(dp), allocatable :: energy(:)
        real(dp) :: norm_sq
        integer :: n, h, j
        logical :: do_scale

        allocate(dual(0))
        if (present(alpha)) alpha = (0.0_dp,0.0_dp)
        if (present(status)) status = stft_dual_status_invalid_input
        n = size(win)
        if (n < 1 .or. hop < 1 .or. hop > n) return
        if (.not. all(ieee_is_finite(real(win,dp))) .or. &
            .not. all(ieee_is_finite(aimag(win)))) return
        if (present(desired_dual)) then
            if (size(desired_dual) /= n) return
            if (.not. all(ieee_is_finite(real(desired_dual,dp))) .or. &
                .not. all(ieee_is_finite(aimag(desired_dual)))) return
            wanted = desired_dual
        else
            allocate(wanted(n))
            wanted = (1.0_dp,0.0_dp)
        end if
        allocate(energy(hop))
        energy = 0.0_dp
        do j = 1, n
            h = mod(j-1,hop)+1
            energy(h) = energy(h) + abs(win(j))**2
        end do
        if (.not. all(ieee_is_finite(energy))) return
        if (any(energy <= tiny(1.0_dp))) then
            if (present(status)) status = stft_dual_status_noninvertible
            return
        end if
        allocate(canonical(n), projected(n))
        canonical = win
        do j = 1,n
            canonical(j) = win(j)/energy(mod(j-1,hop)+1)
        end do
        do h = 1,hop
            inner = (0.0_dp,0.0_dp)
            do j = h,n,hop
                inner = inner + conjg(win(j))*wanted(j)
            end do
            do j = h,n,hop
                projected(j) = canonical(j)*inner
            end do
        end do
        do_scale = .true.
        if (present(scaled)) do_scale = scaled
        scaling = (1.0_dp,0.0_dp)
        if (do_scale) then
            norm_sq = sum(abs(projected)**2)
            if (.not. ieee_is_finite(norm_sq)) return
            if (norm_sq <= epsilon(1.0_dp)**2) then
                if (present(status)) status = stft_dual_status_noninvertible
                return
            end if
            scaling = sum(conjg(projected)*canonical)/norm_sq
            if (.not. ieee_is_finite(real(scaling,dp)) .or. &
                .not. ieee_is_finite(aimag(scaling))) return
            if (abs(scaling) <= tiny(1.0_dp)) then
                if (present(status)) status = stft_dual_status_noninvertible
                return
            end if
        end if
        deallocate(dual)
        dual = canonical + scaling*(wanted-projected)
        if (.not. all(ieee_is_finite(real(dual,dp))) .or. &
            .not. all(ieee_is_finite(aimag(dual)))) then
            deallocate(dual)
            allocate(dual(0))
            return
        end if
        if (present(alpha)) alpha = scaling
        if (present(status)) status = stft_dual_status_success
    end subroutine closest_stft_dual_complex

    subroutine closest_stft_dual_real(win, hop, dual, alpha, desired_dual, scaled, status)
        real(dp), intent(in) :: win(:) !! real STFT analysis window
        integer, intent(in) :: hop !! hop size
        real(dp), allocatable, intent(out) :: dual(:) !! real closest synthesis window
        real(dp), intent(out), optional :: alpha !! real scale factor
        real(dp), intent(in), optional :: desired_dual(:) !! real target window
        logical, intent(in), optional :: scaled !! scale optimal dual if true
        integer, intent(out), optional :: status !! zero indicates success
        complex(dp), allocatable :: result(:), target(:)
        complex(dp) :: factor
        integer :: code

        if (present(desired_dual)) then
            target = cmplx(desired_dual,0.0_dp,dp)
            call closest_stft_dual_complex(cmplx(win,0.0_dp,dp), hop, result, &
                factor, desired_dual=target, scaled=scaled, status=code)
        else
            call closest_stft_dual_complex(cmplx(win,0.0_dp,dp), hop, result, &
                factor, scaled=scaled, status=code)
        end if
        dual = real(result,dp)
        if (present(alpha)) alpha = real(factor,dp)
        if (present(status)) status = code
    end subroutine closest_stft_dual_real
end module scifort_signal_stft_dual
