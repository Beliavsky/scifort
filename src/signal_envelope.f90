! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Independent FFT-domain implementation of scipy.signal.envelope.
! Frequency bounds are integer FFT bins; axes are one-based in Fortran.
module scifort_signal_envelope
    use, intrinsic :: iso_fortran_env, only : int64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_signal_spectral, only : forward_fft
    implicit none
    private

    integer, parameter, public :: envelope_status_success = 0
    integer, parameter, public :: envelope_status_invalid_input = 1

    interface envelope
        module procedure envelope_real_1d
        module procedure envelope_complex_1d
        module procedure envelope_real_2d
        module procedure envelope_complex_2d
        module procedure envelope_real_3d
        module procedure envelope_complex_3d
    end interface envelope
    public :: envelope

    interface envelope_into
        module procedure envelope_into_real
        module procedure envelope_into_complex
    end interface envelope_into
    public :: envelope_into

contains

    function envelope_real_1d(z, bp_low, bp_high, n_out, squared, residual, remainder, status) result(y)
        real(dp), intent(in) :: z(:) !! real samples
        integer, intent(in), optional :: bp_low !! inclusive lower FFT bin (default 1)
        integer, intent(in), optional :: bp_high !! exclusive upper FFT bin
        integer, intent(in), optional :: n_out !! resampled length (default input length)
        logical, intent(in), optional :: squared !! return squared envelope
        character(*), intent(in), optional :: residual !! 'lowpass' (default), 'all', or 'none'
        real(dp), allocatable, intent(out), optional :: remainder(:) !! residual signal
        integer, intent(out), optional :: status !! success or invalid input
        real(dp), allocatable :: y(:)
        complex(dp), allocatable :: r(:)
        integer :: state

        call envelope_kernel(cmplx(z,0.0_dp,dp), .true., y, r, state, &
            bp_low, bp_high, n_out, squared, residual)
        if (present(remainder)) remainder = real(r,dp)
        if (present(status)) status = state
    end function envelope_real_1d

    function envelope_complex_1d(z, bp_low, bp_high, n_out, squared, residual, remainder, status) result(y)
        complex(dp), intent(in) :: z(:) !! complex samples
        integer, intent(in), optional :: bp_low !! inclusive lower FFT bin (default 1)
        integer, intent(in), optional :: bp_high !! exclusive upper FFT bin
        integer, intent(in), optional :: n_out !! resampled length (default input length)
        logical, intent(in), optional :: squared !! return squared envelope
        character(*), intent(in), optional :: residual !! 'lowpass' (default), 'all', or 'none'
        complex(dp), allocatable, intent(out), optional :: remainder(:) !! residual signal
        integer, intent(out), optional :: status !! success or invalid input
        real(dp), allocatable :: y(:)
        complex(dp), allocatable :: r(:)
        integer :: state

        call envelope_kernel(z, .false., y, r, state, bp_low, bp_high, n_out, squared, residual)
        if (present(remainder)) remainder = r
        if (present(status)) status = state
    end function envelope_complex_1d

    function envelope_real_2d(z, bp_low, bp_high, n_out, squared, residual, axis, remainder, status) result(y)
        real(dp), intent(in) :: z(:,:) !! rank-2 real samples
        integer, intent(in), optional :: bp_low !! inclusive lower FFT bin
        integer, intent(in), optional :: bp_high !! exclusive upper FFT bin
        integer, intent(in), optional :: n_out !! resampled length along axis
        logical, intent(in), optional :: squared !! return squared envelope
        character(*), intent(in), optional :: residual !! 'lowpass', 'all', or 'none'
        integer, intent(in), optional :: axis !! filtering dimension, 1 or 2 (default 2)
        real(dp), allocatable, intent(out), optional :: remainder(:,:) !! residual samples
        integer, intent(out), optional :: status !! success or invalid input
        real(dp), allocatable :: y(:,:)
        real(dp), allocatable :: line(:)
        complex(dp), allocatable :: r(:)
        integer :: selected, length_out, state, i

        allocate(y(0,0))
        if (present(remainder)) allocate(remainder(0,0))
        if (present(status)) status = envelope_status_invalid_input
        selected = 2
        if (present(axis)) selected = axis
        if (selected < 1 .or. selected > 2) return
        if (.not. envelope_settings_valid(size(z,selected), bp_low, bp_high, n_out, residual)) return
        length_out = size(z,selected)
        if (present(n_out)) length_out = n_out
        if (length_out < 1 .or. size(z,selected) < 1) return
        if (selected == 1) then
            deallocate(y)
            allocate(y(length_out,size(z,2)))
            if (present(remainder)) then
                deallocate(remainder)
                allocate(remainder(length_out,size(z,2)))
            end if
            do i = 1, size(z,2)
                call envelope_kernel(cmplx(z(:,i),0.0_dp,dp), .true., line, r, state, &
                    bp_low, bp_high, n_out, squared, residual)
                if (state /= envelope_status_success) exit
                y(:,i) = line
                if (present(remainder)) then
                    if (size(r) == length_out) remainder(:,i) = real(r,dp)
                end if
            end do
        else
            deallocate(y)
            allocate(y(size(z,1),length_out))
            if (present(remainder)) then
                deallocate(remainder)
                allocate(remainder(size(z,1),length_out))
            end if
            do i = 1, size(z,1)
                call envelope_kernel(cmplx(z(i,:),0.0_dp,dp), .true., line, r, state, &
                    bp_low, bp_high, n_out, squared, residual)
                if (state /= envelope_status_success) exit
                y(i,:) = line
                if (present(remainder)) then
                    if (size(r) == length_out) remainder(i,:) = real(r,dp)
                end if
            end do
        end if
        if (present(status)) status = envelope_status_success
        ! No residual is returned in the 'none' mode even if the caller
        ! provided an optional residual output argument.
        if (present(remainder)) then
            if (present(residual)) then
                if (trim(residual) == 'none') then
                    deallocate(remainder)
                    allocate(remainder(0,0))
                end if
            end if
        end if
        if (size(z,3-selected) == 0) return
        if (state /= envelope_status_success) then
            deallocate(y)
            allocate(y(0,0))
            if (present(remainder)) then
                deallocate(remainder)
                allocate(remainder(0,0))
            end if
            if (present(status)) status = envelope_status_invalid_input
        end if
    end function envelope_real_2d

    function envelope_complex_2d(z, bp_low, bp_high, n_out, squared, residual, axis, remainder, status) result(y)
        complex(dp), intent(in) :: z(:,:) !! rank-2 complex samples
        integer, intent(in), optional :: bp_low !! inclusive lower FFT bin
        integer, intent(in), optional :: bp_high !! exclusive upper FFT bin
        integer, intent(in), optional :: n_out !! resampled length along axis
        logical, intent(in), optional :: squared !! return squared envelope
        character(*), intent(in), optional :: residual !! 'lowpass', 'all', or 'none'
        integer, intent(in), optional :: axis !! filtering dimension, 1 or 2 (default 2)
        complex(dp), allocatable, intent(out), optional :: remainder(:,:) !! residual samples
        integer, intent(out), optional :: status !! success or invalid input
        real(dp), allocatable :: y(:,:)
        real(dp), allocatable :: line(:)
        complex(dp), allocatable :: r(:)
        integer :: selected, length_out, state, i

        allocate(y(0,0))
        if (present(remainder)) allocate(remainder(0,0))
        if (present(status)) status = envelope_status_invalid_input
        selected = 2
        if (present(axis)) selected = axis
        if (selected < 1 .or. selected > 2) return
        if (.not. envelope_settings_valid(size(z,selected), bp_low, bp_high, n_out, residual)) return
        length_out = size(z,selected)
        if (present(n_out)) length_out = n_out
        if (length_out < 1 .or. size(z,selected) < 1) return
        if (selected == 1) then
            deallocate(y)
            allocate(y(length_out,size(z,2)))
            if (present(remainder)) then
                deallocate(remainder)
                allocate(remainder(length_out,size(z,2)))
            end if
            do i = 1, size(z,2)
                call envelope_kernel(z(:,i), .false., line, r, state, &
                    bp_low, bp_high, n_out, squared, residual)
                if (state /= envelope_status_success) exit
                y(:,i) = line
                if (present(remainder)) then
                    if (size(r) == length_out) remainder(:,i) = r
                end if
            end do
        else
            deallocate(y)
            allocate(y(size(z,1),length_out))
            if (present(remainder)) then
                deallocate(remainder)
                allocate(remainder(size(z,1),length_out))
            end if
            do i = 1, size(z,1)
                call envelope_kernel(z(i,:), .false., line, r, state, &
                    bp_low, bp_high, n_out, squared, residual)
                if (state /= envelope_status_success) exit
                y(i,:) = line
                if (present(remainder)) then
                    if (size(r) == length_out) remainder(i,:) = r
                end if
            end do
        end if
        if (present(status)) status = envelope_status_success
        ! No residual is returned in the 'none' mode even if the caller
        ! provided an optional residual output argument.
        if (present(remainder)) then
            if (present(residual)) then
                if (trim(residual) == 'none') then
                    deallocate(remainder)
                    allocate(remainder(0,0))
                end if
            end if
        end if
        if (size(z,3-selected) == 0) return
        if (state /= envelope_status_success) then
            deallocate(y)
            allocate(y(0,0))
            if (present(remainder)) then
                deallocate(remainder)
                allocate(remainder(0,0))
            end if
            if (present(status)) status = envelope_status_invalid_input
        end if
    end function envelope_complex_2d


    ! Common argument checks also apply if the array has zero batch channels.
    logical function envelope_settings_valid(n, bp_low, bp_high, n_out, residual) result(valid)
        integer, intent(in) :: n !! signal length along the filtering axis
        integer, intent(in), optional :: bp_low !! inclusive signed FFT bin
        integer, intent(in), optional :: bp_high !! exclusive signed FFT bin
        integer, intent(in), optional :: n_out !! resampled signal length
        character(*), intent(in), optional :: residual !! lowpass, all, or none
        integer :: lo, hi

        valid = .false.
        if (n < 1) return
        if (present(n_out)) then
            if (n_out < 1) return
        end if
        lo = 1
        if (present(bp_low)) lo = bp_low
        hi = (n+1)/2
        if (present(bp_high)) hi = bp_high
        if (lo < -(n/2) .or. lo >= hi .or. hi > (n+1)/2) return
        if (present(residual)) then
            select case(trim(residual))
            case ('lowpass','all','none')
            case default
                return
            end select
        end if
        valid = .true.
    end function envelope_settings_valid

    function envelope_real_3d(z, bp_low, bp_high, n_out, squared, residual, axis, remainder, status) result(y)
        real(dp), intent(in) :: z(:,:,:) !! rank-three real samples
        integer, intent(in), optional :: bp_low !! inclusive signed FFT bin
        integer, intent(in), optional :: bp_high !! exclusive signed FFT bin
        integer, intent(in), optional :: n_out !! output length along axis
        logical, intent(in), optional :: squared !! return squared envelope
        character(*), intent(in), optional :: residual !! lowpass, all, or none
        integer, intent(in), optional :: axis !! one-based transform axis, default 3
        real(dp), allocatable, intent(out), optional :: remainder(:,:,:) !! residual samples
        integer, intent(out), optional :: status !! success or invalid input
        real(dp), allocatable :: y(:,:,:)
        real(dp), allocatable :: line(:)
        complex(dp), allocatable :: r(:)
        integer :: selected, shp(3), i, j, k, state
        logical :: keep_residual

        allocate(y(0,0,0))
        if (present(remainder)) allocate(remainder(0,0,0))
        if (present(status)) status = envelope_status_invalid_input
        selected = 3
        if (present(axis)) selected = axis
        if (selected < 1 .or. selected > 3) return
        if (.not. envelope_settings_valid(size(z,selected), bp_low, bp_high, n_out, residual)) return
        shp = shape(z)
        if (present(n_out)) shp(selected) = n_out
        keep_residual = .true.
        if (present(residual)) keep_residual = trim(residual) /= 'none'
        deallocate(y)
        allocate(y(shp(1),shp(2),shp(3)))
        if (present(remainder)) then
            if (keep_residual) then
                deallocate(remainder)
                allocate(remainder(shp(1),shp(2),shp(3)))
            end if
        end if
        state = envelope_status_success
        select case(selected)
        case(1)
            do k = 1, size(z,3)
                do j = 1, size(z,2)
                    call envelope_kernel(cmplx(z(:,j,k),0.0_dp,dp), .true., line, r, state, &
                        bp_low, bp_high, n_out, squared, residual)
                    if (state /= envelope_status_success) exit
                    y(:,j,k) = line
                    if (present(remainder)) then
                        if (keep_residual) remainder(:,j,k) = real(r,dp)
                    end if
                end do
                if (state /= envelope_status_success) exit
            end do
        case(2)
            do k = 1, size(z,3)
                do i = 1, size(z,1)
                    call envelope_kernel(cmplx(z(i,:,k),0.0_dp,dp), .true., line, r, state, &
                        bp_low, bp_high, n_out, squared, residual)
                    if (state /= envelope_status_success) exit
                    y(i,:,k) = line
                    if (present(remainder)) then
                        if (keep_residual) remainder(i,:,k) = real(r,dp)
                    end if
                end do
                if (state /= envelope_status_success) exit
            end do
        case(3)
            do j = 1, size(z,2)
                do i = 1, size(z,1)
                    call envelope_kernel(cmplx(z(i,j,:),0.0_dp,dp), .true., line, r, state, &
                        bp_low, bp_high, n_out, squared, residual)
                    if (state /= envelope_status_success) exit
                    y(i,j,:) = line
                    if (present(remainder)) then
                        if (keep_residual) remainder(i,j,:) = real(r,dp)
                    end if
                end do
                if (state /= envelope_status_success) exit
            end do
        end select
        if (state /= envelope_status_success) then
            deallocate(y)
            allocate(y(0,0,0))
            if (present(remainder)) then
                deallocate(remainder)
                allocate(remainder(0,0,0))
            end if
            return
        end if
        if (present(status)) status = envelope_status_success
    end function envelope_real_3d

    function envelope_complex_3d(z, bp_low, bp_high, n_out, squared, residual, axis, remainder, status) result(y)
        complex(dp), intent(in) :: z(:,:,:) !! rank-three complex samples
        integer, intent(in), optional :: bp_low !! inclusive signed FFT bin
        integer, intent(in), optional :: bp_high !! exclusive signed FFT bin
        integer, intent(in), optional :: n_out !! output length along axis
        logical, intent(in), optional :: squared !! return squared envelope
        character(*), intent(in), optional :: residual !! lowpass, all, or none
        integer, intent(in), optional :: axis !! one-based transform axis, default 3
        complex(dp), allocatable, intent(out), optional :: remainder(:,:,:) !! residual samples
        integer, intent(out), optional :: status !! success or invalid input
        real(dp), allocatable :: y(:,:,:)
        real(dp), allocatable :: line(:)
        complex(dp), allocatable :: r(:)
        integer :: selected, shp(3), i, j, k, state
        logical :: keep_residual

        allocate(y(0,0,0))
        if (present(remainder)) allocate(remainder(0,0,0))
        if (present(status)) status = envelope_status_invalid_input
        selected = 3
        if (present(axis)) selected = axis
        if (selected < 1 .or. selected > 3) return
        if (.not. envelope_settings_valid(size(z,selected), bp_low, bp_high, n_out, residual)) return
        shp = shape(z)
        if (present(n_out)) shp(selected) = n_out
        keep_residual = .true.
        if (present(residual)) keep_residual = trim(residual) /= 'none'
        deallocate(y)
        allocate(y(shp(1),shp(2),shp(3)))
        if (present(remainder)) then
            if (keep_residual) then
                deallocate(remainder)
                allocate(remainder(shp(1),shp(2),shp(3)))
            end if
        end if
        state = envelope_status_success
        select case(selected)
        case(1)
            do k = 1, size(z,3)
                do j = 1, size(z,2)
                    call envelope_kernel(z(:,j,k), .false., line, r, state, &
                        bp_low, bp_high, n_out, squared, residual)
                    if (state /= envelope_status_success) exit
                    y(:,j,k) = line
                    if (present(remainder)) then
                        if (keep_residual) remainder(:,j,k) = r
                    end if
                end do
                if (state /= envelope_status_success) exit
            end do
        case(2)
            do k = 1, size(z,3)
                do i = 1, size(z,1)
                    call envelope_kernel(z(i,:,k), .false., line, r, state, &
                        bp_low, bp_high, n_out, squared, residual)
                    if (state /= envelope_status_success) exit
                    y(i,:,k) = line
                    if (present(remainder)) then
                        if (keep_residual) remainder(i,:,k) = r
                    end if
                end do
                if (state /= envelope_status_success) exit
            end do
        case(3)
            do j = 1, size(z,2)
                do i = 1, size(z,1)
                    call envelope_kernel(z(i,j,:), .false., line, r, state, &
                        bp_low, bp_high, n_out, squared, residual)
                    if (state /= envelope_status_success) exit
                    y(i,j,:) = line
                    if (present(remainder)) then
                        if (keep_residual) remainder(i,j,:) = r
                    end if
                end do
                if (state /= envelope_status_success) exit
            end do
        end select
        if (state /= envelope_status_success) then
            deallocate(y)
            allocate(y(0,0,0))
            if (present(remainder)) then
                deallocate(remainder)
                allocate(remainder(0,0,0))
            end if
            return
        end if
        if (present(status)) status = envelope_status_success
    end function envelope_complex_3d

    subroutine envelope_kernel(z, real_input, y, remainder, status, bp_low, bp_high, n_out, squared, residual)
        complex(dp), intent(in) :: z(:) !! samples (real samples supplied with zero imaginary part)
        logical, intent(in) :: real_input !! input was real-valued
        real(dp), allocatable, intent(out) :: y(:) !! signal amplitude envelope
        complex(dp), allocatable, intent(out) :: remainder(:) !! residual samples, if selected
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: bp_low !! inclusive FFT frequency bin
        integer, intent(in), optional :: bp_high !! exclusive FFT frequency bin
        integer, intent(in), optional :: n_out !! output length
        logical, intent(in), optional :: squared !! return squared envelope
        character(*), intent(in), optional :: residual !! 'lowpass', 'all', or 'none'
        complex(dp), allocatable :: spectrum(:), band(:), baseband(:), retained(:)
        integer :: n, nf, lo, hi, f, j, k
        character(:), allocatable :: residual_mode
        logical :: square

        allocate(y(0), remainder(0))
        status = envelope_status_invalid_input
        n = size(z)
        nf = n
        if (present(n_out)) nf = n_out
        if (n < 1 .or. nf < 1) return
        if (.not. all(ieee_is_finite(real(z,dp))) .or. &
            .not. all(ieee_is_finite(aimag(z)))) return
        lo = 1
        if (present(bp_low)) lo = bp_low
        hi = (n+1)/2
        if (present(bp_high)) hi = bp_high
        if (lo < -(n/2) .or. lo >= hi .or. hi > (n+1)/2) return
        residual_mode = 'lowpass'
        if (present(residual)) residual_mode = trim(residual)
        if (residual_mode /= 'lowpass' .and. residual_mode /= 'all' .and. &
            residual_mode /= 'none') return
        square = .false.
        if (present(squared)) square = squared

        spectrum = z
        call forward_fft(spectrum)
        if (real_input) then
            ! Construct the nonnegative-frequency spectrum used by rfft.
            spectrum(n/2+2:n) = cmplx(0.0_dp,0.0_dp,dp)
            if (lo > 0) then
                spectrum(lo+1:hi) = 2.0_dp*spectrum(lo+1:hi)
            else if (hi > 1) then
                spectrum(2:hi) = 2.0_dp*spectrum(2:hi)
            end if
        end if

        ! Shifting a contiguous band to baseband does not change its envelope.
        allocate(band(hi-lo))
        do j = 1, hi-lo
            f = lo+j-1
            k = modulo(f,n)+1
            band(j) = spectrum(k)
        end do
        ! SciPy 1.17 uses a raw NumPy slice Z[lo:0] for a strictly
        ! negative band ending at zero. It is empty, unlike a half-open
        ! signed-frequency interval. Match that compatibility edge case.
        if (hi == 0) band = cmplx(0.0_dp,0.0_dp,dp)
        allocate(baseband(nf))
        baseband = cmplx(0.0_dp,0.0_dp,dp)
        baseband(1:min(nf,size(band))) = band(1:min(nf,size(band)))
        call inverse_fft(baseband)
        baseband = baseband*(real(nf,dp)/real(n,dp))
        deallocate(y)
        allocate(y(nf))
        if (square) then
            y = real(baseband,dp)**2 + aimag(baseband)**2
        else
            y = abs(baseband)
        end if

        if (residual_mode /= 'none') then
            retained = spectrum
            if (lo == 0) then
                ! SciPy Z[0:] clears every bin in its zero-crossing branch.
                retained = cmplx(0.0_dp,0.0_dp,dp)
            else if (residual_mode == 'all') then
                if (hi /= 0) then
                    do f = lo, hi-1
                        retained(modulo(f,n)+1) = cmplx(0.0_dp,0.0_dp,dp)
                    end do
                end if
            else
                ! The signed FFT frequencies below lo comprise the lowpass residual.
                ! Nyquist is represented as the negative frequency -n/2.
                do j = 0, n-1
                    f = j
                    if (j >= (n+1)/2) f = j-n
                    if (f >= lo) retained(j+1) = cmplx(0.0_dp,0.0_dp,dp)
                end do
            end if
            call spectral_resample(retained, nf, real_input, remainder)
        end if
        status = envelope_status_success
    end subroutine envelope_kernel

    subroutine inverse_fft(z)
        complex(dp), intent(inout) :: z(:) !! unscaled FFT coefficients, replaced by IFFT
        z = conjg(z)
        call forward_fft(z)
        z = conjg(z)/real(size(z),dp)
    end subroutine inverse_fft

    subroutine spectral_resample(x, nout, real_input, y)
        complex(dp), intent(in) :: x(:) !! original FFT bins, nonnegative only if real_input
        integer, intent(in) :: nout !! output length
        logical, intent(in) :: real_input !! reconstruct a real-valued residual
        complex(dp), allocatable, intent(out) :: y(:) !! resampled residual
        complex(dp), allocatable :: bins(:)
        integer :: n, m, m2, k, negative_count

        n = size(x)
        m = min(n,nout)
        m2 = m/2+1
        allocate(bins(nout))
        bins = cmplx(0.0_dp,0.0_dp,dp)
        if (real_input) then
            bins(:m2) = x(:m2)
            if (mod(m,2) == 0 .and. nout /= n) then
                if (nout < n) then
                    bins(m2) = bins(m2)*2.0_dp
                else
                    bins(m2) = bins(m2)*0.5_dp
                end if
            end if
            do k = 1, (nout-1)/2
                if (k < m2) bins(nout-k+1) = conjg(bins(k+1))
            end do
        else
            bins(:m2) = x(:m2)
            negative_count = m-m2
            if (negative_count > 0) bins(nout-negative_count+1:nout) = x(n-negative_count+1:n)
            if (mod(m,2) == 0 .and. nout /= n) then
                if (nout < n) then
                    bins(nout-m/2+1) = bins(nout-m/2+1) + x(n-m/2+1)
                else
                    bins(m/2+1) = bins(m/2+1)*0.5_dp
                    bins(nout-m/2+1) = bins(m/2+1)
                end if
            end if
        end if
        call inverse_fft(bins)
        y = bins*(real(nout,dp)/real(n,dp))
        if (real_input) y = cmplx(real(y,dp),0.0_dp,dp)
    end subroutine spectral_resample


    ! Caller-output envelope analysis across ranks 1..15. This SciFort
    ! interface complements SciPy's array-returning result and returns the
    ! optional residual through a separate argument of the input's kind.

    subroutine envelope_into_real(z, y, bp_low, bp_high, n_out, squared, residual, axis, remainder, status)
        real(dp), intent(in) :: z(..) !! input signal, rank 1..15
        real(dp), intent(inout) :: y(..) !! caller-owned envelope, input shape with modified time axis
        integer, intent(in), optional :: bp_low !! inclusive signed FFT bin, default 1
        integer, intent(in), optional :: bp_high !! exclusive signed FFT bin, default (n+1)/2
        integer, intent(in), optional :: n_out !! length of transformed axis, default n
        logical, intent(in), optional :: squared !! return squared envelope if true
        character(*), intent(in), optional :: residual !! 'lowpass' (default), 'all', or 'none'
        integer, intent(in), optional :: axis !! one-based axis or negative from end, default last
        real(dp), intent(inout), optional :: remainder(..) !! optional residual, same shape as y
        integer, intent(out), optional :: status !! success or invalid input
        complex(dp), allocatable :: samples(:), residual_flat(:)
        real(dp), allocatable :: envelope_flat(:)
        integer, allocatable :: remainder_shape(:)
        integer :: selected, state, remainder_rank
        logical :: want_residual, valid

        if (present(status)) status = envelope_status_invalid_input
        remainder_rank = -1
        allocate(remainder_shape(0))
        if (present(remainder)) then
            remainder_shape = shape(remainder)
            remainder_rank = rank(remainder)
        end if
        call envelope_nd_check(shape(z), shape(y), remainder_shape, remainder_rank, &
            bp_low, bp_high, n_out, residual, axis, selected, want_residual, valid)
        if (.not. valid) return
        call envelope_nd_flatten_real(z, samples)
        call envelope_nd_lanes(samples, .true., shape(z), selected, &
            size(y, kind=int64), envelope_flat, residual_flat, state, want_residual, &
            bp_low, bp_high, n_out, squared, residual)
        if (state /= envelope_status_success) return
        call envelope_nd_scatter_real(y, envelope_flat)
        if (want_residual) then
            call envelope_nd_scatter_real(remainder, real(residual_flat,dp))
        end if
        if (present(status)) status = envelope_status_success
    end subroutine envelope_into_real

    subroutine envelope_into_complex(z, y, bp_low, bp_high, n_out, squared, residual, axis, remainder, status)
        complex(dp), intent(in) :: z(..) !! input signal, rank 1..15
        real(dp), intent(inout) :: y(..) !! caller-owned envelope, input shape with modified time axis
        integer, intent(in), optional :: bp_low !! inclusive signed FFT bin, default 1
        integer, intent(in), optional :: bp_high !! exclusive signed FFT bin, default (n+1)/2
        integer, intent(in), optional :: n_out !! length of transformed axis, default n
        logical, intent(in), optional :: squared !! return squared envelope if true
        character(*), intent(in), optional :: residual !! 'lowpass' (default), 'all', or 'none'
        integer, intent(in), optional :: axis !! one-based axis or negative from end, default last
        complex(dp), intent(inout), optional :: remainder(..) !! optional residual, same shape as y
        integer, intent(out), optional :: status !! success or invalid input
        complex(dp), allocatable :: samples(:), residual_flat(:)
        real(dp), allocatable :: envelope_flat(:)
        integer, allocatable :: remainder_shape(:)
        integer :: selected, state, remainder_rank
        logical :: want_residual, valid

        if (present(status)) status = envelope_status_invalid_input
        remainder_rank = -1
        allocate(remainder_shape(0))
        if (present(remainder)) then
            remainder_shape = shape(remainder)
            remainder_rank = rank(remainder)
        end if
        call envelope_nd_check(shape(z), shape(y), remainder_shape, remainder_rank, &
            bp_low, bp_high, n_out, residual, axis, selected, want_residual, valid)
        if (.not. valid) return
        call envelope_nd_flatten_complex(z, samples)
        call envelope_nd_lanes(samples, .false., shape(z), selected, &
            size(y, kind=int64), envelope_flat, residual_flat, state, want_residual, &
            bp_low, bp_high, n_out, squared, residual)
        if (state /= envelope_status_success) return
        call envelope_nd_scatter_real(y, envelope_flat)
        if (want_residual) then
            call envelope_nd_scatter_complex(remainder, residual_flat)
        end if
        if (present(status)) status = envelope_status_success
    end subroutine envelope_into_complex

    subroutine envelope_nd_check(dims, y_dims, r_dims, r_rank, bp_low, bp_high, &
        n_out, residual, axis, selected, want_residual, valid)
        integer, intent(in) :: dims(:) !! input shape, ranks 1..15
        integer, intent(in) :: y_dims(:) !! caller envelope output shape
        integer, intent(in) :: r_dims(:) !! caller residual output shape, if supplied
        integer, intent(in) :: r_rank !! residual rank, -1 when not supplied
        integer, intent(in), optional :: bp_low !! lower signed frequency index
        integer, intent(in), optional :: bp_high !! upper signed frequency index
        integer, intent(in), optional :: n_out !! resampled signal length
        character(*), intent(in), optional :: residual !! residual mode
        integer, intent(in), optional :: axis !! one-based or negative time axis
        integer, intent(out) :: selected !! normalized positive time axis
        logical, intent(out) :: want_residual !! whether caller requests a residual
        logical, intent(out) :: valid !! all dimensions and parameters are valid
        integer, allocatable :: expected(:)
        character(:), allocatable :: mode

        valid = .false.
        want_residual = .false.
        selected = size(dims)
        if (size(dims) < 1 .or. size(dims) > 15) return
        if (size(y_dims) /= size(dims)) return
        if (present(axis)) selected = axis
        if (selected < 0) selected = selected + size(dims) + 1
        if (selected < 1 .or. selected > size(dims)) return
        if (.not. envelope_settings_valid(dims(selected), bp_low, bp_high, n_out, residual)) return
        expected = dims
        if (present(n_out)) expected(selected) = n_out
        if (any(y_dims /= expected)) return
        mode = 'lowpass'
        if (present(residual)) mode = trim(residual)
        want_residual = r_rank >= 0 .and. mode /= 'none'
        if (want_residual) then
            if (r_rank /= size(dims)) return
            if (size(r_dims) /= size(dims)) return
            if (any(r_dims /= expected)) return
        end if
        valid = .true.
    end subroutine envelope_nd_check

    subroutine envelope_nd_lanes(src, real_input, dims, axis, output_size, &
        y, r, status, want_residual, bp_low, bp_high, n_out, squared, residual)
        use, intrinsic :: iso_fortran_env, only : int64
        complex(dp), intent(in) :: src(:) !! contiguous Fortran-order signal samples
        logical, intent(in) :: real_input !! real input uses one-sided FFT conventions
        integer, intent(in) :: dims(:) !! input shape
        integer, intent(in) :: axis !! normalized one-based sample dimension
        integer(int64), intent(in) :: output_size !! number of output elements
        real(dp), allocatable, intent(out) :: y(:) !! fully staged envelope output
        complex(dp), allocatable, intent(out) :: r(:) !! fully staged residual output
        integer, intent(out) :: status !! success or invalid input
        logical, intent(in) :: want_residual !! compute residual output if true
        integer, intent(in), optional :: bp_low !! inclusive signed FFT bin
        integer, intent(in), optional :: bp_high !! exclusive signed FFT bin
        integer, intent(in), optional :: n_out !! output sample count
        logical, intent(in), optional :: squared !! return squared magnitude
        character(*), intent(in), optional :: residual !! requested residual mode
        real(dp), allocatable :: line_env(:)
        complex(dp), allocatable :: line(:), line_res(:)
        character(:), allocatable :: mode
        integer(int64) :: inner, outer, nx, ny, batch, channel, j, ii, io

        status = envelope_status_invalid_input
        nx = int(dims(axis),int64)
        ny = nx
        if (present(n_out)) ny = int(n_out,int64)
        inner = product(int(dims(:axis-1),int64))
        outer = product(int(dims(axis+1:),int64))
        if (nx < 1 .or. ny < 1 .or. inner < 0 .or. outer < 0) return
        if (size(src,kind=int64) /= inner*nx*outer) return
        if (output_size /= inner*ny*outer) return
        allocate(y(output_size), r(0))
        if (want_residual) then
            deallocate(r)
            allocate(r(output_size))
        end if
        if (output_size == 0_int64) then
            status = envelope_status_success
            return
        end if
        mode = 'none'
        if (want_residual) then
            mode = 'lowpass'
            if (present(residual)) mode = trim(residual)
        end if
        allocate(line(nx))
        do batch = 0_int64, outer-1_int64
            do channel = 0_int64, inner-1_int64
                do j = 0_int64, nx-1_int64
                    ii = batch*inner*nx + j*inner + channel + 1_int64
                    line(j+1) = src(ii)
                end do
                call envelope_kernel(line, real_input, line_env, line_res, status, &
                    bp_low, bp_high, n_out, squared, mode)
                if (status /= envelope_status_success) return
                if (size(line_env,kind=int64) /= ny) then
                    status = envelope_status_invalid_input
                    return
                end if
                do j = 0_int64, ny-1_int64
                    io = batch*inner*ny + j*inner + channel + 1_int64
                    y(io) = line_env(j+1)
                    if (want_residual) r(io) = line_res(j+1)
                end do
            end do
        end do
        status = envelope_status_success
    end subroutine envelope_nd_lanes

    subroutine envelope_nd_flatten_real(z, flat)
        real(dp), intent(in) :: z(..) !! possibly noncontiguous input, rank 1..15
        complex(dp), allocatable, intent(out) :: flat(:) !! contiguous Fortran-order copy
        select rank (z)
        rank (1)
            flat = cmplx(reshape(z,[size(z)]),0.0_dp,dp)
        rank (2)
            flat = cmplx(reshape(z,[size(z)]),0.0_dp,dp)
        rank (3)
            flat = cmplx(reshape(z,[size(z)]),0.0_dp,dp)
        rank (4)
            flat = cmplx(reshape(z,[size(z)]),0.0_dp,dp)
        rank (5)
            flat = cmplx(reshape(z,[size(z)]),0.0_dp,dp)
        rank (6)
            flat = cmplx(reshape(z,[size(z)]),0.0_dp,dp)
        rank (7)
            flat = cmplx(reshape(z,[size(z)]),0.0_dp,dp)
        rank (8)
            flat = cmplx(reshape(z,[size(z)]),0.0_dp,dp)
        rank (9)
            flat = cmplx(reshape(z,[size(z)]),0.0_dp,dp)
        rank (10)
            flat = cmplx(reshape(z,[size(z)]),0.0_dp,dp)
        rank (11)
            flat = cmplx(reshape(z,[size(z)]),0.0_dp,dp)
        rank (12)
            flat = cmplx(reshape(z,[size(z)]),0.0_dp,dp)
        rank (13)
            flat = cmplx(reshape(z,[size(z)]),0.0_dp,dp)
        rank (14)
            flat = cmplx(reshape(z,[size(z)]),0.0_dp,dp)
        rank (15)
            flat = cmplx(reshape(z,[size(z)]),0.0_dp,dp)
        rank default
            allocate(flat(0))
        end select
    end subroutine envelope_nd_flatten_real

    subroutine envelope_nd_flatten_complex(z, flat)
        complex(dp), intent(in) :: z(..) !! possibly noncontiguous input, rank 1..15
        complex(dp), allocatable, intent(out) :: flat(:) !! contiguous Fortran-order copy
        select rank (z)
        rank (1)
            flat = reshape(z,[size(z)])
        rank (2)
            flat = reshape(z,[size(z)])
        rank (3)
            flat = reshape(z,[size(z)])
        rank (4)
            flat = reshape(z,[size(z)])
        rank (5)
            flat = reshape(z,[size(z)])
        rank (6)
            flat = reshape(z,[size(z)])
        rank (7)
            flat = reshape(z,[size(z)])
        rank (8)
            flat = reshape(z,[size(z)])
        rank (9)
            flat = reshape(z,[size(z)])
        rank (10)
            flat = reshape(z,[size(z)])
        rank (11)
            flat = reshape(z,[size(z)])
        rank (12)
            flat = reshape(z,[size(z)])
        rank (13)
            flat = reshape(z,[size(z)])
        rank (14)
            flat = reshape(z,[size(z)])
        rank (15)
            flat = reshape(z,[size(z)])
        rank default
            allocate(flat(0))
        end select
    end subroutine envelope_nd_flatten_complex

    subroutine envelope_nd_scatter_real(out, flat)
        real(dp), intent(inout) :: out(..) !! caller-owned array, possibly strided
        real(dp), intent(in) :: flat(:) !! complete Fortran-order output
        select rank (out)
        rank (1)
            out = reshape(flat, shape(out))
        rank (2)
            out = reshape(flat, shape(out))
        rank (3)
            out = reshape(flat, shape(out))
        rank (4)
            out = reshape(flat, shape(out))
        rank (5)
            out = reshape(flat, shape(out))
        rank (6)
            out = reshape(flat, shape(out))
        rank (7)
            out = reshape(flat, shape(out))
        rank (8)
            out = reshape(flat, shape(out))
        rank (9)
            out = reshape(flat, shape(out))
        rank (10)
            out = reshape(flat, shape(out))
        rank (11)
            out = reshape(flat, shape(out))
        rank (12)
            out = reshape(flat, shape(out))
        rank (13)
            out = reshape(flat, shape(out))
        rank (14)
            out = reshape(flat, shape(out))
        rank (15)
            out = reshape(flat, shape(out))
        rank default
            return
        end select
    end subroutine envelope_nd_scatter_real

    subroutine envelope_nd_scatter_complex(out, flat)
        complex(dp), intent(inout) :: out(..) !! caller-owned array, possibly strided
        complex(dp), intent(in) :: flat(:) !! complete Fortran-order output
        select rank (out)
        rank (1)
            out = reshape(flat, shape(out))
        rank (2)
            out = reshape(flat, shape(out))
        rank (3)
            out = reshape(flat, shape(out))
        rank (4)
            out = reshape(flat, shape(out))
        rank (5)
            out = reshape(flat, shape(out))
        rank (6)
            out = reshape(flat, shape(out))
        rank (7)
            out = reshape(flat, shape(out))
        rank (8)
            out = reshape(flat, shape(out))
        rank (9)
            out = reshape(flat, shape(out))
        rank (10)
            out = reshape(flat, shape(out))
        rank (11)
            out = reshape(flat, shape(out))
        rank (12)
            out = reshape(flat, shape(out))
        rank (13)
            out = reshape(flat, shape(out))
        rank (14)
            out = reshape(flat, shape(out))
        rank (15)
            out = reshape(flat, shape(out))
        rank default
            return
        end select
    end subroutine envelope_nd_scatter_complex

end module scifort_signal_envelope
