! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Rank-independent inverse of the legacy scipy.signal.stft representation.
! The time and frequency axes can be in any two distinct positions.
! The selected frequency axis becomes samples when time is removed.
! No destination is modified unless every channel reconstructs successfully.
module scifort_signal_istft_nd
    use scifort_kinds, only : dp
    use scifort_signal_istft, only : istft, istft_status_success, &
        istft_status_invalid_input
    use scifort_signal_spectral_nd, only : flatten_complex, scatter_real, scatter_complex
    implicit none
    private
    public :: istft_into

    interface istft_into
        module procedure istft_into_real, istft_into_complex
    end interface

contains

    subroutine validate_inverse_shape(d, dout, axis, time_axis, freq_axis, a, fa, ta, status)
        integer, intent(in) :: d(:) !! STFT input shape (rank 2 to 15)
        integer, intent(in) :: dout(:) !! reconstructed signal shape (rank 1 to 14)
        integer, intent(in), optional :: axis !! one-based/negative output sample axis
        integer, intent(in), optional :: time_axis !! one-based/negative input frame axis
        integer, intent(in), optional :: freq_axis !! one-based/negative input frequency axis
        integer, intent(out) :: a !! resolved output sample axis
        integer, intent(out) :: fa !! resolved input frequency axis
        integer, intent(out) :: ta !! resolved input frame axis
        integer, intent(out) :: status !! zero if dimensions and axes are valid
        integer :: nd, j, mapped

        status = istft_status_invalid_input
        nd = size(dout)
        a = nd
        fa = nd
        ta = nd+1
        if (nd < 1 .or. nd > 14 .or. size(d) /= nd+1) return
        if (any(d <= 0) .or. any(dout <= 0)) return

        if (present(time_axis)) ta = time_axis
        if (ta < 0) ta = nd+2+ta
        if (ta < 1 .or. ta > nd+1) return

        ! If freq_axis is given, it identifies the output sample dimension
        ! after the time axis has been removed. The legacy axis argument
        ! can additionally specify that same output dimension.
        if (present(freq_axis)) then
            fa = freq_axis
            if (fa < 0) fa = nd+2+fa
            if (fa < 1 .or. fa > nd+1 .or. fa == ta) return
            a = fa
            if (fa > ta) a = fa-1
        end if
        if (present(axis)) then
            a = axis
            if (a < 0) a = nd+1+a
        end if
        if (a < 1 .or. a > nd) return

        ! Map the output dimensions to their input positions, skipping time.
        mapped = a
        if (mapped >= ta) mapped = mapped+1
        if (.not. present(freq_axis)) fa = mapped
        if (fa /= mapped) return
        do j = 1,nd
            mapped = j
            if (j >= ta) mapped = j+1
            if (j /= a .and. d(mapped) /= dout(j)) return
        end do
        status = istft_status_success
    end subroutine validate_inverse_shape

    subroutine process_inverse(x, d, dout, a, fa, ta, real_output, staged, times, status, &
        fs, window, nperseg, noverlap, nfft, boundary, scaling, window_values, input_onesided, &
        window_parameter, window_parameter2, window_coefficients)
        complex(dp), intent(in) :: x(:) !! flattened frequency/time coefficients
        integer, intent(in) :: d(:) !! STFT dimensions
        integer, intent(in) :: dout(:) !! reconstructed signal dimensions
        integer, intent(in) :: a !! output sample dimension
        integer, intent(in) :: fa !! input frequency dimension
        integer, intent(in) :: ta !! input frame dimension
        logical, intent(in) :: real_output !! reconstruct real rather than complex samples
        complex(dp), intent(out) :: staged(:) !! private output buffer
        real(dp), allocatable, intent(out) :: times(:) !! sample time coordinates
        integer, intent(out) :: status !! success only if all channels are invertible
        real(dp), intent(in), optional :: fs !! samples per unit time, default 1
        character(*), intent(in), optional :: window !! periodic window, default hann
        integer, intent(in), optional :: nperseg !! segment size, inferred if omitted
        integer, intent(in), optional :: noverlap !! overlap count, default half a segment
        integer, intent(in), optional :: nfft !! FFT length, inferred if omitted
        logical, intent(in), optional :: boundary !! discard extensions, default true
        character(*), intent(in), optional :: scaling !! spectrum (default) or psd
        real(dp), intent(in), optional :: window_values(:) !! explicit analysis window
        logical, intent(in), optional :: input_onesided !! default true for real output
        real(dp), intent(in), optional :: window_parameter !! first named window shape parameter
        real(dp), intent(in), optional :: window_parameter2 !! second named window parameter
        real(dp), intent(in), optional :: window_coefficients(:) !! general cosine weights

        complex(dp), allocatable :: lane(:,:), reconstructed(:)
        real(dp), allocatable :: real_samples(:), lane_times(:)
        integer :: inner, outer, nf, nt, nout, i, j, k, b, q, index, st
        integer :: base, tmp, mapped, nd
        integer, allocatable :: stride(:), output_to_input(:)

        status = istft_status_invalid_input
        allocate(times(0))
        inner = product(dout(:a-1))
        outer = product(dout(a+1:))
        nf = d(fa)
        nt = d(ta)
        nout = dout(a)
        nd = size(dout)
        allocate(stride(nd+1), output_to_input(nd))
        stride(1) = 1
        do j = 2,nd+1
            stride(j) = stride(j-1)*d(j-1)
        end do
        do j = 1,nd
            mapped = j
            if (j >= ta) mapped = j+1
            output_to_input(j) = mapped
        end do
        allocate(lane(nf,nt))
        do k = 0,outer-1
            do i = 0,inner-1
                ! Determine the start of this channel in the flattened input.
                ! The two transform axes may be interleaved with batch axes.
                base = 0
                tmp = i
                do j = 1,a-1
                    base = base + mod(tmp,dout(j))*stride(output_to_input(j))
                    tmp = tmp / dout(j)
                end do
                tmp = k
                do j = a+1,nd
                    base = base + mod(tmp,dout(j))*stride(output_to_input(j))
                    tmp = tmp / dout(j)
                end do
                do b = 0,nt-1
                    do j = 0,nf-1
                        index = base + j*stride(fa) + b*stride(ta) + 1
                        lane(j+1,b+1) = x(index)
                    end do
                end do
                if (real_output) then
                    call istft(lane, real_samples, st, fs=fs, window=window, &
                        nperseg=nperseg, noverlap=noverlap, nfft=nfft, &
                        boundary=boundary, scaling=scaling, times=lane_times, &
                        window_values=window_values, input_onesided=input_onesided, &
                        window_parameter=window_parameter, window_parameter2=window_parameter2, &
                        window_coefficients=window_coefficients)
                    if (st /= istft_status_success) then
                        status = st
                        return
                    end if
                    if (size(real_samples) /= nout) return
                else
                    call istft(lane, reconstructed, st, fs=fs, window=window, &
                        nperseg=nperseg, noverlap=noverlap, nfft=nfft, &
                        boundary=boundary, scaling=scaling, times=lane_times, &
                        window_values=window_values, input_onesided=input_onesided, &
                        window_parameter=window_parameter, window_parameter2=window_parameter2, &
                        window_coefficients=window_coefficients)
                    if (st /= istft_status_success) then
                        status = st
                        return
                    end if
                    if (size(reconstructed) /= nout) return
                end if
                if (size(lane_times) /= nout) return
                if (k == 0 .and. i == 0) then
                    times = lane_times
                else
                    if (any(lane_times /= times)) return
                end if
                do q = 0,nout-1
                    index = k*inner*nout+q*inner+i+1
                    if (real_output) then
                        staged(index) = cmplx(real_samples(q+1),0.0_dp,dp)
                    else
                        staged(index) = reconstructed(q+1)
                    end if
                end do
            end do
        end do
        status = istft_status_success
    end subroutine process_inverse

    subroutine istft_into_real(zxx, result, times, fs, window, nperseg, noverlap, nfft, &
        boundary, scaling, axis, status, window_values, input_onesided, time_axis, freq_axis, &
        window_parameter, window_parameter2, window_coefficients)
        complex(dp), intent(in) :: zxx(..) !! input rank 2..15, configurable transform axes
        real(dp), intent(inout) :: result(..) !! caller-sized real signal, rank 1..14
        real(dp), allocatable, intent(out) :: times(:) !! sample times, empty on error
        real(dp), intent(in), optional :: fs !! sampling frequency, default 1
        character(*), intent(in), optional :: window !! periodic window, default hann
        integer, intent(in), optional :: nperseg !! samples per segment
        integer, intent(in), optional :: noverlap !! overlap samples
        integer, intent(in), optional :: nfft !! FFT size
        logical, intent(in), optional :: boundary !! strip extensions, default true
        character(*), intent(in), optional :: scaling !! spectrum or psd
        integer, intent(in), optional :: axis !! one-based/negative output sample dimension
        integer, intent(out), optional :: status !! zero if every channel succeeds
        real(dp), intent(in), optional :: window_values(:) !! explicit synthesis window
        logical, intent(in), optional :: input_onesided !! must be true, default true
        real(dp), intent(in), optional :: window_parameter !! first named window shape parameter
        real(dp), intent(in), optional :: window_parameter2 !! second named window parameter
        real(dp), intent(in), optional :: window_coefficients(:) !! general cosine weights
        integer, intent(in), optional :: time_axis !! one-based/negative input frame axis
        integer, intent(in), optional :: freq_axis !! one-based/negative input frequency axis
        complex(dp), allocatable :: flat(:), staged(:)
        integer :: a, fa, ta, st

        call validate_inverse_shape(shape(zxx),shape(result),axis,time_axis,freq_axis,a,fa,ta,st)
        if (st == istft_status_success) then
            call flatten_complex(zxx,flat)
            allocate(staged(size(result)))
            call process_inverse(flat,shape(zxx),shape(result),a,fa,ta,.true.,staged,times,st, &
                fs,window,nperseg,noverlap,nfft,boundary,scaling,window_values,input_onesided,window_parameter, &
                window_parameter2,window_coefficients)
            if (st == istft_status_success) call scatter_real(result,staged)
        else
            allocate(times(0))
        end if
        if (st /= istft_status_success) then
            if (allocated(times)) deallocate(times)
            allocate(times(0))
        end if
        if (present(status)) status = st
    end subroutine istft_into_real

    subroutine istft_into_complex(zxx, result, times, fs, window, nperseg, noverlap, nfft, &
        boundary, scaling, axis, status, window_values, input_onesided, time_axis, freq_axis, &
        window_parameter, window_parameter2, window_coefficients)
        complex(dp), intent(in) :: zxx(..) !! input rank 2..15, configurable transform axes
        complex(dp), intent(inout) :: result(..) !! caller-sized complex signal, rank 1..14
        real(dp), allocatable, intent(out) :: times(:) !! sample times, empty on error
        real(dp), intent(in), optional :: fs !! sampling frequency, default 1
        character(*), intent(in), optional :: window !! periodic window, default hann
        integer, intent(in), optional :: nperseg !! samples per segment
        integer, intent(in), optional :: noverlap !! overlap samples
        integer, intent(in), optional :: nfft !! FFT size
        logical, intent(in), optional :: boundary !! strip extensions, default true
        character(*), intent(in), optional :: scaling !! spectrum or psd
        integer, intent(in), optional :: axis !! one-based/negative output sample dimension
        integer, intent(out), optional :: status !! zero if every channel succeeds
        real(dp), intent(in), optional :: window_values(:) !! explicit synthesis window
        logical, intent(in), optional :: input_onesided !! default false for complex output
        real(dp), intent(in), optional :: window_parameter !! first named window shape parameter
        real(dp), intent(in), optional :: window_parameter2 !! second named window parameter
        real(dp), intent(in), optional :: window_coefficients(:) !! general cosine weights
        integer, intent(in), optional :: time_axis !! one-based/negative input frame axis
        integer, intent(in), optional :: freq_axis !! one-based/negative input frequency axis
        complex(dp), allocatable :: flat(:), staged(:)
        integer :: a, fa, ta, st

        call validate_inverse_shape(shape(zxx),shape(result),axis,time_axis,freq_axis,a,fa,ta,st)
        if (st == istft_status_success) then
            call flatten_complex(zxx,flat)
            allocate(staged(size(result)))
            call process_inverse(flat,shape(zxx),shape(result),a,fa,ta,.false.,staged,times,st, &
                fs,window,nperseg,noverlap,nfft,boundary,scaling,window_values,input_onesided,window_parameter, &
                window_parameter2,window_coefficients)
            if (st == istft_status_success) call scatter_complex(result,staged)
        else
            allocate(times(0))
        end if
        if (st /= istft_status_success) then
            if (allocated(times)) deallocate(times)
            allocate(times(0))
        end if
        if (present(status)) status = st
    end subroutine istft_into_complex

end module scifort_signal_istft_nd
