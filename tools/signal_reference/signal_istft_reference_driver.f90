! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
! Binary stream driver for independent scipy.signal.istft comparisons.
program signal_istft_reference_driver
    use, intrinsic :: iso_fortran_env, only : int32
    use scifort_kinds, only : dp
    use scifort_signal_istft, only : istft
    implicit none
    integer(int32) :: header(10), out_header(2)
    integer :: in_unit, out_unit, ios, nw, i, status
    character(len=512) :: inpath, outpath
    real(dp) :: fs
    real(dp), allocatable :: real_samples(:), w(:)
    complex(dp), allocatable :: complex_samples(:), spectra(:,:)
    complex(dp), allocatable :: output(:)
    character(len=16) :: window_name, scaling
    logical :: boundary

    if (command_argument_count() /= 2) error stop 'Usage: driver input.bin output.bin'
    call get_command_argument(1, inpath)
    call get_command_argument(2, outpath)
    open(newunit=in_unit, file=trim(inpath), status='old', action='read', access='stream', form='unformatted')
    open(newunit=out_unit, file=trim(outpath), status='replace', action='write', access='stream', form='unformatted')
    do
        read(in_unit, iostat=ios) header
        if (ios /= 0) exit
        nw = header(3)
        read(in_unit) fs
        allocate(spectra(header(1), header(2)))
        read(in_unit) spectra
        boundary = (header(8) /= 0)
        scaling = 'spectrum'
        if (header(9) == 1) scaling = 'psd'
        window_name = 'hann'
        select case (header(10))
        case (1)
            window_name = 'hamming'
        case (2)
            window_name = 'boxcar'
        case (3)
            window_name = 'blackman'
        case (4)
            allocate(w(nw))
            do i = 1, nw
                w(i) = 0.5_dp - 0.5_dp*cos(2.0_dp*acos(-1.0_dp)*real(i-1,dp)/real(nw-1,dp))
            end do
        case (5)
            allocate(w(nw))
            w = 1.0_dp
            w(1:nw:2) = 0.0_dp
        end select
        if (header(7) /= 0) then
            if (allocated(w)) then
                call istft(spectra, complex_samples, status, fs=fs, nperseg=nw, &
                    noverlap=int(header(4)), nfft=int(header(5)), boundary=boundary, &
                    scaling=trim(scaling), window_values=w)
            else
                call istft(spectra, complex_samples, status, fs=fs, nperseg=nw, &
                    noverlap=int(header(4)), nfft=int(header(5)), boundary=boundary, &
                    scaling=trim(scaling), window=trim(window_name))
            end if
            output = complex_samples
            deallocate(complex_samples)
        else
            if (allocated(w)) then
                call istft(spectra, real_samples, status, fs=fs, nperseg=nw, &
                    noverlap=int(header(4)), nfft=int(header(5)), boundary=boundary, &
                    scaling=trim(scaling), window_values=w)
            else
                call istft(spectra, real_samples, status, fs=fs, nperseg=nw, &
                    noverlap=int(header(4)), nfft=int(header(5)), boundary=boundary, &
                    scaling=trim(scaling), window=trim(window_name))
            end if
            output = cmplx(real_samples, 0.0_dp, dp)
            deallocate(real_samples)
        end if
        out_header = [int(status,int32),int(size(output),int32)]
        write(out_unit) out_header
        write(out_unit) output
        deallocate(spectra, output)
        if (allocated(w)) deallocate(w)
    end do
    close(in_unit)
    close(out_unit)
end program signal_istft_reference_driver
