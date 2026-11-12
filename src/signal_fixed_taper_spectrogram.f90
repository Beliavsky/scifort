! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Direct finite-window multichannel spectral matrices for a deterministic
! bank of tapers. Point frequencies and continuous frequency-band integrals
! use the same finite Fourier transform definition as the exact stationary
! VAR(p) population targets; no sampled FFT grid is used for the bands.
!
! SciFort extension, not a scipy.signal public API. No data-adaptive weights.
module scifort_signal_fixed_taper_spectrogram
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_signal_spectral, only : spectral_status_success
    use scifort_signal_multitaper, only : multitaper_spectral_matrix_spectrogram_result
    implicit none
    private

    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: twopi = 2.0_dp*pi

    interface fixed_taper_spectral_matrix_spectrogram
        module procedure fixed_taper_real
        module procedure fixed_taper_complex
    end interface
    public :: fixed_taper_spectral_matrix_spectrogram

contains

    pure real(dp) function sinc_pi(x) result(y)
        real(dp), intent(in) :: x !! real argument in units of pi
        if (abs(x) < 1.0e-5_dp) then
            y = 1.0_dp-(pi*x)**2/6.0_dp+(pi*x)**4/120.0_dp
        else
            y = sin(pi*x)/(pi*x)
        end if
    end function sinc_pi

    function fixed_taper_core(x,tapers,frequencies,bands,taper_weights, &
            normalized,average,noverlap) result(out)
        complex(dp), intent(in) :: x(:,:) !! finite samples by channels, time in first dimension
        real(dp), intent(in) :: tapers(:,:) !! finite deterministic windows (length,taper)
        real(dp), intent(in) :: frequencies(:) !! point frequencies in cycles per sample
        real(dp), intent(in), optional :: bands(:,:) !! band edges (2,n), cycles per sample
        real(dp), intent(in), optional :: taper_weights(:) !! fixed nonnegative weights, default equal
        logical, intent(in), optional :: normalized !! divide EACH taper by squared energy, default true
        logical, intent(in), optional :: average !! divide band integrals by width, default true
        integer, intent(in), optional :: noverlap !! shared samples between windows, default zero
        type(multitaper_spectral_matrix_spectrogram_result) :: out
        real(dp), allocatable :: weights(:), energies(:)
        complex(dp), allocatable :: samples(:,:), kernel(:,:), transform(:), phase(:)
        complex(dp), allocatable :: work(:,:)
        real(dp) :: largest, total, factor, width, midpoint, lag
        integer :: length, kmax, nc, nf, nb, hop, nw, i, j, k, t, s, f, w, start
        logical :: norm, avg

        length = size(tapers,1)
        kmax = size(tapers,2)
        nc = size(x,2)
        nf = size(frequencies)
        nb = 0
        if (present(bands)) then
            if (size(bands,1) /= 2) return
            nb = size(bands,2)
        end if
        if (length < 1 .or. kmax < 1 .or. nc < 1) return
        if (nf+nb < 1 .or. size(x,1) < length) return
        if (.not. all(ieee_is_finite(real(x,dp))) .or. &
            .not. all(ieee_is_finite(aimag(x)))) return
        if (.not. all(ieee_is_finite(tapers))) return
        if (.not. all(ieee_is_finite(frequencies))) return
        if (present(bands)) then
            if (.not. all(ieee_is_finite(bands))) return
            if (any(bands(2,:) <= bands(1,:))) return
        end if
        hop = length
        if (present(noverlap)) then
            if (noverlap < 0 .or. noverlap >= length) return
            hop = length-noverlap
        end if
        norm = .true.
        if (present(normalized)) norm = normalized
        avg = .true.
        if (present(average)) avg = average
        allocate(weights(kmax), energies(kmax))
        weights = 1.0_dp
        if (present(taper_weights)) then
            if (size(taper_weights) /= kmax) return
            if (.not. all(ieee_is_finite(taper_weights))) return
            if (any(taper_weights < 0.0_dp)) return
            weights = taper_weights
        end if
        largest = maxval(weights)
        if (largest <= 0.0_dp) return
        weights = weights/largest
        total = sum(weights)
        if (.not. ieee_is_finite(total) .or. total <= 0.0_dp) return
        weights = weights/total
        energies = 1.0_dp
        do k = 1, kmax
            if (weights(k) <= 0.0_dp) cycle
            if (maxval(abs(tapers(:,k))) <= 0.0_dp) return
            energies(k) = sum(tapers(:,k)**2)
            if (.not. ieee_is_finite(energies(k)) .or. energies(k) <= 0.0_dp) return
        end do

        nw = 1+(size(x,1)-length)/hop
        allocate(out%frequencies(nf+nb),out%times(nw),out%power(nf+nb,nw,nc,nc))
        out%frequencies(1:nf) = frequencies
        if (present(bands)) out%frequencies(nf+1:) = 0.5_dp*(bands(1,:)+bands(2,:))
        out%times = [(real((w-1)*hop,dp)+0.5_dp*real(length,dp),w=1,nw)]
        out%power = (0.0_dp,0.0_dp)
        out%effective_tapers = 1.0_dp/sum(weights**2)
        allocate(samples(length,nc),transform(nc),phase(length),kernel(length,length), &
            work(length,nc))

        ! The band kernel K(t,s) integrates exp(-2*pi*i*f*(t-s)).
        ! Its exact sinc form avoids a frequency quadrature or FFT-grid bias.
        do w = 1, nw
            start = 1+(w-1)*hop
            do k = 1, kmax
                if (weights(k) <= 0.0_dp) cycle
                factor = weights(k)
                if (norm) factor = factor/energies(k)
                do j = 1, nc
                    samples(:,j) = x(start:start+length-1,j)*tapers(:,k)
                end do
                do f = 1, nf
                    phase = [(exp(cmplx(0.0_dp,-twopi*frequencies(f)*real(t-1,dp),dp)), &
                        t=1,length)]
                    do j = 1, nc
                        transform(j) = sum(samples(:,j)*phase)
                    end do
                    do j = 1, nc
                        do i = 1, nc
                            out%power(f,w,i,j) = out%power(f,w,i,j) + &
                                factor*transform(i)*conjg(transform(j))
                        end do
                    end do
                end do
                do f = 1, nb
                    width = bands(2,f)-bands(1,f)
                    midpoint = 0.5_dp*(bands(2,f)+bands(1,f))
                    do s = 1, length
                        do t = 1, length
                            lag = real(t-s,dp)
                            kernel(t,s) = sinc_pi(width*lag)* &
                                exp(cmplx(0.0_dp,-twopi*midpoint*lag,dp))
                            if (.not. avg) kernel(t,s) = width*kernel(t,s)
                        end do
                    end do
                    work = matmul(kernel,conjg(samples))
                    out%power(nf+f,w,:,:) = out%power(nf+f,w,:,:) + &
                        factor*matmul(transpose(samples),work)
                end do
            end do
        end do
        out%status = spectral_status_success
    end function fixed_taper_core

    function fixed_taper_real(x,tapers,frequencies,bands,taper_weights, &
            normalized,average,noverlap) result(out)
        real(dp), intent(in) :: x(:,:) !! finite real samples, time by channel
        real(dp), intent(in) :: tapers(:,:) !! sample by taper deterministic bank
        real(dp), intent(in) :: frequencies(:) !! point frequencies in cycles/sample
        real(dp), intent(in), optional :: bands(:,:) !! continuous bands, shape (2,nband)
        real(dp), intent(in), optional :: taper_weights(:) !! fixed relative taper weights
        logical, intent(in), optional :: normalized !! per-taper energy normalization, default true
        logical, intent(in), optional :: average !! average, not integral, band spectra, default true
        integer, intent(in), optional :: noverlap !! number of overlapping samples, default zero
        type(multitaper_spectral_matrix_spectrogram_result) :: out

        out = fixed_taper_core(cmplx(x,0.0_dp,dp),tapers,frequencies,bands, &
            taper_weights,normalized,average,noverlap)
    end function fixed_taper_real

    function fixed_taper_complex(x,tapers,frequencies,bands,taper_weights, &
            normalized,average,noverlap) result(out)
        complex(dp), intent(in) :: x(:,:) !! finite complex samples, time by channel
        real(dp), intent(in) :: tapers(:,:) !! sample by taper deterministic bank
        real(dp), intent(in) :: frequencies(:) !! point frequencies in cycles/sample
        real(dp), intent(in), optional :: bands(:,:) !! continuous bands, shape (2,nband)
        real(dp), intent(in), optional :: taper_weights(:) !! fixed relative taper weights
        logical, intent(in), optional :: normalized !! per-taper energy normalization, default true
        logical, intent(in), optional :: average !! average, not integral, band spectra, default true
        integer, intent(in), optional :: noverlap !! number of overlapping samples, default zero
        type(multitaper_spectral_matrix_spectrogram_result) :: out

        out = fixed_taper_core(x,tapers,frequencies,bands,taper_weights, &
            normalized,average,noverlap)
    end function fixed_taper_complex

end module scifort_signal_fixed_taper_spectrogram
