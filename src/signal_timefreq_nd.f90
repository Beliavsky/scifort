! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Rank-independent wrappers for legacy scipy.signal.stft/spectrogram.
! The selected signal axis becomes frequency, and time is appended last.
! All output writes are transactional: destination unchanged on failure.
module scifort_signal_timefreq_nd
    use scifort_kinds, only : dp
    use scifort_signal_spectral, only : stft_result, spectrogram_result, &
        spectral_status_success, spectral_status_invalid_input, &
        stft_real => stft, spectrogram_real => spectrogram, &
        spectral_real_detrend_callback
    use scifort_signal_spectral_complex, only : stft_complex, spectrogram_complex, &
        spectral_complex_detrend_callback
    use scifort_signal_spectral_nd, only : flatten_real, flatten_complex, &
        scatter_real, scatter_complex
    implicit none
    private
    public :: stft_into, spectrogram_into

    integer, parameter :: method_stft=1, method_psd=2, method_magnitude=3, &
        method_angle=4, method_phase=5, method_complex=6

    type :: timefreq_job
        integer :: method = 0
        logical :: input_real = .true.
        real(dp), allocatable :: fs
        character(:), allocatable :: window, detrend, scaling, boundary
        real(dp), allocatable :: window_parameter, window_parameter2
        real(dp), allocatable :: window_coefficients(:), window_values(:)
        integer, allocatable :: nperseg, noverlap, nfft
        logical, allocatable :: padded
        logical, allocatable :: return_onesided !! false for full FFT bins; complex always two-sided
        procedure(spectral_real_detrend_callback), pointer, nopass :: real_callback => null()
        procedure(spectral_complex_detrend_callback), pointer, nopass :: complex_callback => null()
    end type timefreq_job

    interface stft_into
        module procedure stft_into_real, stft_into_complex
    end interface
    interface spectrogram_into
        module procedure spectrogram_into_rr, spectrogram_into_cr, &
            spectrogram_into_rc, spectrogram_into_cc
    end interface
contains
    subroutine validate_shapes(d, dout, axis, a, status)
        integer, intent(in) :: d(:) !! input extents
        integer, intent(in) :: dout(:) !! output extents
        integer, intent(in), optional :: axis !! one-based or negative signal axis
        integer, intent(out) :: a !! resolved one-based signal axis
        integer, intent(out) :: status !! validation result
        integer :: nd
        status = spectral_status_invalid_input
        nd = size(d)
        a = nd
        if (nd < 1 .or. nd > 14 .or. size(dout) /= nd+1) return
        if (present(axis)) a = axis
        if (a < 0) a = nd+a+1
        if (a < 1 .or. a > nd) return
        if (any(d <= 0) .or. any(dout <= 0)) return
        if (any(d(:a-1) /= dout(:a-1))) return
        if (any(d(a+1:) /= dout(a+1:nd))) return
        status = spectral_status_success
    end subroutine validate_shapes

    subroutine configure_mode(job, mode, complex_output)
        type(timefreq_job), intent(inout) :: job !! spectrogram job configuration
        character(*), intent(in), optional :: mode !! psd/complex/magnitude/angle/phase
        logical, intent(in) :: complex_output !! true for complex output argument
        character(:), allocatable :: actual
        actual='psd'
        if (complex_output) actual='complex'
        if (present(mode)) actual=trim(mode)
        job%method=0
        select case(actual)
        case('psd')
            if (.not. complex_output) job%method=method_psd
        case('complex')
            if (complex_output) job%method=method_complex
        case('magnitude')
            if (.not. complex_output) job%method=method_magnitude
        case('angle')
            if (.not. complex_output) job%method=method_angle
        case('phase')
            if (.not. complex_output) job%method=method_phase
        end select
    end subroutine configure_mode

    subroutine compute_channel(x,job,values,frequencies,times,status)
        complex(dp), intent(in) :: x(:) !! signal channel
        type(timefreq_job), intent(in) :: job !! fixed transform options
        complex(dp), allocatable, intent(out) :: values(:,:) !! channel output
        real(dp), allocatable, intent(out) :: frequencies(:),times(:) !! coordinate vectors
        integer, intent(out) :: status !! zero for success
        type(stft_result) :: s
        type(spectrogram_result) :: p
        character(:), allocatable :: output_mode
        allocate(values(0,0), frequencies(0), times(0))
        status=spectral_status_invalid_input
        if (job%method < method_stft .or. job%method > method_complex) return
        select case(job%method)
        case(method_stft)
            if (job%input_real) then
                if (associated(job%real_callback)) then
                    s=stft_real(real(x,dp),job%fs,job%window,job%nperseg,job%noverlap, &
                        job%nfft,job%detrend,job%boundary,job%padded,job%scaling, &
                        window_parameter=job%window_parameter, &
                        window_parameter2=job%window_parameter2, &
                        window_coefficients=job%window_coefficients, &
                        window_values=job%window_values, &
                        return_onesided=job%return_onesided, &
                        detrend_callback=job%real_callback)
                else
                    s=stft_real(real(x,dp),job%fs,job%window,job%nperseg,job%noverlap, &
                        job%nfft,job%detrend,job%boundary,job%padded,job%scaling, &
                        window_parameter=job%window_parameter, &
                        window_parameter2=job%window_parameter2, &
                        window_coefficients=job%window_coefficients, &
                        window_values=job%window_values, &
                        return_onesided=job%return_onesided)
                end if
            else
                if (associated(job%complex_callback)) then
                    s=stft_complex(x,job%fs,job%window,job%nperseg,job%noverlap, &
                        job%nfft,job%detrend,job%boundary,job%padded,job%scaling, &
                        window_parameter=job%window_parameter, &
                        window_parameter2=job%window_parameter2, &
                        window_coefficients=job%window_coefficients, &
                        window_values=job%window_values, &
                        return_onesided=job%return_onesided, &
                        detrend_callback=job%complex_callback)
                else
                    s=stft_complex(x,job%fs,job%window,job%nperseg,job%noverlap, &
                        job%nfft,job%detrend,job%boundary,job%padded,job%scaling, &
                        window_parameter=job%window_parameter, &
                        window_parameter2=job%window_parameter2, &
                        window_coefficients=job%window_coefficients, &
                        window_values=job%window_values, &
                        return_onesided=job%return_onesided)
                end if
            end if
            if (s%status /= spectral_status_success) return
            values=s%values
            frequencies=s%frequencies
            times=s%times
        case(method_psd,method_complex,method_magnitude,method_angle,method_phase)
            ! Share all spectrogram defaults, scaling and FFT signed-zero
            ! conventions with the 1-D public implementation.  Previous
            ! adapters independently reconstructed STFTs, leading to subtle
            ! phase-branch and window/default discrepancies.
            select case(job%method)
            case(method_psd)
                output_mode='psd'
            case(method_complex)
                output_mode='complex'
            case(method_magnitude)
                output_mode='magnitude'
            case(method_angle)
                output_mode='angle'
            case(method_phase)
                output_mode='phase'
            end select
            if (job%input_real) then
                if (associated(job%real_callback)) then
                    p=spectrogram_real(real(x,dp),job%fs,job%window,job%nperseg, &
                        job%noverlap,job%nfft,job%detrend,job%scaling, &
                        window_parameter=job%window_parameter, &
                        window_parameter2=job%window_parameter2, &
                        window_coefficients=job%window_coefficients, &
                        window_values=job%window_values, &
                        return_onesided=job%return_onesided,mode=output_mode, &
                        detrend_callback=job%real_callback)
                else
                    p=spectrogram_real(real(x,dp),job%fs,job%window,job%nperseg, &
                        job%noverlap,job%nfft,job%detrend,job%scaling, &
                        window_parameter=job%window_parameter, &
                        window_parameter2=job%window_parameter2, &
                        window_coefficients=job%window_coefficients, &
                        window_values=job%window_values, &
                        return_onesided=job%return_onesided,mode=output_mode)
                end if
            else
                if (associated(job%complex_callback)) then
                    p=spectrogram_complex(x,job%fs,job%window,job%nperseg, &
                        job%noverlap,job%nfft,job%detrend,job%scaling, &
                        window_parameter=job%window_parameter, &
                        window_parameter2=job%window_parameter2, &
                        window_coefficients=job%window_coefficients, &
                        window_values=job%window_values, &
                        return_onesided=job%return_onesided,mode=output_mode, &
                        detrend_callback=job%complex_callback)
                else
                    p=spectrogram_complex(x,job%fs,job%window,job%nperseg, &
                        job%noverlap,job%nfft,job%detrend,job%scaling, &
                        window_parameter=job%window_parameter, &
                        window_parameter2=job%window_parameter2, &
                        window_coefficients=job%window_coefficients, &
                        window_values=job%window_values, &
                        return_onesided=job%return_onesided,mode=output_mode)
                end if
            end if
            if (p%status /= spectral_status_success) return
            if (job%method == method_complex) then
                values=p%values
            else
                values=cmplx(p%power,0.0_dp,dp)
            end if
            frequencies=p%frequencies
            times=p%times
        end select
        status=spectral_status_success
    end subroutine compute_channel

    subroutine process_lanes(x,d,dout,a,job,staged,frequencies,times,status)
        complex(dp), intent(in) :: x(:) !! flattened Fortran input
        integer, intent(in) :: d(:) !! signal extents
        integer, intent(in) :: dout(:) !! output extents
        integer, intent(in) :: a !! selected signal dimension
        type(timefreq_job), intent(in) :: job !! options
        complex(dp), intent(out) :: staged(:) !! staged output (no partial writes to caller)
        real(dp), allocatable, intent(out) :: frequencies(:),times(:) !! coordinates
        integer, intent(out) :: status !! result code
        complex(dp), allocatable :: lane(:), z(:,:)
        real(dp), allocatable :: f(:), t(:)
        integer :: inner,outer,n,nf,nt,i,j,k,b,idx,st
        status=spectral_status_invalid_input
        allocate(frequencies(0),times(0))
        inner=product(d(:a-1))
        outer=product(d(a+1:))
        n=d(a)
        nf=dout(a)
        nt=dout(size(dout))
        allocate(lane(n))
        do k=0,outer-1
            do i=0,inner-1
                do j=0,n-1
                    idx=k*inner*n+j*inner+i+1
                    lane(j+1)=x(idx)
                end do
                call compute_channel(lane,job,z,f,t,st)
                if (st/=spectral_status_success) return
                if (size(z,1)/=nf .or. size(z,2)/=nt) return
                if (size(f)/=nf .or. size(t)/=nt) return
                if (k==0 .and. i==0) then
                    frequencies=f
                    times=t
                else
                    if (any(f/=frequencies) .or. any(t/=times)) return
                end if
                do b=0,nt-1
                    do j=0,nf-1
                        idx=b*inner*nf*outer+k*inner*nf+j*inner+i+1
                        staged(idx)=z(j+1,b+1)
                    end do
                end do
            end do
        end do
        status=spectral_status_success
    end subroutine process_lanes

    subroutine apply_complex(x,d,result,frequencies,times,job,axis,status)
        complex(dp), intent(in) :: x(:) !! column-major flat samples
        integer, intent(in) :: d(:) !! input extents
        complex(dp), intent(inout) :: result(..) !! rank+1 destination array
        real(dp), allocatable, intent(out) :: frequencies(:),times(:) !! coordinates
        type(timefreq_job), intent(in) :: job !! transform options
        integer, intent(in), optional :: axis !! one-based or negative
        integer, intent(out) :: status !! zero if all lanes succeeded
        complex(dp), allocatable :: staged(:)
        integer :: a
        call validate_shapes(d,shape(result),axis,a,status)
        if (status==spectral_status_success) then
            allocate(staged(size(result)))
            call process_lanes(x,d,shape(result),a,job,staged,frequencies,times,status)
            if (status==spectral_status_success) call scatter_complex(result,staged)
        else
            allocate(frequencies(0),times(0))
        end if
        if (status/=spectral_status_success) then
            if (allocated(frequencies)) deallocate(frequencies)
            if (allocated(times)) deallocate(times)
            allocate(frequencies(0),times(0))
        end if
    end subroutine apply_complex

    subroutine apply_real(x,d,result,frequencies,times,job,axis,status)
        complex(dp), intent(in) :: x(:) !! column-major flat samples
        integer, intent(in) :: d(:) !! input extents
        real(dp), intent(inout) :: result(..) !! rank+1 destination array
        real(dp), allocatable, intent(out) :: frequencies(:),times(:) !! coordinates
        type(timefreq_job), intent(in) :: job !! transform options
        integer, intent(in), optional :: axis !! one-based or negative
        integer, intent(out) :: status !! zero if all lanes succeeded
        complex(dp), allocatable :: staged(:)
        integer :: a
        call validate_shapes(d,shape(result),axis,a,status)
        if (status==spectral_status_success) then
            allocate(staged(size(result)))
            call process_lanes(x,d,shape(result),a,job,staged,frequencies,times,status)
            if (status==spectral_status_success) call scatter_real(result,staged)
        else
            allocate(frequencies(0),times(0))
        end if
        if (status/=spectral_status_success) then
            if (allocated(frequencies)) deallocate(frequencies)
            if (allocated(times)) deallocate(times)
            allocate(frequencies(0),times(0))
        end if
    end subroutine apply_real

    subroutine stft_into_real(x,result,frequencies,times,fs,window,nperseg,noverlap,nfft,detrend, &
        boundary,padded,scaling,axis,status, &
        window_parameter, window_parameter2, window_coefficients, window_values, return_onesided, detrend_callback)
        real(dp), intent(in) :: x(..) !! rank 1..14 input samples
        complex(dp), intent(inout) :: result(..) !! rank+1 output
        real(dp), allocatable, intent(out) :: frequencies(:),times(:) !! frequency/time coordinates
        real(dp), intent(in), optional :: fs !! sampling frequency
        character(*), intent(in), optional :: window !! periodic window name
        integer, intent(in), optional :: nperseg !! samples in each segment
        integer, intent(in), optional :: noverlap !! overlapped samples
        integer, intent(in), optional :: nfft !! FFT size
        character(*), intent(in), optional :: detrend !! constant, linear, or none
        character(*), intent(in), optional :: boundary !! zeros, even, odd, constant, or none
        logical, intent(in), optional :: padded !! pad final segment
        character(*), intent(in), optional :: scaling !! spectrum/psd for STFT; density/spectrum for spectrogram
        integer, intent(in), optional :: axis !! one-based or negative signal axis
        integer, intent(out), optional :: status !! zero on success
        real(dp), intent(in), optional :: window_parameter !! named window shape parameter
        real(dp), intent(in), optional :: window_parameter2 !! second named window parameter
        real(dp), intent(in), optional :: window_coefficients(:) !! general cosine coefficients
        real(dp), intent(in), optional :: window_values(:) !! explicit real window samples
        logical, intent(in), optional :: return_onesided !! false computes complete two-sided FFT
        procedure(spectral_real_detrend_callback), optional :: detrend_callback !! segment callback; excludes detrend=
        type(timefreq_job) :: job
        complex(dp), allocatable :: flat(:)
        integer :: st
        job%method=method_stft
        job%input_real=.true.
        if (present(fs)) job%fs=fs
        if (present(window)) job%window=window
        if (present(window_parameter)) job%window_parameter = window_parameter
        if (present(window_parameter2)) job%window_parameter2 = window_parameter2
        if (present(window_coefficients)) job%window_coefficients = window_coefficients
        if (present(window_values)) job%window_values = window_values
        if (present(return_onesided)) job%return_onesided = return_onesided
        if (present(nperseg)) job%nperseg=nperseg
        if (present(noverlap)) job%noverlap=noverlap
        if (present(nfft)) job%nfft=nfft
        if (present(detrend)) job%detrend=detrend
        if (present(detrend_callback)) job%real_callback => detrend_callback
        if (present(scaling)) job%scaling=scaling
        if (present(boundary)) job%boundary=boundary
        if (present(padded)) job%padded=padded
        call flatten_real(x,flat)
        call apply_complex(flat,shape(x),result,frequencies,times,job,axis,st)
        if (present(status)) status=st
    end subroutine stft_into_real

    subroutine stft_into_complex(x,result,frequencies,times,fs,window,nperseg,noverlap,nfft,detrend, &
        boundary,padded,scaling,axis,status, &
        window_parameter, window_parameter2, window_coefficients, window_values, return_onesided, detrend_callback)
        complex(dp), intent(in) :: x(..) !! rank 1..14 input samples
        complex(dp), intent(inout) :: result(..) !! rank+1 output
        real(dp), allocatable, intent(out) :: frequencies(:),times(:) !! frequency/time coordinates
        real(dp), intent(in), optional :: fs !! sampling frequency
        character(*), intent(in), optional :: window !! periodic window name
        integer, intent(in), optional :: nperseg !! samples in each segment
        integer, intent(in), optional :: noverlap !! overlapped samples
        integer, intent(in), optional :: nfft !! FFT size
        character(*), intent(in), optional :: detrend !! constant, linear, or none
        character(*), intent(in), optional :: boundary !! zeros, even, odd, constant, or none
        logical, intent(in), optional :: padded !! pad final segment
        character(*), intent(in), optional :: scaling !! spectrum/psd for STFT; density/spectrum for spectrogram
        integer, intent(in), optional :: axis !! one-based or negative signal axis
        integer, intent(out), optional :: status !! zero on success
        real(dp), intent(in), optional :: window_parameter !! named window shape parameter
        real(dp), intent(in), optional :: window_parameter2 !! second named window parameter
        real(dp), intent(in), optional :: window_coefficients(:) !! general cosine coefficients
        real(dp), intent(in), optional :: window_values(:) !! explicit real window samples
        logical, intent(in), optional :: return_onesided !! false computes complete two-sided FFT
        procedure(spectral_complex_detrend_callback), optional :: detrend_callback !! segment callback; excludes detrend=
        type(timefreq_job) :: job
        complex(dp), allocatable :: flat(:)
        integer :: st
        job%method=method_stft
        job%input_real=.false.
        if (present(fs)) job%fs=fs
        if (present(window)) job%window=window
        if (present(window_parameter)) job%window_parameter = window_parameter
        if (present(window_parameter2)) job%window_parameter2 = window_parameter2
        if (present(window_coefficients)) job%window_coefficients = window_coefficients
        if (present(window_values)) job%window_values = window_values
        if (present(return_onesided)) job%return_onesided = return_onesided
        if (present(nperseg)) job%nperseg=nperseg
        if (present(noverlap)) job%noverlap=noverlap
        if (present(nfft)) job%nfft=nfft
        if (present(detrend)) job%detrend=detrend
        if (present(detrend_callback)) job%complex_callback => detrend_callback
        if (present(scaling)) job%scaling=scaling
        if (present(boundary)) job%boundary=boundary
        if (present(padded)) job%padded=padded
        call flatten_complex(x,flat)
        call apply_complex(flat,shape(x),result,frequencies,times,job,axis,st)
        if (present(status)) status=st
    end subroutine stft_into_complex

    subroutine spectrogram_into_rr(x,result,frequencies,times,fs,window,nperseg,noverlap,nfft,detrend,scaling,mode,axis,status, &
        window_parameter, window_parameter2, window_coefficients, window_values, return_onesided, detrend_callback)
        real(dp), intent(in) :: x(..) !! rank 1..14 input samples
        real(dp), intent(inout) :: result(..) !! rank+1 output
        real(dp), allocatable, intent(out) :: frequencies(:),times(:) !! frequency/time coordinates
        real(dp), intent(in), optional :: fs !! sampling frequency
        character(*), intent(in), optional :: window !! periodic window name
        integer, intent(in), optional :: nperseg !! samples in each segment
        integer, intent(in), optional :: noverlap !! overlapped samples
        integer, intent(in), optional :: nfft !! FFT size
        character(*), intent(in), optional :: detrend !! constant, linear, or none
        character(*), intent(in), optional :: scaling !! spectrum/psd for STFT; density/spectrum for spectrogram
        character(*), intent(in), optional :: mode !! psd/magnitude/angle/phase or complex
        integer, intent(in), optional :: axis !! one-based or negative signal axis
        integer, intent(out), optional :: status !! zero on success
        real(dp), intent(in), optional :: window_parameter !! named window shape parameter
        real(dp), intent(in), optional :: window_parameter2 !! second named window parameter
        real(dp), intent(in), optional :: window_coefficients(:) !! general cosine coefficients
        real(dp), intent(in), optional :: window_values(:) !! explicit real window samples
        logical, intent(in), optional :: return_onesided !! false computes complete two-sided FFT
        procedure(spectral_real_detrend_callback), optional :: detrend_callback !! segment callback; excludes detrend=
        type(timefreq_job) :: job
        complex(dp), allocatable :: flat(:)
        integer :: st
        call configure_mode(job,mode,.false.)
        job%input_real=.true.
        if (present(fs)) job%fs=fs
        if (present(window)) job%window=window
        if (present(window_parameter)) job%window_parameter = window_parameter
        if (present(window_parameter2)) job%window_parameter2 = window_parameter2
        if (present(window_coefficients)) job%window_coefficients = window_coefficients
        if (present(window_values)) job%window_values = window_values
        if (present(return_onesided)) job%return_onesided = return_onesided
        if (present(nperseg)) job%nperseg=nperseg
        if (present(noverlap)) job%noverlap=noverlap
        if (present(nfft)) job%nfft=nfft
        if (present(detrend)) job%detrend=detrend
        if (present(detrend_callback)) job%real_callback => detrend_callback
        if (present(scaling)) job%scaling=scaling
        call flatten_real(x,flat)
        call apply_real(flat,shape(x),result,frequencies,times,job,axis,st)
        if (present(status)) status=st
    end subroutine spectrogram_into_rr

    subroutine spectrogram_into_cr(x,result,frequencies,times,fs,window,nperseg,noverlap,nfft,detrend,scaling,mode,axis,status, &
        window_parameter, window_parameter2, window_coefficients, window_values, return_onesided, detrend_callback)
        complex(dp), intent(in) :: x(..) !! rank 1..14 input samples
        real(dp), intent(inout) :: result(..) !! rank+1 output
        real(dp), allocatable, intent(out) :: frequencies(:),times(:) !! frequency/time coordinates
        real(dp), intent(in), optional :: fs !! sampling frequency
        character(*), intent(in), optional :: window !! periodic window name
        integer, intent(in), optional :: nperseg !! samples in each segment
        integer, intent(in), optional :: noverlap !! overlapped samples
        integer, intent(in), optional :: nfft !! FFT size
        character(*), intent(in), optional :: detrend !! constant, linear, or none
        character(*), intent(in), optional :: scaling !! spectrum/psd for STFT; density/spectrum for spectrogram
        character(*), intent(in), optional :: mode !! psd/magnitude/angle/phase or complex
        integer, intent(in), optional :: axis !! one-based or negative signal axis
        integer, intent(out), optional :: status !! zero on success
        real(dp), intent(in), optional :: window_parameter !! named window shape parameter
        real(dp), intent(in), optional :: window_parameter2 !! second named window parameter
        real(dp), intent(in), optional :: window_coefficients(:) !! general cosine coefficients
        real(dp), intent(in), optional :: window_values(:) !! explicit real window samples
        logical, intent(in), optional :: return_onesided !! false computes complete two-sided FFT
        procedure(spectral_complex_detrend_callback), optional :: detrend_callback !! segment callback; excludes detrend=
        type(timefreq_job) :: job
        complex(dp), allocatable :: flat(:)
        integer :: st
        call configure_mode(job,mode,.false.)
        job%input_real=.false.
        if (present(fs)) job%fs=fs
        if (present(window)) job%window=window
        if (present(window_parameter)) job%window_parameter = window_parameter
        if (present(window_parameter2)) job%window_parameter2 = window_parameter2
        if (present(window_coefficients)) job%window_coefficients = window_coefficients
        if (present(window_values)) job%window_values = window_values
        if (present(return_onesided)) job%return_onesided = return_onesided
        if (present(nperseg)) job%nperseg=nperseg
        if (present(noverlap)) job%noverlap=noverlap
        if (present(nfft)) job%nfft=nfft
        if (present(detrend)) job%detrend=detrend
        if (present(detrend_callback)) job%complex_callback => detrend_callback
        if (present(scaling)) job%scaling=scaling
        call flatten_complex(x,flat)
        call apply_real(flat,shape(x),result,frequencies,times,job,axis,st)
        if (present(status)) status=st
    end subroutine spectrogram_into_cr

    subroutine spectrogram_into_rc(x,result,frequencies,times,fs,window,nperseg,noverlap,nfft,detrend,scaling,mode,axis,status, &
        window_parameter, window_parameter2, window_coefficients, window_values, return_onesided, detrend_callback)
        real(dp), intent(in) :: x(..) !! rank 1..14 input samples
        complex(dp), intent(inout) :: result(..) !! rank+1 output
        real(dp), allocatable, intent(out) :: frequencies(:),times(:) !! frequency/time coordinates
        real(dp), intent(in), optional :: fs !! sampling frequency
        character(*), intent(in), optional :: window !! periodic window name
        integer, intent(in), optional :: nperseg !! samples in each segment
        integer, intent(in), optional :: noverlap !! overlapped samples
        integer, intent(in), optional :: nfft !! FFT size
        character(*), intent(in), optional :: detrend !! constant, linear, or none
        character(*), intent(in), optional :: scaling !! spectrum/psd for STFT; density/spectrum for spectrogram
        character(*), intent(in), optional :: mode !! psd/magnitude/angle/phase or complex
        integer, intent(in), optional :: axis !! one-based or negative signal axis
        integer, intent(out), optional :: status !! zero on success
        real(dp), intent(in), optional :: window_parameter !! named window shape parameter
        real(dp), intent(in), optional :: window_parameter2 !! second named window parameter
        real(dp), intent(in), optional :: window_coefficients(:) !! general cosine coefficients
        real(dp), intent(in), optional :: window_values(:) !! explicit real window samples
        logical, intent(in), optional :: return_onesided !! false computes complete two-sided FFT
        procedure(spectral_real_detrend_callback), optional :: detrend_callback !! segment callback; excludes detrend=
        type(timefreq_job) :: job
        complex(dp), allocatable :: flat(:)
        integer :: st
        call configure_mode(job,mode,.true.)
        job%input_real=.true.
        if (present(fs)) job%fs=fs
        if (present(window)) job%window=window
        if (present(window_parameter)) job%window_parameter = window_parameter
        if (present(window_parameter2)) job%window_parameter2 = window_parameter2
        if (present(window_coefficients)) job%window_coefficients = window_coefficients
        if (present(window_values)) job%window_values = window_values
        if (present(return_onesided)) job%return_onesided = return_onesided
        if (present(nperseg)) job%nperseg=nperseg
        if (present(noverlap)) job%noverlap=noverlap
        if (present(nfft)) job%nfft=nfft
        if (present(detrend)) job%detrend=detrend
        if (present(detrend_callback)) job%real_callback => detrend_callback
        if (present(scaling)) job%scaling=scaling
        call flatten_real(x,flat)
        call apply_complex(flat,shape(x),result,frequencies,times,job,axis,st)
        if (present(status)) status=st
    end subroutine spectrogram_into_rc

    subroutine spectrogram_into_cc(x,result,frequencies,times,fs,window,nperseg,noverlap,nfft,detrend,scaling,mode,axis,status, &
        window_parameter, window_parameter2, window_coefficients, window_values, return_onesided, detrend_callback)
        complex(dp), intent(in) :: x(..) !! rank 1..14 input samples
        complex(dp), intent(inout) :: result(..) !! rank+1 output
        real(dp), allocatable, intent(out) :: frequencies(:),times(:) !! frequency/time coordinates
        real(dp), intent(in), optional :: fs !! sampling frequency
        character(*), intent(in), optional :: window !! periodic window name
        integer, intent(in), optional :: nperseg !! samples in each segment
        integer, intent(in), optional :: noverlap !! overlapped samples
        integer, intent(in), optional :: nfft !! FFT size
        character(*), intent(in), optional :: detrend !! constant, linear, or none
        character(*), intent(in), optional :: scaling !! spectrum/psd for STFT; density/spectrum for spectrogram
        character(*), intent(in), optional :: mode !! psd/magnitude/angle/phase or complex
        integer, intent(in), optional :: axis !! one-based or negative signal axis
        integer, intent(out), optional :: status !! zero on success
        real(dp), intent(in), optional :: window_parameter !! named window shape parameter
        real(dp), intent(in), optional :: window_parameter2 !! second named window parameter
        real(dp), intent(in), optional :: window_coefficients(:) !! general cosine coefficients
        real(dp), intent(in), optional :: window_values(:) !! explicit real window samples
        logical, intent(in), optional :: return_onesided !! false computes complete two-sided FFT
        procedure(spectral_complex_detrend_callback), optional :: detrend_callback !! segment callback; excludes detrend=
        type(timefreq_job) :: job
        complex(dp), allocatable :: flat(:)
        integer :: st
        call configure_mode(job,mode,.true.)
        job%input_real=.false.
        if (present(fs)) job%fs=fs
        if (present(window)) job%window=window
        if (present(window_parameter)) job%window_parameter = window_parameter
        if (present(window_parameter2)) job%window_parameter2 = window_parameter2
        if (present(window_coefficients)) job%window_coefficients = window_coefficients
        if (present(window_values)) job%window_values = window_values
        if (present(return_onesided)) job%return_onesided = return_onesided
        if (present(nperseg)) job%nperseg=nperseg
        if (present(noverlap)) job%noverlap=noverlap
        if (present(nfft)) job%nfft=nfft
        if (present(detrend)) job%detrend=detrend
        if (present(detrend_callback)) job%complex_callback => detrend_callback
        if (present(scaling)) job%scaling=scaling
        call flatten_complex(x,flat)
        call apply_complex(flat,shape(x),result,frequencies,times,job,axis,st)
        if (present(status)) status=st
    end subroutine spectrogram_into_cc

end module scifort_signal_timefreq_nd
