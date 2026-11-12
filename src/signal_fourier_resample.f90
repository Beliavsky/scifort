! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Fourier-domain resampling of periodic one-dimensional signals.
! Implements SciPy's positive/negative frequency truncation and Nyquist
! splitting/combining when the shorter transform has even length.
module scifort_signal_fourier_resample
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_signal_spectral, only : forward_fft, spectral_named_window, spectral_status_success
    implicit none
    private
    integer, parameter, public :: fourier_resample_status_success = 0
    integer, parameter, public :: fourier_resample_status_invalid_input = 1
    interface resample
        module procedure resample_nd_r2
        module procedure resample_nd_c2
        module procedure resample_nd_r3
        module procedure resample_nd_c3
        module procedure resample_real, resample_complex
    end interface resample
    public :: resample, real_frequency_window_callback, complex_frequency_window_callback
    public :: evaluate_real_frequency_window, evaluate_complex_frequency_window

    abstract interface
        function real_frequency_window_callback(frequencies) result(weights)
            import dp
            real(dp), intent(in) :: frequencies(:) !! FFT frequencies in cycles per sample
            real(dp) :: weights(size(frequencies)) !! one real weight per frequency bin
        end function real_frequency_window_callback
        function complex_frequency_window_callback(frequencies) result(weights)
            import dp
            real(dp), intent(in) :: frequencies(:) !! FFT frequencies in cycles per sample
            complex(dp) :: weights(size(frequencies)) !! complex FFT-order weights
        end function complex_frequency_window_callback
    end interface
contains
    subroutine fft_frequency_grid(nx, frequencies)
        integer, intent(in) :: nx !! positive FFT length
        real(dp), intent(out) :: frequencies(nx) !! FFT-order frequencies in cycles per sample
        integer :: k
        do k = 0,nx-1
            if (k <= (nx-1)/2) then
                frequencies(k+1) = real(k,dp)/real(nx,dp)
            else
                frequencies(k+1) = real(k-nx,dp)/real(nx,dp)
            end if
        end do
    end subroutine fft_frequency_grid

    subroutine evaluate_real_frequency_window(nx, callback, weights, status)
        integer, intent(in) :: nx !! positive FFT length
        procedure(real_frequency_window_callback) :: callback
        real(dp), allocatable, intent(out) :: weights(:) !! sampled real response
        integer, intent(out) :: status !! zero if all sampled weights are finite
        real(dp), allocatable :: frequencies(:)
        status = fourier_resample_status_invalid_input
        if (nx < 1) return
        allocate(frequencies(nx))
        call fft_frequency_grid(nx,frequencies)
        weights = callback(frequencies)
        if (.not. all(ieee_is_finite(weights))) return
        status = fourier_resample_status_success
    end subroutine evaluate_real_frequency_window

    subroutine evaluate_complex_frequency_window(nx, callback, weights, status)
        integer, intent(in) :: nx !! positive FFT length
        procedure(complex_frequency_window_callback) :: callback
        complex(dp), allocatable, intent(out) :: weights(:) !! sampled complex response
        integer, intent(out) :: status !! zero if all sampled weights are finite
        real(dp), allocatable :: frequencies(:)
        status = fourier_resample_status_invalid_input
        if (nx < 1) return
        allocate(frequencies(nx))
        call fft_frequency_grid(nx,frequencies)
        weights = callback(frequencies)
        if (.not. all(ieee_is_finite(real(weights,dp)))) return
        if (.not. all(ieee_is_finite(aimag(weights)))) return
        status = fourier_resample_status_success
    end subroutine evaluate_complex_frequency_window

    function resample_complex(x, num, window, t, output_times, domain, status, &
            frequency_window, complex_window, real_input, window_function, complex_window_function, &
            window_parameter, window_parameter2, window_coefficients) result(y)
        complex(dp), intent(in) :: x(:) !! periodic input samples or Fourier coefficients
        integer, intent(in) :: num !! desired positive output sample count
        character(*), intent(in), optional :: window !! periodic spectral taper, optional
        real(dp), intent(in), optional :: t(:) !! input time grid (first two values used)
        real(dp), allocatable, intent(out), optional :: output_times(:) !! output grid
        character(*), intent(in), optional :: domain !! time (default) or freq
        integer, intent(out), optional :: status !! zero on success
        real(dp), intent(in), optional :: frequency_window(:) !! FFT-order real bin weights
        complex(dp), intent(in), optional :: complex_window(:) !! FFT-order complex weights
        logical, intent(in), optional :: real_input !! internal flag for real-input folding
        procedure(real_frequency_window_callback), optional :: window_function
        procedure(complex_frequency_window_callback), optional :: complex_window_function
        real(dp), intent(in), optional :: window_parameter !! first named window shape parameter
        real(dp), intent(in), optional :: window_parameter2 !! second named window shape parameter
        real(dp), intent(in), optional :: window_coefficients(:) !! general-cosine series coefficients
        complex(dp), allocatable :: y(:), spectrum(:), retained(:), bin_weights(:)
        real(dp), allocatable :: win(:), evaluated_real(:)
        complex(dp), allocatable :: evaluated_complex(:)
        integer :: nx, m, half, nlow, k, st, partner
        logical :: fold_real
        character(:), allocatable :: data_domain

        if (present(status)) status = fourier_resample_status_invalid_input
        allocate(y(0))
        if (present(output_times)) allocate(output_times(0))
        nx = size(x)
        if (num < 1 .or. nx < 1) return
        if (.not. all(ieee_is_finite(real(x,dp)))) return
        if (.not. all(ieee_is_finite(aimag(x)))) return
        if (present(t)) then
            if (size(t) < 2) return
            if (.not. all(ieee_is_finite(t(:2)))) return
        end if
        data_domain = 'time'
        if (present(domain)) data_domain = trim(domain)
        if (data_domain /= 'time' .and. data_domain /= 'freq') return
        fold_real = .false.
        if (present(real_input)) fold_real = real_input
        if (fold_real .and. data_domain /= 'time') return
        ! Window specifications are mutually exclusive, including callbacks.
        if (present(window_parameter) .or. present(window_parameter2) .or. &
                present(window_coefficients)) then
            if (.not. present(window)) return
        end if
        if (present(window)) then
            if (present(frequency_window) .or. present(complex_window) .or. &
                present(window_function) .or. present(complex_window_function)) return
        end if
        if (present(frequency_window)) then
            if (present(complex_window) .or. present(window_function) .or. &
                present(complex_window_function)) return
        end if
        if (present(complex_window)) then
            if (present(window_function) .or. present(complex_window_function)) return
        end if
        if (present(window_function) .and. present(complex_window_function)) return
        if (fold_real .and. present(complex_window_function)) return
        if (present(frequency_window)) then
            if (size(frequency_window) /= nx) return
            if (.not. all(ieee_is_finite(frequency_window))) return
        end if
        if (present(complex_window)) then
            if (fold_real .or. size(complex_window) /= nx) return
            if (.not. all(ieee_is_finite(real(complex_window,dp)))) return
            if (.not. all(ieee_is_finite(aimag(complex_window)))) return
        end if
        if (present(window_function)) then
            call evaluate_real_frequency_window(nx,window_function,evaluated_real,st)
            if (st /= fourier_resample_status_success) return
        end if
        if (present(complex_window_function)) then
            call evaluate_complex_frequency_window(nx,complex_window_function,evaluated_complex,st)
            if (st /= fourier_resample_status_success) return
        end if
        spectrum = x
        if (data_domain == 'time') call forward_fft(spectrum)
        if (present(window) .or. present(frequency_window) .or. present(complex_window) .or. &
                present(window_function) .or. present(complex_window_function)) then
            allocate(bin_weights(nx))
            if (present(window)) then
                win = spectral_named_window(trim(window),nx,st, &
                    parameter=window_parameter, parameter2=window_parameter2, &
                    coefficients=window_coefficients, fftbins=.true.)
                if (st /= spectral_status_success) return
                ! SciPy's fftshift of the periodic window, FFT bin order.
                do k = 0, nx-1
                    bin_weights(k+1) = cmplx(win(modulo(k+(nx+1)/2,nx)+1),0.0_dp,dp)
                end do
            else if (present(frequency_window)) then
                bin_weights = cmplx(frequency_window,0.0_dp,dp)
            else if (present(complex_window)) then
                bin_weights = complex_window
            else if (present(window_function)) then
                bin_weights = cmplx(evaluated_real,0.0_dp,dp)
            else
                bin_weights = evaluated_complex
            end if
            if (fold_real) then
                ! Real-input resampling uses a one-sided transform. SciPy
                ! averages conjugate-frequency weights before applying it.
                ! Read from unmodified entries: k and (-k modulo nx).
                win = real(bin_weights,dp)
                do k = 0, nx-1
                    partner = modulo(nx-k,nx)
                    bin_weights(k+1) = cmplx(0.5_dp*(win(k+1)+win(partner+1)),0.0_dp,dp)
                end do
            end if
            spectrum = spectrum*bin_weights
        end if
        allocate(retained(num))
        retained = cmplx(0.0_dp,0.0_dp,dp)
        m = min(nx,num)
        half = m/2
        nlow = (m-1)/2
        retained(1:nlow+1) = spectrum(1:nlow+1)
        do k = 1, nlow
            retained(num-k+1) = spectrum(nx-k+1)
        end do
        if (mod(m,2) == 0) then
            if (num > nx) then
                retained(half+1) = spectrum(half+1)/2.0_dp
                retained(num-half+1) = spectrum(half+1)/2.0_dp
            else if (num < nx) then
                retained(half+1) = spectrum(half+1) + spectrum(nx-half+1)
            else
                retained(half+1) = spectrum(half+1)
            end if
        end if
        ! Unnormalized inverse FFT divided by the input length, not the output length.
        retained = conjg(retained)
        call forward_fft(retained)
        y = conjg(retained)/real(nx,dp)
        if (present(output_times)) then
            if (present(t)) then
                deallocate(output_times)
                allocate(output_times(num))
                do k = 1,num
                    output_times(k) = t(1) + (t(2)-t(1)) * real(nx,dp)/real(num,dp) * real(k-1,dp)
                end do
            end if
        end if
        if (present(status)) status = fourier_resample_status_success
    end function resample_complex

    function resample_real(x, num, window, t, output_times, status, frequency_window, window_function, &
            window_parameter, window_parameter2, window_coefficients) result(y)
        real(dp), intent(in) :: x(:) !! periodic real samples
        integer, intent(in) :: num !! desired number of samples
        character(*), intent(in), optional :: window !! optional spectral taper
        real(dp), intent(in), optional :: t(:) !! input time grid, first two entries used
        real(dp), allocatable, intent(out), optional :: output_times(:) !! output timestamps
        integer, intent(out), optional :: status !! zero on success
        real(dp), intent(in), optional :: frequency_window(:) !! FFT-order real bin weights
        procedure(real_frequency_window_callback), optional :: window_function
        real(dp), intent(in), optional :: window_parameter !! first named window shape parameter
        real(dp), intent(in), optional :: window_parameter2 !! second named window shape parameter
        real(dp), intent(in), optional :: window_coefficients(:) !! general-cosine series coefficients
        real(dp), allocatable :: y(:)
        complex(dp), allocatable :: z(:)
        z = resample_complex(cmplx(x,0.0_dp,dp),num,window,t,output_times,status=status, &
            frequency_window=frequency_window,real_input=.true.,window_function=window_function, &
            window_parameter=window_parameter,window_parameter2=window_parameter2, &
            window_coefficients=window_coefficients)
        y = real(z,dp)
    end function resample_real

    ! BEGIN GENERATED MULTIRATE ND WRAPPERS
    function resample_nd_r2(x, num, window, t, output_times, axis, status, frequency_window, window_function, &
            window_parameter, window_parameter2, window_coefficients) result(y)
        real(dp), intent(in) :: x(:,:) !! input samples
        integer, intent(in) :: num !! output length along the transform axis
        character(*), intent(in), optional :: window !! periodic spectral window
        real(dp), intent(in), optional :: t(:) !! original sample times
        real(dp), allocatable, intent(out), optional :: output_times(:) !! resampled times
        integer, intent(in), optional :: axis !! one-based or negative axis
        integer, intent(out), optional :: status !! zero on success
        real(dp), intent(in), optional :: frequency_window(:) !! FFT-order real weights
        procedure(real_frequency_window_callback), optional :: window_function
        real(dp), intent(in), optional :: window_parameter !! first named window shape parameter
        real(dp), intent(in), optional :: window_parameter2 !! second named window shape parameter
        real(dp), intent(in), optional :: window_coefficients(:) !! general-cosine series coefficients
        real(dp), allocatable :: cached_real_window(:)
        real(dp), allocatable :: y(:,:), line(:)
        integer :: a, code, sizes(2)
        integer :: i, j
        allocate(y(0,0))
        if (present(status)) status = fourier_resample_status_invalid_input
        if (any(shape(x) == 0)) return
        a = 1
        if (present(axis)) a = axis
        if (a < 0) a = a + 2 + 1
        if (a < 1 .or. a > 2) return
        if (present(window_function)) then
            if (present(window) .or. present(frequency_window)) return
            call evaluate_real_frequency_window(size(x,a),window_function,cached_real_window,code)
            if (code /= fourier_resample_status_success) return
        else if (present(frequency_window)) then
            cached_real_window = frequency_window
        end if
        select case (a)
        case (1)
            line = resample_real(x(:,1), num, window, t, output_times, status=code, frequency_window=cached_real_window, &
                            window_parameter=window_parameter, window_parameter2=window_parameter2, &
                            window_coefficients=window_coefficients)
        case (2)
            line = resample_real(x(1,:), num, window, t, output_times, status=code, frequency_window=cached_real_window, &
                            window_parameter=window_parameter, window_parameter2=window_parameter2, &
                            window_coefficients=window_coefficients)
        end select
        if (code /= fourier_resample_status_success .or. size(line) < 1) return
        sizes = shape(x)
        sizes(a) = size(line)
        deallocate(y)
        allocate(y(sizes(1),sizes(2)))
        select case (a)
        case (1)
            do i = 1, size(x,2)
                    line = resample_real(x(:,i), num, window=window, status=code, frequency_window=cached_real_window, &
                            window_parameter=window_parameter, window_parameter2=window_parameter2, &
                            window_coefficients=window_coefficients)
                    if (code /= fourier_resample_status_success .or. size(line) /= sizes(a)) then
                        deallocate(y)
                        allocate(y(0,0))
                        return
                    end if
                    y(:,i) = line
            end do
        case (2)
            do i = 1, size(x,1)
                    line = resample_real(x(i,:), num, window=window, status=code, frequency_window=cached_real_window, &
                            window_parameter=window_parameter, window_parameter2=window_parameter2, &
                            window_coefficients=window_coefficients)
                    if (code /= fourier_resample_status_success .or. size(line) /= sizes(a)) then
                        deallocate(y)
                        allocate(y(0,0))
                        return
                    end if
                    y(i,:) = line
            end do
        end select
        if (present(status)) status = fourier_resample_status_success
    end function resample_nd_r2

    function resample_nd_c2(x, num, window, t, output_times, domain, axis, status, frequency_window, complex_window, &
            window_function, complex_window_function, &
            window_parameter, window_parameter2, window_coefficients) result(y)
        complex(dp), intent(in) :: x(:,:) !! input samples
        integer, intent(in) :: num !! output length along the transform axis
        character(*), intent(in), optional :: window !! periodic spectral window
        real(dp), intent(in), optional :: t(:) !! original sample times
        real(dp), allocatable, intent(out), optional :: output_times(:) !! resampled times
        character(*), intent(in), optional :: domain !! time or freq
        integer, intent(in), optional :: axis !! one-based or negative axis
        integer, intent(out), optional :: status !! zero on success
        real(dp), intent(in), optional :: frequency_window(:) !! FFT-order real weights
        complex(dp), intent(in), optional :: complex_window(:) !! FFT-order complex weights
        procedure(real_frequency_window_callback), optional :: window_function
        procedure(complex_frequency_window_callback), optional :: complex_window_function
        real(dp), intent(in), optional :: window_parameter !! first named window shape parameter
        real(dp), intent(in), optional :: window_parameter2 !! second named window shape parameter
        real(dp), intent(in), optional :: window_coefficients(:) !! general-cosine series coefficients
        real(dp), allocatable :: cached_real_window(:)
        complex(dp), allocatable :: cached_complex_window(:)
        complex(dp), allocatable :: y(:,:), line(:)
        integer :: a, code, sizes(2)
        integer :: i, j
        allocate(y(0,0))
        if (present(status)) status = fourier_resample_status_invalid_input
        if (any(shape(x) == 0)) return
        a = 1
        if (present(axis)) a = axis
        if (a < 0) a = a + 2 + 1
        if (a < 1 .or. a > 2) return
        if (present(window_function)) then
            if (present(window) .or. present(frequency_window) .or. &
                present(complex_window) .or. present(complex_window_function)) return
            call evaluate_real_frequency_window(size(x,a),window_function,cached_real_window,code)
            if (code /= fourier_resample_status_success) return
        else if (present(complex_window_function)) then
            if (present(window) .or. present(frequency_window) .or. present(complex_window)) return
            call evaluate_complex_frequency_window(size(x,a),complex_window_function,cached_complex_window,code)
            if (code /= fourier_resample_status_success) return
        else
            if (present(frequency_window)) cached_real_window = frequency_window
            if (present(complex_window)) cached_complex_window = complex_window
        end if
        select case (a)
        case (1)
            line = resample_complex(x(:,1), num, window, t, output_times, domain, status=code, &
                            frequency_window=cached_real_window, complex_window=cached_complex_window, &
                            window_parameter=window_parameter, window_parameter2=window_parameter2, &
                            window_coefficients=window_coefficients)
        case (2)
            line = resample_complex(x(1,:), num, window, t, output_times, domain, status=code, &
                            frequency_window=cached_real_window, complex_window=cached_complex_window, &
                            window_parameter=window_parameter, window_parameter2=window_parameter2, &
                            window_coefficients=window_coefficients)
        end select
        if (code /= fourier_resample_status_success .or. size(line) < 1) return
        sizes = shape(x)
        sizes(a) = size(line)
        deallocate(y)
        allocate(y(sizes(1),sizes(2)))
        select case (a)
        case (1)
            do i = 1, size(x,2)
                    line = resample_complex(x(:,i), num, window=window, domain=domain, status=code, &
                            frequency_window=cached_real_window, complex_window=cached_complex_window, &
                            window_parameter=window_parameter, window_parameter2=window_parameter2, &
                            window_coefficients=window_coefficients)
                    if (code /= fourier_resample_status_success .or. size(line) /= sizes(a)) then
                        deallocate(y)
                        allocate(y(0,0))
                        return
                    end if
                    y(:,i) = line
            end do
        case (2)
            do i = 1, size(x,1)
                    line = resample_complex(x(i,:), num, window=window, domain=domain, status=code, &
                            frequency_window=cached_real_window, complex_window=cached_complex_window, &
                            window_parameter=window_parameter, window_parameter2=window_parameter2, &
                            window_coefficients=window_coefficients)
                    if (code /= fourier_resample_status_success .or. size(line) /= sizes(a)) then
                        deallocate(y)
                        allocate(y(0,0))
                        return
                    end if
                    y(i,:) = line
            end do
        end select
        if (present(status)) status = fourier_resample_status_success
    end function resample_nd_c2

    function resample_nd_r3(x, num, window, t, output_times, axis, status, frequency_window, window_function, &
            window_parameter, window_parameter2, window_coefficients) result(y)
        real(dp), intent(in) :: x(:,:,:) !! input samples
        integer, intent(in) :: num !! output length along the transform axis
        character(*), intent(in), optional :: window !! periodic spectral window
        real(dp), intent(in), optional :: t(:) !! original sample times
        real(dp), allocatable, intent(out), optional :: output_times(:) !! resampled times
        integer, intent(in), optional :: axis !! one-based or negative axis
        integer, intent(out), optional :: status !! zero on success
        real(dp), intent(in), optional :: frequency_window(:) !! FFT-order real weights
        procedure(real_frequency_window_callback), optional :: window_function
        real(dp), intent(in), optional :: window_parameter !! first named window shape parameter
        real(dp), intent(in), optional :: window_parameter2 !! second named window shape parameter
        real(dp), intent(in), optional :: window_coefficients(:) !! general-cosine series coefficients
        real(dp), allocatable :: cached_real_window(:)
        real(dp), allocatable :: y(:,:,:), line(:)
        integer :: a, code, sizes(3)
        integer :: i, j, k
        allocate(y(0,0,0))
        if (present(status)) status = fourier_resample_status_invalid_input
        if (any(shape(x) == 0)) return
        a = 1
        if (present(axis)) a = axis
        if (a < 0) a = a + 3 + 1
        if (a < 1 .or. a > 3) return
        if (present(window_function)) then
            if (present(window) .or. present(frequency_window)) return
            call evaluate_real_frequency_window(size(x,a),window_function,cached_real_window,code)
            if (code /= fourier_resample_status_success) return
        else if (present(frequency_window)) then
            cached_real_window = frequency_window
        end if
        select case (a)
        case (1)
            line = resample_real(x(:,1,1), num, window, t, output_times, status=code, frequency_window=cached_real_window, &
                            window_parameter=window_parameter, window_parameter2=window_parameter2, &
                            window_coefficients=window_coefficients)
        case (2)
            line = resample_real(x(1,:,1), num, window, t, output_times, status=code, frequency_window=cached_real_window, &
                            window_parameter=window_parameter, window_parameter2=window_parameter2, &
                            window_coefficients=window_coefficients)
        case (3)
            line = resample_real(x(1,1,:), num, window, t, output_times, status=code, frequency_window=cached_real_window, &
                            window_parameter=window_parameter, window_parameter2=window_parameter2, &
                            window_coefficients=window_coefficients)
        end select
        if (code /= fourier_resample_status_success .or. size(line) < 1) return
        sizes = shape(x)
        sizes(a) = size(line)
        deallocate(y)
        allocate(y(sizes(1),sizes(2),sizes(3)))
        select case (a)
        case (1)
            do i = 1, size(x,2)
                do j = 1, size(x,3)
                        line = resample_real(x(:,i,j), num, window=window, status=code, frequency_window=cached_real_window, &
                            window_parameter=window_parameter, window_parameter2=window_parameter2, &
                            window_coefficients=window_coefficients)
                        if (code /= fourier_resample_status_success .or. size(line) /= sizes(a)) then
                            deallocate(y)
                            allocate(y(0,0,0))
                            return
                        end if
                        y(:,i,j) = line
                end do
            end do
        case (2)
            do i = 1, size(x,1)
                do j = 1, size(x,3)
                        line = resample_real(x(i,:,j), num, window=window, status=code, frequency_window=cached_real_window, &
                            window_parameter=window_parameter, window_parameter2=window_parameter2, &
                            window_coefficients=window_coefficients)
                        if (code /= fourier_resample_status_success .or. size(line) /= sizes(a)) then
                            deallocate(y)
                            allocate(y(0,0,0))
                            return
                        end if
                        y(i,:,j) = line
                end do
            end do
        case (3)
            do i = 1, size(x,1)
                do j = 1, size(x,2)
                        line = resample_real(x(i,j,:), num, window=window, status=code, frequency_window=cached_real_window, &
                            window_parameter=window_parameter, window_parameter2=window_parameter2, &
                            window_coefficients=window_coefficients)
                        if (code /= fourier_resample_status_success .or. size(line) /= sizes(a)) then
                            deallocate(y)
                            allocate(y(0,0,0))
                            return
                        end if
                        y(i,j,:) = line
                end do
            end do
        end select
        if (present(status)) status = fourier_resample_status_success
    end function resample_nd_r3

    function resample_nd_c3(x, num, window, t, output_times, domain, axis, status, frequency_window, complex_window, &
            window_function, complex_window_function, &
            window_parameter, window_parameter2, window_coefficients) result(y)
        complex(dp), intent(in) :: x(:,:,:) !! input samples
        integer, intent(in) :: num !! output length along the transform axis
        character(*), intent(in), optional :: window !! periodic spectral window
        real(dp), intent(in), optional :: t(:) !! original sample times
        real(dp), allocatable, intent(out), optional :: output_times(:) !! resampled times
        character(*), intent(in), optional :: domain !! time or freq
        integer, intent(in), optional :: axis !! one-based or negative axis
        integer, intent(out), optional :: status !! zero on success
        real(dp), intent(in), optional :: frequency_window(:) !! FFT-order real weights
        complex(dp), intent(in), optional :: complex_window(:) !! FFT-order complex weights
        procedure(real_frequency_window_callback), optional :: window_function
        procedure(complex_frequency_window_callback), optional :: complex_window_function
        real(dp), intent(in), optional :: window_parameter !! first named window shape parameter
        real(dp), intent(in), optional :: window_parameter2 !! second named window shape parameter
        real(dp), intent(in), optional :: window_coefficients(:) !! general-cosine series coefficients
        real(dp), allocatable :: cached_real_window(:)
        complex(dp), allocatable :: cached_complex_window(:)
        complex(dp), allocatable :: y(:,:,:), line(:)
        integer :: a, code, sizes(3)
        integer :: i, j, k
        allocate(y(0,0,0))
        if (present(status)) status = fourier_resample_status_invalid_input
        if (any(shape(x) == 0)) return
        a = 1
        if (present(axis)) a = axis
        if (a < 0) a = a + 3 + 1
        if (a < 1 .or. a > 3) return
        if (present(window_function)) then
            if (present(window) .or. present(frequency_window) .or. &
                present(complex_window) .or. present(complex_window_function)) return
            call evaluate_real_frequency_window(size(x,a),window_function,cached_real_window,code)
            if (code /= fourier_resample_status_success) return
        else if (present(complex_window_function)) then
            if (present(window) .or. present(frequency_window) .or. present(complex_window)) return
            call evaluate_complex_frequency_window(size(x,a),complex_window_function,cached_complex_window,code)
            if (code /= fourier_resample_status_success) return
        else
            if (present(frequency_window)) cached_real_window = frequency_window
            if (present(complex_window)) cached_complex_window = complex_window
        end if
        select case (a)
        case (1)
            line = resample_complex(x(:,1,1), num, window, t, output_times, domain, status=code, &
                            frequency_window=cached_real_window, complex_window=cached_complex_window, &
                            window_parameter=window_parameter, window_parameter2=window_parameter2, &
                            window_coefficients=window_coefficients)
        case (2)
            line = resample_complex(x(1,:,1), num, window, t, output_times, domain, status=code, &
                            frequency_window=cached_real_window, complex_window=cached_complex_window, &
                            window_parameter=window_parameter, window_parameter2=window_parameter2, &
                            window_coefficients=window_coefficients)
        case (3)
            line = resample_complex(x(1,1,:), num, window, t, output_times, domain, status=code, &
                            frequency_window=cached_real_window, complex_window=cached_complex_window, &
                            window_parameter=window_parameter, window_parameter2=window_parameter2, &
                            window_coefficients=window_coefficients)
        end select
        if (code /= fourier_resample_status_success .or. size(line) < 1) return
        sizes = shape(x)
        sizes(a) = size(line)
        deallocate(y)
        allocate(y(sizes(1),sizes(2),sizes(3)))
        select case (a)
        case (1)
            do i = 1, size(x,2)
                do j = 1, size(x,3)
                        line = resample_complex(x(:,i,j), num, window=window, domain=domain, status=code, &
                            frequency_window=cached_real_window, complex_window=cached_complex_window, &
                            window_parameter=window_parameter, window_parameter2=window_parameter2, &
                            window_coefficients=window_coefficients)
                        if (code /= fourier_resample_status_success .or. size(line) /= sizes(a)) then
                            deallocate(y)
                            allocate(y(0,0,0))
                            return
                        end if
                        y(:,i,j) = line
                end do
            end do
        case (2)
            do i = 1, size(x,1)
                do j = 1, size(x,3)
                        line = resample_complex(x(i,:,j), num, window=window, domain=domain, status=code, &
                            frequency_window=cached_real_window, complex_window=cached_complex_window, &
                            window_parameter=window_parameter, window_parameter2=window_parameter2, &
                            window_coefficients=window_coefficients)
                        if (code /= fourier_resample_status_success .or. size(line) /= sizes(a)) then
                            deallocate(y)
                            allocate(y(0,0,0))
                            return
                        end if
                        y(i,:,j) = line
                end do
            end do
        case (3)
            do i = 1, size(x,1)
                do j = 1, size(x,2)
                        line = resample_complex(x(i,j,:), num, window=window, domain=domain, status=code, &
                            frequency_window=cached_real_window, complex_window=cached_complex_window, &
                            window_parameter=window_parameter, window_parameter2=window_parameter2, &
                            window_coefficients=window_coefficients)
                        if (code /= fourier_resample_status_success .or. size(line) /= sizes(a)) then
                            deallocate(y)
                            allocate(y(0,0,0))
                            return
                        end if
                        y(i,j,:) = line
                end do
            end do
        end select
        if (present(status)) status = fourier_resample_status_success
    end function resample_nd_c3
    ! END GENERATED MULTIRATE ND WRAPPERS
end module scifort_signal_fourier_resample
