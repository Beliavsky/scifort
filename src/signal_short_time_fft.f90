! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Centered, invertible short-time Fourier analysis following scipy.signal.ShortTimeFFT.
! Independent implementation using the existing SciFort FFT kernel.
module scifort_signal_short_time_fft
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_signal_spectral, only : forward_fft, get_window, spectral_status_success
    implicit none
    private

    integer, parameter, public :: short_fft_success = 0
    integer, parameter, public :: short_fft_invalid_input = 1
    integer, parameter, public :: short_fft_noninvertible = 2
    real(dp), parameter :: pi = acos(-1.0_dp)

    ! A callable detrend transforms each unwindowed time slice, as in SciPy's
    ! ShortTimeFFT.stft_detrend. Separate real and complex callbacks preserve
    ! the input signal type and keep the public Fortran interface type safe.
    abstract interface
        subroutine stft_real_detrend_callback(samples)
            import dp
            real(dp), intent(inout) :: samples(:) !! one unwindowed real STFT slice
        end subroutine stft_real_detrend_callback
        subroutine stft_complex_detrend_callback(samples)
            import dp
            complex(dp), intent(inout) :: samples(:) !! one unwindowed complex STFT slice
        end subroutine stft_complex_detrend_callback
    end interface
    public :: stft_real_detrend_callback, stft_complex_detrend_callback

    type, public :: short_time_fft_result
        complex(dp), allocatable :: values(:,:) !! FFT bins by centered time slices
        real(dp), allocatable :: frequencies(:) !! FFT frequencies, in hertz
        real(dp), allocatable :: times(:) !! time centers, in seconds
        integer :: p0 = 0 !! first time-slice index (zero-based)
        integer :: status = short_fft_invalid_input !! zero on success
    end type short_time_fft_result

    type, public :: short_time_spectrogram_result
        real(dp), allocatable :: power(:,:) !! squared magnitude of STFT
        real(dp), allocatable :: frequencies(:) !! frequency sampling points
        real(dp), allocatable :: times(:) !! centered time sampling points
        integer :: p0 = 0 !! first zero-based slice index
        integer :: status = short_fft_invalid_input !! zero on success
    end type short_time_spectrogram_result

    ! For multidimensional input, the selected signal axis is replaced by
    ! a frequency axis and the final output axis contains time slices.
    type, public :: short_time_fft_result_3d
        complex(dp), allocatable :: values(:,:,:) !! frequency and time with one batch axis
        real(dp), allocatable :: frequencies(:) !! FFT frequency bins
        real(dp), allocatable :: times(:) !! time slice centers
        integer :: p0 = 0 !! first zero-based slice index
        integer :: status = short_fft_invalid_input !! zero on success
    end type short_time_fft_result_3d

    type, public :: short_time_fft_result_4d
        complex(dp), allocatable :: values(:,:,:,:) !! frequency and time with two batch axes
        real(dp), allocatable :: frequencies(:) !! FFT frequency bins
        real(dp), allocatable :: times(:) !! time slice centers
        integer :: p0 = 0 !! first zero-based slice index
        integer :: status = short_fft_invalid_input !! zero on success
    end type short_time_fft_result_4d

    type, public :: short_time_spectrogram_result_3d
        real(dp), allocatable :: power(:,:,:) !! multidimensional squared-magnitude STFT
        real(dp), allocatable :: frequencies(:) !! frequency bins
        real(dp), allocatable :: times(:) !! time centers
        integer :: p0 = 0 !! first zero-based slice
        integer :: status = short_fft_invalid_input !! zero on success
    end type short_time_spectrogram_result_3d

    type, public :: short_time_spectrogram_result_4d
        real(dp), allocatable :: power(:,:,:,:) !! multidimensional squared-magnitude STFT
        real(dp), allocatable :: frequencies(:) !! frequency bins
        real(dp), allocatable :: times(:) !! time centers
        integer :: p0 = 0 !! first zero-based slice
        integer :: status = short_fft_invalid_input !! zero on success
    end type short_time_spectrogram_result_4d

    ! Rank-4 signals gain a frequency axis at the selected signal axis and
    ! a final time-slice axis in the resulting rank-5 spectrum.
    type, public :: short_time_fft_result_5d
        complex(dp), allocatable :: values(:,:,:,:,:) !! STFT coefficients
        real(dp), allocatable :: frequencies(:) !! frequency bins
        real(dp), allocatable :: times(:) !! slice centers
        integer :: p0 = 0 !! first zero-based time slice
        integer :: status = short_fft_invalid_input !! zero on success
    end type short_time_fft_result_5d

    type, public :: short_time_spectrogram_result_5d
        real(dp), allocatable :: power(:,:,:,:,:) !! squared-magnitude STFT
        real(dp), allocatable :: frequencies(:) !! frequency bins
        real(dp), allocatable :: times(:) !! slice centers
        integer :: p0 = 0 !! first zero-based time slice
        integer :: status = short_fft_invalid_input !! zero on success
    end type short_time_spectrogram_result_5d

    ! Rank-5 signals gain a frequency axis at the selected signal axis and
    ! a final time-slice axis in the resulting rank-6 spectrum.
    type, public :: short_time_fft_result_6d
        complex(dp), allocatable :: values(:,:,:,:,:,:) !! STFT coefficients
        real(dp), allocatable :: frequencies(:) !! frequency bins
        real(dp), allocatable :: times(:) !! slice centers
        integer :: p0 = 0 !! first zero-based time slice
        integer :: status = short_fft_invalid_input !! zero on success
    end type short_time_fft_result_6d

    type, public :: short_time_spectrogram_result_6d
        real(dp), allocatable :: power(:,:,:,:,:,:) !! squared-magnitude STFT
        real(dp), allocatable :: frequencies(:) !! frequency bins
        real(dp), allocatable :: times(:) !! slice centers
        integer :: p0 = 0 !! first zero-based time slice
        integer :: status = short_fft_invalid_input !! zero on success
    end type short_time_spectrogram_result_6d

    ! Rank-6 signals yield rank-7 STFT coefficients.
    type, public :: short_time_fft_result_7d
        complex(dp), allocatable :: values(:,:,:,:,:,:,:) !! STFT coefficients
        real(dp), allocatable :: frequencies(:) !! frequency bins
        real(dp), allocatable :: times(:) !! time slice centers
        integer :: p0 = 0 !! first zero-based time-slice index
        integer :: status = short_fft_invalid_input !! operation status
    end type short_time_fft_result_7d

    type, public :: short_time_spectrogram_result_7d
        real(dp), allocatable :: power(:,:,:,:,:,:,:) !! squared magnitude of the STFT
        real(dp), allocatable :: frequencies(:) !! frequency bins
        real(dp), allocatable :: times(:) !! time slice centers
        integer :: p0 = 0 !! first zero-based time-slice index
        integer :: status = short_fft_invalid_input !! operation status
    end type short_time_spectrogram_result_7d

    ! Rank-7 signals yield rank-8 STFT coefficients.
    type, public :: short_time_fft_result_8d
        complex(dp), allocatable :: values(:,:,:,:,:,:,:,:) !! STFT coefficients
        real(dp), allocatable :: frequencies(:) !! frequency bins
        real(dp), allocatable :: times(:) !! time slice centers
        integer :: p0 = 0 !! first zero-based time-slice index
        integer :: status = short_fft_invalid_input !! operation status
    end type short_time_fft_result_8d

    type, public :: short_time_spectrogram_result_8d
        real(dp), allocatable :: power(:,:,:,:,:,:,:,:) !! squared magnitude of the STFT
        real(dp), allocatable :: frequencies(:) !! frequency bins
        real(dp), allocatable :: times(:) !! time slice centers
        integer :: p0 = 0 !! first zero-based time-slice index
        integer :: status = short_fft_invalid_input !! operation status
    end type short_time_spectrogram_result_8d

    ! Rank-8 signals yield rank-9 STFT coefficients.
    type, public :: short_time_fft_result_9d
        complex(dp), allocatable :: values(:,:,:,:,:,:,:,:,:) !! STFT coefficients
        real(dp), allocatable :: frequencies(:) !! frequency bins
        real(dp), allocatable :: times(:) !! time slice centers
        integer :: p0 = 0 !! first zero-based time-slice index
        integer :: status = short_fft_invalid_input !! operation status
    end type short_time_fft_result_9d

    type, public :: short_time_spectrogram_result_9d
        real(dp), allocatable :: power(:,:,:,:,:,:,:,:,:) !! squared magnitude of the STFT
        real(dp), allocatable :: frequencies(:) !! frequency bins
        real(dp), allocatable :: times(:) !! time slice centers
        integer :: p0 = 0 !! first zero-based time-slice index
        integer :: status = short_fft_invalid_input !! operation status
    end type short_time_spectrogram_result_9d

    ! Rank-9 signals yield rank-10 STFT coefficients.
    type, public :: short_time_fft_result_10d
        complex(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:) !! STFT coefficients
        real(dp), allocatable :: frequencies(:) !! frequency bins
        real(dp), allocatable :: times(:) !! time slice centers
        integer :: p0 = 0 !! first zero-based time-slice index
        integer :: status = short_fft_invalid_input !! operation status
    end type short_time_fft_result_10d

    type, public :: short_time_spectrogram_result_10d
        real(dp), allocatable :: power(:,:,:,:,:,:,:,:,:,:) !! squared magnitude of the STFT
        real(dp), allocatable :: frequencies(:) !! frequency bins
        real(dp), allocatable :: times(:) !! time slice centers
        integer :: p0 = 0 !! first zero-based time-slice index
        integer :: status = short_fft_invalid_input !! operation status
    end type short_time_spectrogram_result_10d

    ! Rank-10 signals yield rank-11 STFT coefficients.
    type, public :: short_time_fft_result_11d
        complex(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:,:) !! STFT coefficients
        real(dp), allocatable :: frequencies(:) !! frequency bins
        real(dp), allocatable :: times(:) !! time slice centers
        integer :: p0 = 0 !! first zero-based time-slice index
        integer :: status = short_fft_invalid_input !! operation status
    end type short_time_fft_result_11d

    type, public :: short_time_spectrogram_result_11d
        real(dp), allocatable :: power(:,:,:,:,:,:,:,:,:,:,:) !! squared magnitude of the STFT
        real(dp), allocatable :: frequencies(:) !! frequency bins
        real(dp), allocatable :: times(:) !! time slice centers
        integer :: p0 = 0 !! first zero-based time-slice index
        integer :: status = short_fft_invalid_input !! operation status
    end type short_time_spectrogram_result_11d

    ! Rank-11 signals yield rank-12 STFT coefficients.
    type, public :: short_time_fft_result_12d
        complex(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:,:,:) !! STFT coefficients
        real(dp), allocatable :: frequencies(:) !! frequency bins
        real(dp), allocatable :: times(:) !! time slice centers
        integer :: p0 = 0 !! first zero-based time-slice index
        integer :: status = short_fft_invalid_input !! operation status
    end type short_time_fft_result_12d

    type, public :: short_time_spectrogram_result_12d
        real(dp), allocatable :: power(:,:,:,:,:,:,:,:,:,:,:,:) !! squared magnitude of the STFT
        real(dp), allocatable :: frequencies(:) !! frequency bins
        real(dp), allocatable :: times(:) !! time slice centers
        integer :: p0 = 0 !! first zero-based time-slice index
        integer :: status = short_fft_invalid_input !! operation status
    end type short_time_spectrogram_result_12d

    ! Rank-12 signals yield rank-13 STFT coefficients.
    type, public :: short_time_fft_result_13d
        complex(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:,:,:,:) !! STFT coefficients
        real(dp), allocatable :: frequencies(:) !! frequency bins
        real(dp), allocatable :: times(:) !! time slice centers
        integer :: p0 = 0 !! first zero-based time-slice index
        integer :: status = short_fft_invalid_input !! operation status
    end type short_time_fft_result_13d

    type, public :: short_time_spectrogram_result_13d
        real(dp), allocatable :: power(:,:,:,:,:,:,:,:,:,:,:,:,:) !! squared magnitude of the STFT
        real(dp), allocatable :: frequencies(:) !! frequency bins
        real(dp), allocatable :: times(:) !! time slice centers
        integer :: p0 = 0 !! first zero-based time-slice index
        integer :: status = short_fft_invalid_input !! operation status
    end type short_time_spectrogram_result_13d

    ! Rank-13 signals yield rank-14 STFT coefficients.
    type, public :: short_time_fft_result_14d
        complex(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:,:,:,:,:) !! STFT coefficients
        real(dp), allocatable :: frequencies(:) !! frequency bins
        real(dp), allocatable :: times(:) !! time slice centers
        integer :: p0 = 0 !! first zero-based time-slice index
        integer :: status = short_fft_invalid_input !! operation status
    end type short_time_fft_result_14d

    type, public :: short_time_spectrogram_result_14d
        real(dp), allocatable :: power(:,:,:,:,:,:,:,:,:,:,:,:,:,:) !! squared magnitude of the STFT
        real(dp), allocatable :: frequencies(:) !! frequency bins
        real(dp), allocatable :: times(:) !! time slice centers
        integer :: p0 = 0 !! first zero-based time-slice index
        integer :: status = short_fft_invalid_input !! operation status
    end type short_time_spectrogram_result_14d

    ! Rank-14 signals yield rank-15 STFT coefficients.
    type, public :: short_time_fft_result_15d
        complex(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:) !! STFT coefficients
        real(dp), allocatable :: frequencies(:) !! frequency bins
        real(dp), allocatable :: times(:) !! time slice centers
        integer :: p0 = 0 !! first zero-based time-slice index
        integer :: status = short_fft_invalid_input !! operation status
    end type short_time_fft_result_15d

    type, public :: short_time_spectrogram_result_15d
        real(dp), allocatable :: power(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:) !! squared magnitude of the STFT
        real(dp), allocatable :: frequencies(:) !! frequency bins
        real(dp), allocatable :: times(:) !! time slice centers
        integer :: p0 = 0 !! first zero-based time-slice index
        integer :: status = short_fft_invalid_input !! operation status
    end type short_time_spectrogram_result_15d

    type, public :: short_time_fft
        complex(dp), allocatable :: win(:) !! scaled analysis window
        complex(dp), allocatable :: dual_win(:) !! canonical synthesis window
        integer :: hop = 0 !! slice increment
        integer :: mfft = 0 !! transform length, at least window length
        real(dp) :: fs = 1.0_dp !! sampling frequency
        character(12) :: fft_mode = 'onesided' !! onesided, onesided2X, twosided, centered
        character(12) :: scaling = 'none' !! none, magnitude, psd
        integer :: phase_shift = 0 !! linear phase in samples
        integer :: p_min = 0 !! earliest window slice overlapping input start
        integer :: status = short_fft_invalid_input !! initialization status
    contains
        procedure, private :: initialize_real, initialize_complex
        generic :: initialize => initialize_real, initialize_complex
        procedure :: initialize_window
        procedure, private :: initialize_from_dual_real, initialize_from_dual_complex
        generic :: initialize_from_dual => initialize_from_dual_real, initialize_from_dual_complex
        procedure, private :: equal_dual_real, equal_dual_complex
        generic :: initialize_equal_dual => equal_dual_real, equal_dual_complex
        procedure :: p_max
        procedure :: p_num
        procedure :: p_range
        procedure :: k_min
        procedure :: k_max
        procedure :: lower_border_end
        procedure :: upper_border_begin
        procedure :: nearest_k_p
        procedure :: time_points
        procedure :: frequency_points
        procedure :: frequency_count
        procedure :: delta_t
        procedure :: delta_f
        procedure :: invertible
        procedure :: extent
        procedure, private :: stft_real, stft_complex
        procedure, private :: stft_nd2_real, stft_nd2_complex
        procedure, private :: stft_nd3_real, stft_nd3_complex
        procedure, private :: stft_nd4_real, stft_nd4_complex
        procedure, private :: stft_nd5_real, stft_nd5_complex
        procedure, private :: stft_nd6_real, stft_nd6_complex
        procedure, private :: stft_nd7_real, stft_nd7_complex
        procedure, private :: stft_nd8_real, stft_nd8_complex
        procedure, private :: stft_nd9_real, stft_nd9_complex
        procedure, private :: stft_nd10_real, stft_nd10_complex
        procedure, private :: stft_nd11_real, stft_nd11_complex
        procedure, private :: stft_nd12_real, stft_nd12_complex
        procedure, private :: stft_nd13_real, stft_nd13_complex
        procedure, private :: stft_nd14_real, stft_nd14_complex
        generic :: stft => stft_real, stft_complex, &
            stft_nd2_real, stft_nd2_complex, &
            stft_nd3_real, stft_nd3_complex, &
            stft_nd4_real, stft_nd4_complex, &
            stft_nd5_real, stft_nd5_complex, &
            stft_nd6_real, stft_nd6_complex, &
            stft_nd7_real, stft_nd7_complex, &
            stft_nd8_real, stft_nd8_complex, &
            stft_nd9_real, stft_nd9_complex, &
            stft_nd10_real, stft_nd10_complex, &
            stft_nd11_real, stft_nd11_complex, &
            stft_nd12_real, stft_nd12_complex, &
            stft_nd13_real, stft_nd13_complex, &
            stft_nd14_real, stft_nd14_complex
        procedure, private :: stft_detrend_real, stft_detrend_complex
        procedure, private :: stft_detrend_callable_real, stft_detrend_callable_complex
        procedure, private :: stft_detrend_callable_nd2_real, stft_detrend_callable_nd2_complex
        procedure, private :: stft_detrend_callable_nd3_real, stft_detrend_callable_nd3_complex
        procedure, private :: stft_detrend_nd2_real, stft_detrend_nd2_complex
        procedure, private :: stft_detrend_nd3_real, stft_detrend_nd3_complex
        procedure, private :: stft_detrend_nd4_real, stft_detrend_nd4_complex
        procedure, private :: stft_detrend_nd5_real, stft_detrend_nd5_complex
        procedure, private :: stft_detrend_nd6_real, stft_detrend_nd6_complex
        procedure, private :: stft_detrend_nd7_real, stft_detrend_nd7_complex
        procedure, private :: stft_detrend_nd8_real, stft_detrend_nd8_complex
        procedure, private :: stft_detrend_nd9_real, stft_detrend_nd9_complex
        procedure, private :: stft_detrend_nd10_real, stft_detrend_nd10_complex
        procedure, private :: stft_detrend_nd11_real, stft_detrend_nd11_complex
        procedure, private :: stft_detrend_nd12_real, stft_detrend_nd12_complex
        procedure, private :: stft_detrend_nd13_real, stft_detrend_nd13_complex
        procedure, private :: stft_detrend_nd14_real, stft_detrend_nd14_complex
        generic :: stft_detrend => stft_detrend_real, stft_detrend_complex, &
            stft_detrend_callable_real, stft_detrend_callable_complex, &
            stft_detrend_callable_nd2_real, stft_detrend_callable_nd2_complex, &
            stft_detrend_callable_nd3_real, stft_detrend_callable_nd3_complex, &
            stft_detrend_nd2_real, stft_detrend_nd2_complex, &
            stft_detrend_nd3_real, stft_detrend_nd3_complex, &
            stft_detrend_nd4_real, stft_detrend_nd4_complex, &
            stft_detrend_nd5_real, stft_detrend_nd5_complex, &
            stft_detrend_nd6_real, stft_detrend_nd6_complex, &
            stft_detrend_nd7_real, stft_detrend_nd7_complex, &
            stft_detrend_nd8_real, stft_detrend_nd8_complex, &
            stft_detrend_nd9_real, stft_detrend_nd9_complex, &
            stft_detrend_nd10_real, stft_detrend_nd10_complex, &
            stft_detrend_nd11_real, stft_detrend_nd11_complex, &
            stft_detrend_nd12_real, stft_detrend_nd12_complex, &
            stft_detrend_nd13_real, stft_detrend_nd13_complex, &
            stft_detrend_nd14_real, stft_detrend_nd14_complex
        procedure, private :: istft_real, istft_complex
        procedure, private :: istft_nd3_real, istft_nd3_complex
        procedure, private :: istft_nd4_real, istft_nd4_complex
        procedure, private :: istft_nd5_real, istft_nd5_complex
        procedure, private :: istft_nd6_real, istft_nd6_complex
        procedure, private :: istft_nd7_real, istft_nd7_complex
        procedure, private :: istft_nd8_real, istft_nd8_complex
        procedure, private :: istft_nd9_real, istft_nd9_complex
        procedure, private :: istft_nd10_real, istft_nd10_complex
        procedure, private :: istft_nd11_real, istft_nd11_complex
        procedure, private :: istft_nd12_real, istft_nd12_complex
        procedure, private :: istft_nd13_real, istft_nd13_complex
        procedure, private :: istft_nd14_real, istft_nd14_complex
        procedure, private :: istft_nd15_real, istft_nd15_complex
        generic :: istft => istft_real, istft_complex, &
            istft_nd3_real, istft_nd3_complex, &
            istft_nd4_real, istft_nd4_complex, &
            istft_nd5_real, istft_nd5_complex, &
            istft_nd6_real, istft_nd6_complex, &
            istft_nd7_real, istft_nd7_complex, &
            istft_nd8_real, istft_nd8_complex, &
            istft_nd9_real, istft_nd9_complex, &
            istft_nd10_real, istft_nd10_complex, &
            istft_nd11_real, istft_nd11_complex, &
            istft_nd12_real, istft_nd12_complex, &
            istft_nd13_real, istft_nd13_complex, &
            istft_nd14_real, istft_nd14_complex, &
            istft_nd15_real, istft_nd15_complex
        procedure, private :: spectrogram_real, spectrogram_complex
        procedure, private :: spectrogram_nd2_real, spectrogram_nd2_complex
        procedure, private :: spectrogram_nd3_real, spectrogram_nd3_complex
        procedure, private :: spectrogram_nd4_real, spectrogram_nd4_complex
        procedure, private :: spectrogram_nd5_real, spectrogram_nd5_complex
        procedure, private :: spectrogram_nd6_real, spectrogram_nd6_complex
        procedure, private :: spectrogram_nd7_real, spectrogram_nd7_complex
        procedure, private :: spectrogram_nd8_real, spectrogram_nd8_complex
        procedure, private :: spectrogram_nd9_real, spectrogram_nd9_complex
        procedure, private :: spectrogram_nd10_real, spectrogram_nd10_complex
        procedure, private :: spectrogram_nd11_real, spectrogram_nd11_complex
        procedure, private :: spectrogram_nd12_real, spectrogram_nd12_complex
        procedure, private :: spectrogram_nd13_real, spectrogram_nd13_complex
        procedure, private :: spectrogram_nd14_real, spectrogram_nd14_complex
        generic :: spectrogram => spectrogram_real, spectrogram_complex, &
            spectrogram_nd2_real, spectrogram_nd2_complex, &
            spectrogram_nd3_real, spectrogram_nd3_complex, &
            spectrogram_nd4_real, spectrogram_nd4_complex, &
            spectrogram_nd5_real, spectrogram_nd5_complex, &
            spectrogram_nd6_real, spectrogram_nd6_complex, &
            spectrogram_nd7_real, spectrogram_nd7_complex, &
            spectrogram_nd8_real, spectrogram_nd8_complex, &
            spectrogram_nd9_real, spectrogram_nd9_complex, &
            spectrogram_nd10_real, spectrogram_nd10_complex, &
            spectrogram_nd11_real, spectrogram_nd11_complex, &
            spectrogram_nd12_real, spectrogram_nd12_complex, &
            spectrogram_nd13_real, spectrogram_nd13_complex, &
            spectrogram_nd14_real, spectrogram_nd14_complex
        procedure, private :: cross_spectrogram_real, cross_spectrogram_complex
        procedure, private :: cross_spectrogram_real_complex, cross_spectrogram_complex_real
        procedure, private :: cross_spectrogram_nd2_real, cross_spectrogram_nd2_complex
        procedure, private :: cross_spectrogram_nd3_real, cross_spectrogram_nd3_complex
        procedure, private :: cross_spectrogram_nd4_real, cross_spectrogram_nd4_complex
        procedure, private :: cross_spectrogram_nd5_real, cross_spectrogram_nd5_complex
        procedure, private :: cross_spectrogram_nd6_real, cross_spectrogram_nd6_complex
        procedure, private :: cross_spectrogram_nd7_real, cross_spectrogram_nd7_complex
        procedure, private :: cross_spectrogram_nd8_real, cross_spectrogram_nd8_complex
        procedure, private :: cross_spectrogram_nd9_real, cross_spectrogram_nd9_complex
        procedure, private :: cross_spectrogram_nd10_real, cross_spectrogram_nd10_complex
        procedure, private :: cross_spectrogram_nd11_real, cross_spectrogram_nd11_complex
        procedure, private :: cross_spectrogram_nd12_real, cross_spectrogram_nd12_complex
        procedure, private :: cross_spectrogram_nd13_real, cross_spectrogram_nd13_complex
        procedure, private :: cross_spectrogram_nd14_real, cross_spectrogram_nd14_complex
        generic :: cross_spectrogram => cross_spectrogram_real, cross_spectrogram_complex, &
            cross_spectrogram_real_complex, cross_spectrogram_complex_real, &
            cross_spectrogram_nd2_real, cross_spectrogram_nd2_complex, &
            cross_spectrogram_nd3_real, cross_spectrogram_nd3_complex, &
            cross_spectrogram_nd4_real, cross_spectrogram_nd4_complex, &
            cross_spectrogram_nd5_real, cross_spectrogram_nd5_complex, &
            cross_spectrogram_nd6_real, cross_spectrogram_nd6_complex, &
            cross_spectrogram_nd7_real, cross_spectrogram_nd7_complex, &
            cross_spectrogram_nd8_real, cross_spectrogram_nd8_complex, &
            cross_spectrogram_nd9_real, cross_spectrogram_nd9_complex, &
            cross_spectrogram_nd10_real, cross_spectrogram_nd10_complex, &
            cross_spectrogram_nd11_real, cross_spectrogram_nd11_complex, &
            cross_spectrogram_nd12_real, cross_spectrogram_nd12_complex, &
            cross_spectrogram_nd13_real, cross_spectrogram_nd13_complex, &
            cross_spectrogram_nd14_real, cross_spectrogram_nd14_complex
    end type short_time_fft

contains

    subroutine initialize_window(self, name, nperseg, hop, fs, mfft, fft_mode, scale_to, phase_shift, &
            parameter, parameter2, symmetric)
        class(short_time_fft), intent(inout) :: self !! analysis plan
        character(*), intent(in) :: name !! window name accepted by get_window
        integer, intent(in) :: nperseg !! window length
        integer, intent(in) :: hop !! samples per step
        real(dp), intent(in) :: fs !! sampling rate
        integer, intent(in), optional :: mfft !! transform length
        character(*), intent(in), optional :: fft_mode !! FFT spectrum layout
        character(*), intent(in), optional :: scale_to !! none, psd, magnitude
        integer, intent(in), optional :: phase_shift !! phase offset
        real(dp), intent(in), optional :: parameter !! first window shape parameter
        real(dp), intent(in), optional :: parameter2 !! second window shape parameter
        logical, intent(in), optional :: symmetric !! symmetric, rather than periodic
        real(dp), allocatable :: w(:)
        integer :: ierr
        logical :: periodic

        self%status = short_fft_invalid_input
        if (nperseg < 1) return
        periodic = .true.
        if (present(symmetric)) periodic = .not. symmetric
        w = get_window(name, nperseg, fftbins=periodic, status=ierr, &
            parameter=parameter, parameter2=parameter2)
        if (ierr /= spectral_status_success) return
        call self%initialize_real(w, hop, fs, mfft, fft_mode, scale_to, phase_shift)
    end subroutine initialize_window

    subroutine equal_dual_real(self, desired_window, hop, fs, mfft, fft_mode, scale_to, phase_shift)
        class(short_time_fft), intent(inout) :: self !! STFT plan
        real(dp), intent(in) :: desired_window(:) !! approximate common window and dual
        integer, intent(in) :: hop !! sample step
        real(dp), intent(in) :: fs !! sample rate
        integer, intent(in), optional :: mfft !! FFT size
        character(*), intent(in), optional :: fft_mode !! spectral FFT mode
        character(*), intent(in), optional :: scale_to !! none, magnitude, psd, unitary
        integer, intent(in), optional :: phase_shift !! phase convention
        call self%equal_dual_complex(cmplx(desired_window,0.0_dp,dp), hop, fs, mfft, &
            fft_mode, scale_to, phase_shift)
    end subroutine equal_dual_real

    subroutine equal_dual_complex(self, desired_window, hop, fs, mfft, fft_mode, scale_to, phase_shift)
        class(short_time_fft), intent(inout) :: self !! STFT plan
        complex(dp), intent(in) :: desired_window(:) !! target equal-dual window
        integer, intent(in) :: hop !! sample step
        real(dp), intent(in) :: fs !! sampling frequency
        integer, intent(in), optional :: mfft !! FFT length
        character(*), intent(in), optional :: fft_mode !! spectral layout
        character(*), intent(in), optional :: scale_to !! none, magnitude, psd, unitary
        integer, intent(in), optional :: phase_shift !! samples
        complex(dp), allocatable :: tight(:)
        real(dp), allocatable :: energy(:)
        real(dp) :: factor
        integer :: j, n, nf
        character(:), allocatable :: scale

        self%status = short_fft_invalid_input
        n = size(desired_window)
        nf = n
        if (present(mfft)) nf = mfft
        if (hop < 1 .or. n < 1 .or. hop > n .or. nf < n) return
        if (.not. all(ieee_is_finite(real(desired_window,dp))) .or. &
            .not. all(ieee_is_finite(aimag(desired_window)))) return
        scale = 'none'
        if (present(scale_to)) scale = trim(scale_to)
        select case (scale)
        case ('none','magnitude','psd','unitary')
        case default
            return
        end select
        allocate(energy(hop),tight(n))
        energy = 0.0_dp
        do j=1,n
            energy(mod(j-1,hop)+1) = energy(mod(j-1,hop)+1)+abs(desired_window(j))**2
        end do
        if (any(energy <= tiny(1.0_dp))) then
            self%status = short_fft_noninvertible
            return
        end if
        do j=1,n
            tight(j) = desired_window(j)/sqrt(energy(mod(j-1,hop)+1))
        end do
        if (scale == 'unitary') then
            factor = sqrt(real(nf,dp))
            call self%initialize_complex(tight/factor,hop,fs,mfft,fft_mode, &
                phase_shift=phase_shift,dual_window=tight*factor)
            if (self%status == short_fft_success) self%scaling = 'unitary'
        else
            call self%initialize_complex(tight,hop,fs,mfft,fft_mode,scale_to,phase_shift)
        end if
    end subroutine equal_dual_complex

    subroutine initialize_from_dual_real(self, dual_window, hop, fs, mfft, fft_mode, phase_shift)
        class(short_time_fft), intent(inout) :: self !! transform plan
        real(dp), intent(in) :: dual_window(:) !! synthesis window
        integer, intent(in) :: hop !! samples between time slices
        real(dp), intent(in) :: fs !! sample rate
        integer, intent(in), optional :: mfft !! FFT size
        character(*), intent(in), optional :: fft_mode !! one-/two-sided FFT mode
        integer, intent(in), optional :: phase_shift !! integer linear phase shift
        complex(dp), allocatable :: win(:)
        integer :: n, j
        real(dp), allocatable :: energy(:)

        n = size(dual_window)
        self%status = short_fft_invalid_input
        if (hop < 1 .or. n < 1 .or. hop > n) return
        allocate(energy(hop),win(n))
        energy = 0.0_dp
        do j=1,n
            energy(mod(j-1,hop)+1) = energy(mod(j-1,hop)+1)+dual_window(j)**2
        end do
        if (any(energy <= tiny(1.0_dp))) then
            self%status = short_fft_noninvertible
            return
        end if
        do j=1,n
            win(j) = dual_window(j)/energy(mod(j-1,hop)+1)
        end do
        call self%initialize_complex(win,hop,fs,mfft,fft_mode,phase_shift=phase_shift, &
            dual_window=cmplx(dual_window,0.0_dp,dp))
    end subroutine initialize_from_dual_real

    subroutine initialize_from_dual_complex(self, dual_window, hop, fs, mfft, fft_mode, phase_shift)
        class(short_time_fft), intent(inout) :: self !! transform plan
        complex(dp), intent(in) :: dual_window(:) !! synthesis window
        integer, intent(in) :: hop !! time-slice increment
        real(dp), intent(in) :: fs !! sampling rate
        integer, intent(in), optional :: mfft !! FFT size
        character(*), intent(in), optional :: fft_mode !! FFT mode
        integer, intent(in), optional :: phase_shift !! phase shift
        complex(dp), allocatable :: win(:)
        real(dp), allocatable :: energy(:)
        integer :: n, j

        n = size(dual_window)
        self%status = short_fft_invalid_input
        if (hop < 1 .or. n < 1 .or. hop > n) return
        allocate(energy(hop),win(n))
        energy = 0.0_dp
        do j=1,n
            energy(mod(j-1,hop)+1) = energy(mod(j-1,hop)+1)+abs(dual_window(j))**2
        end do
        if (any(energy <= tiny(1.0_dp))) then
            self%status = short_fft_noninvertible
            return
        end if
        do j=1,n
            win(j) = dual_window(j)/energy(mod(j-1,hop)+1)
        end do
        call self%initialize_complex(win,hop,fs,mfft,fft_mode,phase_shift=phase_shift, &
            dual_window=dual_window)
    end subroutine initialize_from_dual_complex

    subroutine initialize_real(self, window, hop, fs, mfft, fft_mode, scale_to, phase_shift, dual_window)
        class(short_time_fft), intent(inout) :: self !! FFT configuration
        real(dp), intent(in) :: window(:) !! analysis window
        integer, intent(in) :: hop !! sample step
        real(dp), intent(in) :: fs !! sample rate
        integer, intent(in), optional :: mfft !! transform length
        character(*), intent(in), optional :: fft_mode !! FFT output type
        character(*), intent(in), optional :: scale_to !! magnitude, psd, none
        integer, intent(in), optional :: phase_shift !! zero by default
        real(dp), intent(in), optional :: dual_window(:) !! user-supplied synthesis window
        complex(dp), allocatable :: dual(:)

        if (present(dual_window)) dual = cmplx(dual_window, 0.0_dp, dp)
        if (present(dual_window)) then
            call self%initialize_complex(cmplx(window, 0.0_dp, dp), hop, fs, mfft, fft_mode, &
                scale_to, phase_shift, dual)
        else
            call self%initialize_complex(cmplx(window, 0.0_dp, dp), hop, fs, mfft, fft_mode, &
                scale_to, phase_shift)
        end if
    end subroutine initialize_real

    subroutine initialize_complex(self, window, hop, fs, mfft, fft_mode, scale_to, phase_shift, dual_window)
        class(short_time_fft), intent(inout) :: self !! FFT configuration
        complex(dp), intent(in) :: window(:) !! real or complex analysis window
        integer, intent(in) :: hop !! hop in samples
        real(dp), intent(in) :: fs !! sample rate
        integer, intent(in), optional :: mfft !! transform length
        character(*), intent(in), optional :: fft_mode !! FFT output layout
        character(*), intent(in), optional :: scale_to !! magnitude, psd, none
        integer, intent(in), optional :: phase_shift !! samples to shift FFT phase
        complex(dp), intent(in), optional :: dual_window(:) !! explicit synthesis window
        integer :: n, j, nf, phase
        real(dp) :: factor
        real(dp), allocatable :: energy(:)
        character(:), allocatable :: mode, scale

        self%status = short_fft_invalid_input
        if (allocated(self%win)) deallocate(self%win)
        if (allocated(self%dual_win)) deallocate(self%dual_win)
        n = size(window)
        if (n < 1 .or. hop < 1 .or. hop > n) return
        if (.not. ieee_is_finite(fs) .or. fs <= 0.0_dp) return
        if (.not. all(ieee_is_finite(real(window,dp))) .or. .not. all(ieee_is_finite(aimag(window)))) return
        nf = n
        if (present(mfft)) nf = mfft
        if (nf < n) return
        mode = 'onesided'
        if (present(fft_mode)) mode = trim(fft_mode)
        select case (mode)
        case ('onesided', 'onesided2X')
            if (any(aimag(window) /= 0.0_dp)) return
        case ('twosided', 'centered')
        case default
            return
        end select
        scale = 'none'
        if (present(scale_to)) scale = trim(scale_to)
        if (mode == 'onesided2X' .and. scale == 'none') return
        select case (scale)
        case ('none')
            factor = 1.0_dp
        case ('magnitude')
            factor = abs(sum(window))
            if (factor <= tiny(factor)) return
            factor = 1.0_dp / factor
        case ('psd')
            factor = sqrt(fs*sum(abs(window)**2))
            if (.not. ieee_is_finite(factor) .or. factor <= tiny(factor)) return
            factor = 1.0_dp / factor
        case default
            return
        end select
        phase = 0
        if (present(phase_shift)) phase = phase_shift
        if (abs(real(phase,dp)) >= real(nf,dp)) return
        allocate(energy(hop))
        energy = 0.0_dp
        do j = 1,n
            energy(mod(j-1,hop)+1) = energy(mod(j-1,hop)+1) + abs(window(j))**2
        end do
        if (any(energy <= tiny(1.0_dp))) then
            self%status = short_fft_noninvertible
            return
        end if
        self%win = factor*window
        self%hop = hop
        self%fs = fs
        self%mfft = nf
        self%fft_mode = mode
        self%scaling = scale
        self%phase_shift = phase
        allocate(self%dual_win(n))
        if (present(dual_window)) then
            if (size(dual_window) /= n) then
                self%status = short_fft_invalid_input
                return
            end if
            if (.not. all(ieee_is_finite(real(dual_window,dp))) .or. &
                .not. all(ieee_is_finite(aimag(dual_window)))) return
            ! A provided dual must satisfy the sum of win*conjg(dual) over each residue.
            do j = 1,hop
                if (abs(sum(conjg(self%win(j:n:hop))*dual_window(j:n:hop)) - (1.0_dp,0.0_dp)) > &
                    1.0e-10_dp) return
            end do
            self%dual_win = dual_window
        else
            do j = 1,n
                self%dual_win(j) = window(j)/(energy(mod(j-1,hop)+1)*factor)
            end do
        end if
        self%p_min = 0
        do j = -1,-n/hop-2,-1
            if ((j*hop-n/2+n) <= 0) exit
            if (all(abs(self%win(max(1,1-j*hop+n/2):n)) == 0.0_dp)) exit
            self%p_min = j
        end do
        self%status = short_fft_success
    end subroutine initialize_complex

    integer function p_max(self, n) result(p1)
        class(short_time_fft), intent(in) :: self !! configuration
        integer, intent(in) :: n !! input length
        integer :: p, next_first, last

        p1 = 0
        if (self%status /= short_fft_success) return
        if (n < size(self%win)-size(self%win)/2) return
        do p = n/self%hop, n/self%hop+size(self%win)/self%hop+3
            next_first = (p+1)*self%hop-size(self%win)/2
            if (next_first >= n) then
                p1 = p+1
                return
            end if
            last = n-next_first
            if (last <= 0) then
                p1 = p+1
                return
            end if
            if (last < size(self%win)) then
                if (all(abs(self%win(1:last)) == 0.0_dp)) then
                    p1 = p+1
                    return
                end if
            end if
        end do
    end function p_max

    integer function p_num(self, n) result(num)
        class(short_time_fft), intent(in) :: self !! plan
        integer, intent(in) :: n !! signal length
        num = 0
        if (self%status /= short_fft_success) return
        if (n < size(self%win)-size(self%win)/2) return
        num = self%p_max(n)-self%p_min
    end function p_num

    ! Coordinate and support helpers corresponding to SciPy's ShortTimeFFT.
    ! Sample and slice indices are zero-based; the returned Fortran arrays
    ! are one-based as usual.  Invalid plans and ranges return an error code.
    subroutine p_range(self, n, first, last, status, p0, p1)
        class(short_time_fft), intent(in) :: self !! STFT plan
        integer, intent(in) :: n !! signal length
        integer, intent(out) :: first !! inclusive first slice
        integer, intent(out) :: last !! exclusive last slice
        integer, intent(out) :: status !! short_fft_success or error
        integer, intent(in), optional :: p0 !! optional inclusive first slice
        integer, intent(in), optional :: p1 !! optional exclusive last slice
        integer :: pm

        first = 0
        last = 0
        status = short_fft_invalid_input
        if (self%status /= short_fft_success) return
        if (n < size(self%win) - size(self%win)/2) return
        pm = self%p_max(n)
        first = self%p_min
        last = pm
        if (present(p0)) first = p0
        if (present(p1)) last = p1
        if (first < self%p_min .or. first >= last .or. last > pm) return
        status = short_fft_success
    end subroutine p_range

    integer function k_min(self) result(k)
        class(short_time_fft), intent(in) :: self !! STFT plan

        k = 0
        if (self%status /= short_fft_success) return
        k = self%p_min*self%hop - size(self%win)/2
    end function k_min

    integer function k_max(self, n) result(k)
        class(short_time_fft), intent(in) :: self !! STFT plan
        integer, intent(in) :: n !! signal length
        integer :: nw, pm

        k = 0
        if (self%status /= short_fft_success) return
        nw = size(self%win)
        if (n < nw-nw/2) return
        pm = self%p_max(n)
        k = (pm-1)*self%hop - nw/2 + nw
    end function k_max

    function lower_border_end(self) result(bound)
        class(short_time_fft), intent(in) :: self !! STFT plan
        integer :: bound(2) !! first sample and slice free of left padding
        integer :: nw, m0, k0, k, q, mid

        bound = 0
        if (self%status /= short_fft_success) return
        nw = size(self%win)
        mid = nw/2
        m0 = findloc(abs(self%win) > 0.0_dp, .true., dim=1)
        if (m0 == 0) return
        k0 = -mid+m0-1
        q = 0
        do k=k0,self%hop,self%hop
            if (k+self%hop >= 0) then
                bound = [k+nw,q+1]
                return
            end if
            q = q+1
        end do
        bound = [0,max(self%p_min,0)]
    end function lower_border_end

    function upper_border_begin(self, n) result(bound)
        class(short_time_fft), intent(in) :: self !! STFT plan
        integer, intent(in) :: n !! signal length
        integer :: bound(2) !! sample and slice where right padding begins
        integer :: nw, mid, q, q1, q2, k, j

        bound = 0
        if (self%status /= short_fft_success) return
        nw = size(self%win)
        mid = nw/2
        if (n < nw-mid) return
        q2 = n/self%hop+1
        q1 = max(floor(real(n-nw,dp)/real(self%hop,dp))-1,-1)
        do q=q2,q1+1,-1
            k = q*self%hop+nw-mid
            if (k <= n) then
                bound = [(q+1)*self%hop-mid,q+1]
                return
            end if
            ! Python's win[n-k:] means the final (k-n) window values.
            j = max(1,nw+n-k+1)
            if (all(abs(self%win(j:nw)) == 0.0_dp)) then
                bound = [(q+1)*self%hop-mid,q+1]
                return
            end if
        end do
    end function upper_border_begin

    integer function nearest_k_p(self, k, left) result(kp)
        class(short_time_fft), intent(in) :: self !! STFT plan
        integer, intent(in) :: k !! requested zero-based sample index
        logical, intent(in), optional :: left !! round left (default) or right
        integer :: q
        logical :: round_left

        kp = 0
        if (self%status /= short_fft_success) return
        round_left = .true.
        if (present(left)) round_left = left
        ! Fortran integer division truncates toward zero: use floor for
        ! negative input indices to reproduce Python's // operator.
        q = floor(real(k,dp)/real(self%hop,dp))
        kp = q*self%hop
        if (.not. round_left .and. kp /= k) kp = kp+self%hop
    end function nearest_k_p

    function time_points(self, n, status, p0, p1, k_offset) result(times)
        class(short_time_fft), intent(in) :: self !! STFT plan
        integer, intent(in) :: n !! signal length
        integer, intent(out) :: status !! zero if successful
        integer, intent(in), optional :: p0 !! first slice
        integer, intent(in), optional :: p1 !! last slice, exclusive
        integer, intent(in), optional :: k_offset !! input sample offset
        real(dp), allocatable :: times(:) !! STFT centers, seconds
        integer :: first, last, p, off

        allocate(times(0))
        call self%p_range(n,first,last,status,p0,p1)
        if (status /= short_fft_success) return
        off = 0
        if (present(k_offset)) off = k_offset
        times = [(real(p*self%hop+off,dp)/self%fs,p=first,last-1)]
    end function time_points

    integer function frequency_count(self) result(n)
        class(short_time_fft), intent(in) :: self !! STFT plan

        n = 0
        if (self%status /= short_fft_success) return
        n = self%mfft
        if (self%fft_mode == 'onesided' .or. self%fft_mode == 'onesided2X') n = self%mfft/2+1
    end function frequency_count

    function frequency_points(self) result(freq)
        class(short_time_fft), intent(in) :: self !! STFT plan
        real(dp), allocatable :: freq(:) !! frequencies in hertz
        integer :: nf, q, k

        allocate(freq(0))
        if (self%status /= short_fft_success) return
        nf = self%frequency_count()
        deallocate(freq)
        allocate(freq(nf))
        do q=0,nf-1
            k = q
            select case (self%fft_mode)
            case ('twosided')
                if (q >= (self%mfft+1)/2) k = q-self%mfft
            case ('centered')
                k = q-self%mfft/2
            end select
            freq(q+1) = real(k,dp)*self%fs/real(self%mfft,dp)
        end do
    end function frequency_points

    real(dp) function delta_t(self) result(dt)
        class(short_time_fft), intent(in) :: self !! STFT plan
        dt = 0.0_dp
        if (self%status == short_fft_success) dt = real(self%hop,dp)/self%fs
    end function delta_t

    real(dp) function delta_f(self) result(df)
        class(short_time_fft), intent(in) :: self !! STFT plan
        df = 0.0_dp
        if (self%status == short_fft_success) df = self%fs/real(self%mfft,dp)
    end function delta_f

    logical function invertible(self) result(ok)
        class(short_time_fft), intent(in) :: self !! STFT plan
        ok = self%status == short_fft_success .and. allocated(self%dual_win)
    end function invertible

    subroutine extent(self, n, bounds, status, axes_seq, center_bins)
        class(short_time_fft), intent(in) :: self !! STFT plan
        integer, intent(in) :: n !! number of input samples
        real(dp), intent(out) :: bounds(4) !! tf or ft plot extent
        integer, intent(out) :: status !! zero on success
        character(*), intent(in), optional :: axes_seq !! tf or ft
        logical, intent(in), optional :: center_bins !! center bins on coordinates
        real(dp) :: t0, t1, f0, f1, shift
        integer :: q0, q1, pm
        character(:), allocatable :: order

        bounds = 0.0_dp
        status = short_fft_invalid_input
        if (self%status /= short_fft_success) return
        if (n < size(self%win)-size(self%win)/2) return
        order = 'tf'
        if (present(axes_seq)) order = trim(axes_seq)
        if (order /= 'tf' .and. order /= 'ft') return
        select case (self%fft_mode)
        case ('onesided','onesided2X')
            q0 = 0
            q1 = self%frequency_count()
        case ('centered')
            q0 = -self%mfft/2
            q1 = self%mfft/2+mod(self%mfft,2)
        case default
            ! SciPy does not assign a monotone imshow extent to two-sided FFTs.
            return
        end select
        shift = 0.0_dp
        if (present(center_bins)) then
            if (center_bins) shift = 0.5_dp
        end if
        pm = self%p_max(n)
        t0 = (real(self%p_min,dp)-shift)*self%delta_t()
        t1 = (real(pm,dp)-shift)*self%delta_t()
        f0 = (real(q0,dp)-shift)*self%delta_f()
        f1 = (real(q1,dp)-shift)*self%delta_f()
        if (order == 'tf') then
            bounds = [t0,t1,f0,f1]
        else
            bounds = [f0,f1,t0,t1]
        end if
        status = short_fft_success
    end subroutine extent

    subroutine stft_real(self, x, result, p0, p1, k_offset, padding)
        class(short_time_fft), intent(in) :: self !! transform settings
        real(dp), intent(in) :: x(:) !! time-domain signal
        type(short_time_fft_result), intent(out) :: result !! STFT spectrum and axes
        integer, intent(in), optional :: p0 !! first slice, zero-based
        integer, intent(in), optional :: p1 !! exclusive final slice
        integer, intent(in), optional :: k_offset !! signal index offset
        character(*), intent(in), optional :: padding !! zeros, edge, even, odd
        call stft_kernel(self, cmplx(x,0.0_dp,dp), result, p0, p1, k_offset, padding, .false.)
    end subroutine stft_real

    subroutine stft_complex(self, x, result, p0, p1, k_offset, padding)
        class(short_time_fft), intent(in) :: self !! transform settings
        complex(dp), intent(in) :: x(:) !! time-domain samples
        type(short_time_fft_result), intent(out) :: result !! STFT spectrum and axes
        integer, intent(in), optional :: p0 !! first slice, zero-based
        integer, intent(in), optional :: p1 !! exclusive final slice
        integer, intent(in), optional :: k_offset !! sample index shift
        character(*), intent(in), optional :: padding !! zeros, edge, even, odd
        call stft_kernel(self, x, result, p0, p1, k_offset, padding, .true.)
    end subroutine stft_complex

    subroutine stft_detrend_real(self, x, result, detrend, p0, p1, k_offset, padding)
        class(short_time_fft), intent(in) :: self !! STFT configuration
        real(dp), intent(in) :: x(:) !! time signal
        type(short_time_fft_result), intent(out) :: result !! frequency by time slices
        character(*), intent(in) :: detrend !! constant, linear, or none
        integer, intent(in), optional :: p0 !! first slice index
        integer, intent(in), optional :: p1 !! last slice index, exclusive
        integer, intent(in), optional :: k_offset !! sample offset
        character(*), intent(in), optional :: padding !! boundary extension
        call stft_kernel(self,cmplx(x,0.0_dp,dp),result,p0,p1,k_offset,padding,.false.,detrend)
    end subroutine stft_detrend_real

    subroutine stft_detrend_complex(self, x, result, detrend, p0, p1, k_offset, padding)
        class(short_time_fft), intent(in) :: self !! transform settings
        complex(dp), intent(in) :: x(:) !! samples
        type(short_time_fft_result), intent(out) :: result !! transformed slices
        character(*), intent(in) :: detrend !! constant, linear, or none
        integer, intent(in), optional :: p0 !! first slice index
        integer, intent(in), optional :: p1 !! last slice index
        integer, intent(in), optional :: k_offset !! sample offset
        character(*), intent(in), optional :: padding !! boundary handling
        call stft_kernel(self,x,result,p0,p1,k_offset,padding,.true.,detrend)
    end subroutine stft_detrend_complex

    subroutine stft_detrend_callable_real(self, x, result, detrend, p0, p1, k_offset, padding)
        class(short_time_fft), intent(in) :: self !! initialized ShortTimeFFT plan
        real(dp), intent(in) :: x(:) !! time-domain real samples
        type(short_time_fft_result), intent(out) :: result !! complex STFT and coordinates
        procedure(stft_real_detrend_callback) :: detrend !! callable on each real slice
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: p1 !! last time slice, exclusive
        integer, intent(in), optional :: k_offset !! sample offset
        character(*), intent(in), optional :: padding !! zeros, edge, even or odd

        call stft_kernel(self, cmplx(x, 0.0_dp, dp), result, p0, p1, &
            k_offset, padding, .false., real_callback=detrend)
    end subroutine stft_detrend_callable_real

    subroutine stft_detrend_callable_complex(self, x, result, detrend, p0, p1, k_offset, padding)
        class(short_time_fft), intent(in) :: self !! initialized ShortTimeFFT plan
        complex(dp), intent(in) :: x(:) !! time-domain complex samples
        type(short_time_fft_result), intent(out) :: result !! complex STFT and coordinates
        procedure(stft_complex_detrend_callback) :: detrend !! callable on each complex slice
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: p1 !! last time slice, exclusive
        integer, intent(in), optional :: k_offset !! sample offset
        character(*), intent(in), optional :: padding !! zeros, edge, even or odd

        call stft_kernel(self, x, result, p0, p1, k_offset, padding, &
            .true., complex_callback=detrend)
    end subroutine stft_detrend_callable_complex

    subroutine remove_trend(x, trend)
        complex(dp), intent(inout) :: x(:) !! time slice to detrend
        character(*), intent(in) :: trend !! constant, linear, none
        integer :: n, j
        complex(dp) :: mu, slope
        real(dp) :: center, denominator

        n = size(x)
        if (n < 1 .or. trend == 'none') return
        mu = sum(x)/real(n,dp)
        if (trend == 'constant') then
            x = x - mu
            return
        end if
        center = real(n+1,dp)*0.5_dp
        denominator = real(n,dp)*(real(n,dp)**2-1.0_dp)/12.0_dp
        slope = (0.0_dp,0.0_dp)
        if (n > 1) then
            do j = 1,n
                slope = slope + (real(j,dp)-center)*(x(j)-mu)
            end do
            slope = slope/denominator
        end if
        do j=1,n
            x(j) = x(j) - mu - slope*(real(j,dp)-center)
        end do
    end subroutine remove_trend

    pure function extended_sample(x, k, pad) result(value)
        complex(dp), intent(in) :: x(:) !! samples, n>=1
        integer, intent(in) :: k !! zero-based possibly out-of-bounds position
        character(*), intent(in) :: pad !! padding rule
        complex(dp) :: value, offset
        integer :: idx, n, signval

        n = size(x)
        if (k>=0 .and. k<n) then
            value = x(k+1)
            return
        end if
        select case (pad)
        case ('zeros')
            value = (0.0_dp,0.0_dp)
        case ('edge')
            value = x(max(1,min(n,k+1)))
        case ('even','odd')
            if (n == 1) then
                value = x(1)
                return
            end if
            idx = k
            signval = 1
            offset = (0.0_dp,0.0_dp)
            do while (idx < 0 .or. idx >= n)
                if (idx < 0) then
                    if (pad == 'odd') then
                        offset = offset + real(signval,dp)*2.0_dp*x(1)
                        signval = -signval
                    end if
                    idx = -idx
                else
                    if (pad == 'odd') then
                        offset = offset + real(signval,dp)*2.0_dp*x(n)
                        signval = -signval
                    end if
                    idx = 2*(n-1)-idx
                end if
            end do
            value = offset + real(signval,dp)*x(idx+1)
        case default
            value = (0.0_dp,0.0_dp)
        end select
    end function extended_sample

    subroutine stft_kernel(self, x, result, p0, p1, k_offset, padding, complex_input, detrend, &
            real_callback, complex_callback)
        class(short_time_fft), intent(in) :: self !! configured FFT plan
        complex(dp), intent(in) :: x(:) !! time-domain data
        type(short_time_fft_result), intent(out) :: result !! spectra and coordinates
        integer, intent(in), optional :: p0 !! zero-based first slice
        integer, intent(in), optional :: p1 !! exclusive last slice
        integer, intent(in), optional :: k_offset !! additional window offset
        character(*), intent(in), optional :: padding !! boundary mode
        logical, intent(in) :: complex_input !! rejects complex in one-sided mode
        character(*), intent(in), optional :: detrend !! none, linear, or constant
        procedure(stft_real_detrend_callback), optional :: real_callback !! per-slice real transform
        procedure(stft_complex_detrend_callback), optional :: complex_callback !! per-slice complex transform
        complex(dp), allocatable :: bins(:), shifted(:), samples(:)
        real(dp), allocatable :: real_samples(:)
        integer :: n, nw, nf, nfreq, first, last, p, j, q, shift, zero_offset, k, ksrc
        character(:), allocatable :: pad, trend
        real(dp) :: fac

        result%status = short_fft_invalid_input
        n = size(x)
        if (self%status /= short_fft_success) return
        nw = size(self%win)
        nf = self%mfft
        if (n < nw-nw/2) return
        if (complex_input .and. (self%fft_mode == 'onesided' .or. self%fft_mode == 'onesided2X')) return
        if (.not. all(ieee_is_finite(real(x,dp))) .or. .not. all(ieee_is_finite(aimag(x)))) return
        pad = 'zeros'
        if (present(padding)) pad = trim(padding)
        select case(pad)
        case ('zeros','edge','even','odd')
        case default
            return
        end select
        trend = 'none'
        if (present(detrend)) trend = trim(detrend)
        if (present(real_callback) .and. present(complex_callback)) return
        if (present(detrend) .and. (present(real_callback) .or. present(complex_callback))) return
        if (complex_input .and. present(real_callback)) return
        if (.not. complex_input .and. present(complex_callback)) return
        select case(trend)
        case ('none','constant','linear')
        case default
            return
        end select
        first = self%p_min
        last = self%p_max(n)
        if (present(p0)) first = p0
        if (present(p1)) last = p1
        if (first < self%p_min .or. last > self%p_max(n) .or. first >= last) return
        zero_offset = 0
        if (present(k_offset)) zero_offset = k_offset
        nfreq = nf
        if (self%fft_mode == 'onesided' .or. self%fft_mode == 'onesided2X') nfreq = nf/2+1
        allocate(result%values(nfreq,last-first),result%frequencies(nfreq),result%times(last-first))
        allocate(bins(nf),shifted(nf),samples(nw))
        if (present(real_callback)) allocate(real_samples(nw))
        shift = modulo(self%phase_shift+nw/2,nw)
        do p = first,last-1
            bins = (0.0_dp,0.0_dp)
            do j = 1,nw
                k = p*self%hop-nw/2+j-1+zero_offset
                samples(j) = extended_sample(x,k,pad)
            end do
            if (trend /= 'none') call remove_trend(samples,trend)
            if (present(real_callback)) then
                real_samples = real(samples,dp)
                call real_callback(real_samples)
                if (.not. all(ieee_is_finite(real_samples))) then
                    deallocate(result%values,result%frequencies,result%times)
                    return
                end if
                samples = cmplx(real_samples,0.0_dp,dp)
            end if
            if (present(complex_callback)) then
                call complex_callback(samples)
                if (.not. all(ieee_is_finite(real(samples,dp))) .or. &
                    .not. all(ieee_is_finite(aimag(samples)))) then
                    deallocate(result%values,result%frequencies,result%times)
                    return
                end if
            end if
            bins(1:nw) = samples*conjg(self%win)
            ! scipy.signal.ShortTimeFFT: pad first, then circularly shift the FFT input.
            do j = 1,nf
                ksrc = modulo(j-1+shift,nf)+1
                shifted(j) = bins(ksrc)
            end do
            call forward_fft(shifted)
            if (self%fft_mode == 'centered') then
                do q = 0,nf-1
                    result%values(q+1,p-first+1) = shifted(modulo(q+(nf+1)/2,nf)+1)
                end do
            else
                result%values(:,p-first+1) = shifted(1:nfreq)
            end if
            if (self%fft_mode == 'onesided2X') then
                fac = 2.0_dp
                if (self%scaling == 'psd') fac = sqrt(2.0_dp)
                do j = 2,nfreq
                    if (mod(nf,2) == 0 .and. j == nfreq) cycle
                    result%values(j,p-first+1) = fac*result%values(j,p-first+1)
                end do
            end if
            result%times(p-first+1) = real(p*self%hop+zero_offset,dp)/self%fs
        end do
        do q = 0,nfreq-1
            k = q
            if (self%fft_mode == 'twosided' .and. q >= (nf+1)/2) k = q-nf
            if (self%fft_mode == 'centered') k = q-nf/2
            result%frequencies(q+1) = real(k,dp)*self%fs/real(nf,dp)
        end do
        result%p0 = first
        result%status = short_fft_success
    end subroutine stft_kernel

    subroutine istft_real(self, s, x, status, k0, k1, p0)
        class(short_time_fft), intent(in) :: self !! transform settings
        complex(dp), intent(in) :: s(:,:) !! STFT coefficients
        real(dp), allocatable, intent(out) :: x(:) !! reconstructed real samples
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: k0 !! first output sample (zero-based)
        integer, intent(in), optional :: k1 !! exclusive upper output index
        integer, intent(in), optional :: p0 !! first input slice (default p_min)
        complex(dp), allocatable :: z(:)

        allocate(x(0))
        status = short_fft_invalid_input
        if (self%fft_mode /= 'onesided' .and. self%fft_mode /= 'onesided2X') return
        call inverse_kernel(self,s,z,status,k0,k1,p0)
        if (status == short_fft_success) x = real(z,dp)
    end subroutine istft_real

    subroutine istft_complex(self, s, x, status, k0, k1, p0)
        class(short_time_fft), intent(in) :: self !! transform settings
        complex(dp), intent(in) :: s(:,:) !! STFT matrix
        complex(dp), allocatable, intent(out) :: x(:) !! reconstructed complex samples
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: k0 !! inclusive first sample
        integer, intent(in), optional :: k1 !! exclusive last sample
        integer, intent(in), optional :: p0 !! first slice index

        status = short_fft_invalid_input
        allocate(x(0))
        if (self%fft_mode == 'onesided' .or. self%fft_mode == 'onesided2X') return
        call inverse_kernel(self,s,x,status,k0,k1,p0)
    end subroutine istft_complex

    subroutine inverse_kernel(self,s,x,status,k0,k1,p0)
        class(short_time_fft), intent(in) :: self !! inverse plan
        complex(dp), intent(in) :: s(:,:) !! bin by time-slice coefficients
        complex(dp), allocatable, intent(out) :: x(:) !! signal reconstruction
        integer, intent(out) :: status !! operation status
        integer, intent(in), optional :: k0 !! inclusive output start index
        integer, intent(in), optional :: k1 !! exclusive output end index
        integer, intent(in), optional :: p0 !! first coefficient slice
        complex(dp), allocatable :: bins(:), reconstructed(:)
        integer :: nw, nf, nfreq, start, finish, first, last, p, i, q, shift, jj
        real(dp) :: fac

        allocate(x(0))
        status = short_fft_invalid_input
        if (self%status /= short_fft_success) return
        nw = size(self%win)
        nf = self%mfft
        nfreq = nf
        if (self%fft_mode == 'onesided' .or. self%fft_mode == 'onesided2X') nfreq = nf/2+1
        if (size(s,1) /= nfreq .or. size(s,2) < 1) return
        if (.not. all(ieee_is_finite(real(s,dp))) .or. .not. all(ieee_is_finite(aimag(s)))) return
        first = self%p_min
        if (present(p0)) first = p0
        start = 0
        if (present(k0)) start = k0
        finish = (first+size(s,2)-1)*self%hop+nw-nw/2
        if (present(k1)) finish = k1
        if (finish <= start) return
        if (start < first*self%hop-nw/2 .or. &
            finish > (first+size(s,2)-1)*self%hop+nw-nw/2) return
        deallocate(x)
        allocate(x(finish-start),bins(nf),reconstructed(nf))
        x = (0.0_dp,0.0_dp)
        shift = modulo(self%phase_shift+nw/2,nw)
        do p = 1,size(s,2)
            bins = (0.0_dp,0.0_dp)
            select case (self%fft_mode)
            case ('centered')
                do q = 0,nf-1
                    bins(modulo(q+(nf+1)/2,nf)+1) = s(q+1,p)
                end do
            case ('twosided')
                bins = s(:,p)
            case ('onesided','onesided2X')
                bins(1:nfreq) = s(:,p)
                if (self%fft_mode == 'onesided2X') then
                    fac = 2.0_dp
                    if (self%scaling == 'psd') fac = sqrt(2.0_dp)
                    do q = 2,nfreq
                        if (mod(nf,2) == 0 .and. q == nfreq) cycle
                        bins(q) = bins(q)/fac
                    end do
                end if
                do q = nfreq+1,nf
                    bins(q) = conjg(bins(nf-q+2))
                end do
            end select
            reconstructed = conjg(bins)
            call forward_fft(reconstructed)
            reconstructed = conjg(reconstructed)/real(nf,dp)
            do i = 0,nw-1
                jj = modulo(i-shift,nf)+1
                q = (first+p-1)*self%hop-nw/2+i
                if (q >= start .and. q < finish) then
                    x(q-start+1) = x(q-start+1)+reconstructed(jj)*self%dual_win(i+1)
                end if
            end do
        end do
        status = short_fft_success
    end subroutine inverse_kernel

    subroutine spectrogram_real(self, x, result, p0, p1, padding)
        class(short_time_fft), intent(in) :: self !! FFT settings
        real(dp), intent(in) :: x(:) !! real time signal
        type(short_time_spectrogram_result), intent(out) :: result !! power and axes
        integer, intent(in), optional :: p0 !! first time slice
        integer, intent(in), optional :: p1 !! last slice, exclusive
        character(*), intent(in), optional :: padding !! boundary mode
        type(short_time_fft_result) :: spectrum
        call self%stft(x,spectrum,p0,p1,padding=padding)
        call from_stft(spectrum,result)
    end subroutine spectrogram_real

    subroutine spectrogram_complex(self, x, result, p0, p1, padding)
        class(short_time_fft), intent(in) :: self !! FFT settings
        complex(dp), intent(in) :: x(:) !! complex signal
        type(short_time_spectrogram_result), intent(out) :: result !! squared magnitudes
        integer, intent(in), optional :: p0 !! first slice index
        integer, intent(in), optional :: p1 !! last slice, exclusive
        character(*), intent(in), optional :: padding !! padding mode
        type(short_time_fft_result) :: spectrum
        call self%stft(x,spectrum,p0,p1,padding=padding)
        call from_stft(spectrum,result)
    end subroutine spectrogram_complex

    subroutine from_stft(spectrum,result)
        type(short_time_fft_result), intent(in) :: spectrum !! complex coefficients
        type(short_time_spectrogram_result), intent(out) :: result !! power spectrum
        result%status = spectrum%status
        if (spectrum%status /= short_fft_success) return
        result%power = abs(spectrum%values)**2
        result%frequencies = spectrum%frequencies
        result%times = spectrum%times
        result%p0 = spectrum%p0
    end subroutine from_stft

    subroutine cross_spectrogram_real(self, x, y, values, status, p0, p1, padding)
        class(short_time_fft), intent(in) :: self !! STFT plan
        real(dp), intent(in) :: x(:) !! first real signal
        real(dp), intent(in) :: y(:) !! second real signal
        complex(dp), allocatable, intent(out) :: values(:,:) !! Sx * conjg(Sy)
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: p0 !! first time slice
        integer, intent(in), optional :: p1 !! end time slice
        character(*), intent(in), optional :: padding !! extension mode
        type(short_time_fft_result) :: sx, sy
        allocate(values(0,0))
        status = short_fft_invalid_input
        if (size(x) /= size(y)) return
        call self%stft(x,sx,p0,p1,padding=padding)
        if (sx%status /= short_fft_success) return
        call self%stft(y,sy,p0,p1,padding=padding)
        if (sy%status /= short_fft_success) return
        values = sx%values*conjg(sy%values)
        status = short_fft_success
    end subroutine cross_spectrogram_real

    subroutine cross_spectrogram_complex(self, x, y, values, status, p0, p1, padding)
        class(short_time_fft), intent(in) :: self !! STFT plan
        complex(dp), intent(in) :: x(:) !! first signal
        complex(dp), intent(in) :: y(:) !! second signal
        complex(dp), allocatable, intent(out) :: values(:,:) !! Sx * conjg(Sy)
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: p0 !! first time slice
        integer, intent(in), optional :: p1 !! final slice, exclusive
        character(*), intent(in), optional :: padding !! extension mode
        type(short_time_fft_result) :: sx, sy
        allocate(values(0,0))
        status = short_fft_invalid_input
        if (size(x) /= size(y)) return
        call self%stft(x,sx,p0,p1,padding=padding)
        if (sx%status /= short_fft_success) return
        call self%stft(y,sy,p0,p1,padding=padding)
        if (sy%status /= short_fft_success) return
        values = sx%values*conjg(sy%values)
        status = short_fft_success
    end subroutine cross_spectrogram_complex

    subroutine cross_spectrogram_real_complex(self, x, y, values, status, p0, p1, padding)
        class(short_time_fft), intent(in) :: self !! initialized analysis plan
        real(dp), intent(in) :: x(:) !! first input signal
        complex(dp), intent(in) :: y(:) !! second input signal
        complex(dp), allocatable, intent(out) :: values(:,:) !! Sx * conjg(Sy)
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: p0 !! first slice index, inclusive
        integer, intent(in), optional :: p1 !! final slice index, exclusive
        character(*), intent(in), optional :: padding !! extension mode
        type(short_time_fft_result) :: sx, sy

        allocate(values(0,0))
        status = short_fft_invalid_input
        if (size(x) /= size(y)) return
        ! Complex input requires the full FFT, regardless of operand order.
        if (self%fft_mode == 'onesided' .or. self%fft_mode == 'onesided2X') return
        call self%stft(x,sx,p0,p1,padding=padding)
        if (sx%status /= short_fft_success) return
        call self%stft(y,sy,p0,p1,padding=padding)
        if (sy%status /= short_fft_success) return
        values = sx%values * conjg(sy%values)
        status = short_fft_success
    end subroutine cross_spectrogram_real_complex

    subroutine cross_spectrogram_complex_real(self, x, y, values, status, p0, p1, padding)
        class(short_time_fft), intent(in) :: self !! initialized analysis plan
        complex(dp), intent(in) :: x(:) !! first input signal
        real(dp), intent(in) :: y(:) !! second input signal
        complex(dp), allocatable, intent(out) :: values(:,:) !! Sx * conjg(Sy)
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: p0 !! first slice index, inclusive
        integer, intent(in), optional :: p1 !! final slice index, exclusive
        character(*), intent(in), optional :: padding !! extension mode
        type(short_time_fft_result) :: sx, sy

        allocate(values(0,0))
        status = short_fft_invalid_input
        if (size(x) /= size(y)) return
        ! Complex input requires the full FFT, regardless of operand order.
        if (self%fft_mode == 'onesided' .or. self%fft_mode == 'onesided2X') return
        call self%stft(x,sx,p0,p1,padding=padding)
        if (sx%status /= short_fft_success) return
        call self%stft(y,sy,p0,p1,padding=padding)
        if (sy%status /= short_fft_success) return
        values = sx%values * conjg(sy%values)
        status = short_fft_success
    end subroutine cross_spectrogram_complex_real

    ! Axis values follow Fortran convention (positive: 1..rank). Negative
    ! axes follow NumPy convention (-1 is the last dimension). Axis 0 is
    ! deliberately invalid to avoid ambiguity between the two conventions.
    pure integer function fft_resolve_axis(input_axis, rank, default_axis) result(axis)
        integer, intent(in), optional :: input_axis !! requested axis
        integer, intent(in) :: rank !! rank of input array
        integer, intent(in) :: default_axis !! default axis, 1-based
        axis = default_axis
        if (present(input_axis)) axis = input_axis
        if (axis < 0) axis = axis + rank + 1
        if (axis < 1 .or. axis > rank) axis = 0
    end function fft_resolve_axis

    subroutine stft_nd2_real(self, x, result, p0, p1, k_offset, padding, axis)
        class(short_time_fft), intent(in) :: self !! transform configuration
        real(dp), intent(in) :: x(:,:) !! input signal batch
        type(short_time_fft_result_3d), intent(out) :: result !! frequency and time coefficients
        integer, intent(in), optional :: axis !! signal axis, 1-based or negative
        integer, intent(in), optional :: p0 !! first slice, zero-based
        integer, intent(in), optional :: p1 !! exclusive upper time slice
        integer, intent(in), optional :: k_offset !! input sample offset
        character(*), intent(in), optional :: padding !! extension mode
        call stft_nd2_kernel(self, cmplx(x,0.0_dp,dp), result, axis, p0, p1, &
            k_offset, padding, .false.)
    end subroutine stft_nd2_real

    subroutine stft_detrend_nd2_real(self, x, result, detrend, p0, p1, k_offset, padding, axis)
        class(short_time_fft), intent(in) :: self !! transform configuration
        real(dp), intent(in) :: x(:,:) !! input signal batch
        type(short_time_fft_result_3d), intent(out) :: result !! frequency and time coefficients
        character(*), intent(in) :: detrend !! detrending mode
        integer, intent(in), optional :: axis !! signal axis, 1-based or negative
        integer, intent(in), optional :: p0 !! first slice, zero-based
        integer, intent(in), optional :: p1 !! exclusive upper time slice
        integer, intent(in), optional :: k_offset !! input sample offset
        character(*), intent(in), optional :: padding !! extension mode
        call stft_nd2_kernel(self, cmplx(x,0.0_dp,dp), result, axis, p0, p1, &
            k_offset, padding, .false., detrend)
    end subroutine stft_detrend_nd2_real

    subroutine stft_nd2_complex(self, x, result, p0, p1, k_offset, padding, axis)
        class(short_time_fft), intent(in) :: self !! transform configuration
        complex(dp), intent(in) :: x(:,:) !! input signal batch
        type(short_time_fft_result_3d), intent(out) :: result !! frequency and time coefficients
        integer, intent(in), optional :: axis !! signal axis, 1-based or negative
        integer, intent(in), optional :: p0 !! first slice, zero-based
        integer, intent(in), optional :: p1 !! exclusive upper time slice
        integer, intent(in), optional :: k_offset !! input sample offset
        character(*), intent(in), optional :: padding !! extension mode
        call stft_nd2_kernel(self, x, result, axis, p0, p1, &
            k_offset, padding, .true.)
    end subroutine stft_nd2_complex

    subroutine stft_detrend_nd2_complex(self, x, result, detrend, p0, p1, k_offset, padding, axis)
        class(short_time_fft), intent(in) :: self !! transform configuration
        complex(dp), intent(in) :: x(:,:) !! input signal batch
        type(short_time_fft_result_3d), intent(out) :: result !! frequency and time coefficients
        character(*), intent(in) :: detrend !! detrending mode
        integer, intent(in), optional :: axis !! signal axis, 1-based or negative
        integer, intent(in), optional :: p0 !! first slice, zero-based
        integer, intent(in), optional :: p1 !! exclusive upper time slice
        integer, intent(in), optional :: k_offset !! input sample offset
        character(*), intent(in), optional :: padding !! extension mode
        call stft_nd2_kernel(self, x, result, axis, p0, p1, &
            k_offset, padding, .true., detrend)
    end subroutine stft_detrend_nd2_complex

    subroutine stft_nd3_real(self, x, result, p0, p1, k_offset, padding, axis)
        class(short_time_fft), intent(in) :: self !! transform configuration
        real(dp), intent(in) :: x(:,:,:) !! input signal batch
        type(short_time_fft_result_4d), intent(out) :: result !! frequency and time coefficients
        integer, intent(in), optional :: axis !! signal axis, 1-based or negative
        integer, intent(in), optional :: p0 !! first slice, zero-based
        integer, intent(in), optional :: p1 !! exclusive upper time slice
        integer, intent(in), optional :: k_offset !! input sample offset
        character(*), intent(in), optional :: padding !! extension mode
        call stft_nd3_kernel(self, cmplx(x,0.0_dp,dp), result, axis, p0, p1, &
            k_offset, padding, .false.)
    end subroutine stft_nd3_real

    subroutine stft_detrend_nd3_real(self, x, result, detrend, p0, p1, k_offset, padding, axis)
        class(short_time_fft), intent(in) :: self !! transform configuration
        real(dp), intent(in) :: x(:,:,:) !! input signal batch
        type(short_time_fft_result_4d), intent(out) :: result !! frequency and time coefficients
        character(*), intent(in) :: detrend !! detrending mode
        integer, intent(in), optional :: axis !! signal axis, 1-based or negative
        integer, intent(in), optional :: p0 !! first slice, zero-based
        integer, intent(in), optional :: p1 !! exclusive upper time slice
        integer, intent(in), optional :: k_offset !! input sample offset
        character(*), intent(in), optional :: padding !! extension mode
        call stft_nd3_kernel(self, cmplx(x,0.0_dp,dp), result, axis, p0, p1, &
            k_offset, padding, .false., detrend)
    end subroutine stft_detrend_nd3_real

    subroutine stft_nd3_complex(self, x, result, p0, p1, k_offset, padding, axis)
        class(short_time_fft), intent(in) :: self !! transform configuration
        complex(dp), intent(in) :: x(:,:,:) !! input signal batch
        type(short_time_fft_result_4d), intent(out) :: result !! frequency and time coefficients
        integer, intent(in), optional :: axis !! signal axis, 1-based or negative
        integer, intent(in), optional :: p0 !! first slice, zero-based
        integer, intent(in), optional :: p1 !! exclusive upper time slice
        integer, intent(in), optional :: k_offset !! input sample offset
        character(*), intent(in), optional :: padding !! extension mode
        call stft_nd3_kernel(self, x, result, axis, p0, p1, &
            k_offset, padding, .true.)
    end subroutine stft_nd3_complex

    subroutine stft_detrend_nd3_complex(self, x, result, detrend, p0, p1, k_offset, padding, axis)
        class(short_time_fft), intent(in) :: self !! transform configuration
        complex(dp), intent(in) :: x(:,:,:) !! input signal batch
        type(short_time_fft_result_4d), intent(out) :: result !! frequency and time coefficients
        character(*), intent(in) :: detrend !! detrending mode
        integer, intent(in), optional :: axis !! signal axis, 1-based or negative
        integer, intent(in), optional :: p0 !! first slice, zero-based
        integer, intent(in), optional :: p1 !! exclusive upper time slice
        integer, intent(in), optional :: k_offset !! input sample offset
        character(*), intent(in), optional :: padding !! extension mode
        call stft_nd3_kernel(self, x, result, axis, p0, p1, &
            k_offset, padding, .true., detrend)
    end subroutine stft_detrend_nd3_complex

    subroutine stft_detrend_callable_nd2_real(self, x, result, detrend, p0, p1, k_offset, padding, axis)
        class(short_time_fft), intent(in) :: self !! initialized analysis plan
        real(dp), intent(in) :: x(:,:) !! rank-2 input samples
        type(short_time_fft_result_3d), intent(out) :: result !! batch STFT coefficients
        procedure(stft_real_detrend_callback) :: detrend !! per-channel, per-slice transform
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: p1 !! exclusive last time slice
        integer, intent(in), optional :: k_offset !! sample offset
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! one-based or negative signal dimension

        call stft_nd2_kernel(self, cmplx(x,0.0_dp,dp), result, axis, p0, p1, &
            k_offset, padding, .false., real_callback=detrend)
    end subroutine stft_detrend_callable_nd2_real

    subroutine stft_detrend_callable_nd2_complex(self, x, result, detrend, p0, p1, k_offset, padding, axis)
        class(short_time_fft), intent(in) :: self !! initialized analysis plan
        complex(dp), intent(in) :: x(:,:) !! rank-2 input samples
        type(short_time_fft_result_3d), intent(out) :: result !! batch STFT coefficients
        procedure(stft_complex_detrend_callback) :: detrend !! per-channel, per-slice transform
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: p1 !! exclusive last time slice
        integer, intent(in), optional :: k_offset !! sample offset
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! one-based or negative signal dimension

        call stft_nd2_kernel(self, x, result, axis, p0, p1, &
            k_offset, padding, .true., complex_callback=detrend)
    end subroutine stft_detrend_callable_nd2_complex

    subroutine stft_detrend_callable_nd3_real(self, x, result, detrend, p0, p1, k_offset, padding, axis)
        class(short_time_fft), intent(in) :: self !! initialized analysis plan
        real(dp), intent(in) :: x(:,:,:) !! rank-3 input samples
        type(short_time_fft_result_4d), intent(out) :: result !! batch STFT coefficients
        procedure(stft_real_detrend_callback) :: detrend !! per-channel, per-slice transform
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: p1 !! exclusive last time slice
        integer, intent(in), optional :: k_offset !! sample offset
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! one-based or negative signal dimension

        call stft_nd3_kernel(self, cmplx(x,0.0_dp,dp), result, axis, p0, p1, &
            k_offset, padding, .false., real_callback=detrend)
    end subroutine stft_detrend_callable_nd3_real

    subroutine stft_detrend_callable_nd3_complex(self, x, result, detrend, p0, p1, k_offset, padding, axis)
        class(short_time_fft), intent(in) :: self !! initialized analysis plan
        complex(dp), intent(in) :: x(:,:,:) !! rank-3 input samples
        type(short_time_fft_result_4d), intent(out) :: result !! batch STFT coefficients
        procedure(stft_complex_detrend_callback) :: detrend !! per-channel, per-slice transform
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: p1 !! exclusive last time slice
        integer, intent(in), optional :: k_offset !! sample offset
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! one-based or negative signal dimension

        call stft_nd3_kernel(self, x, result, axis, p0, p1, &
            k_offset, padding, .true., complex_callback=detrend)
    end subroutine stft_detrend_callable_nd3_complex

    subroutine stft_nd2_kernel(self, x, result, axis, p0, p1, k_offset, padding, complex_input, detrend, &
            real_callback, complex_callback)
        class(short_time_fft), intent(in) :: self !! configured plan
        complex(dp), intent(in) :: x(:,:) !! batch of input signals
        type(short_time_fft_result_3d), intent(out) :: result !! 3D spectrum
        integer, intent(in), optional :: axis !! input signal dimension
        integer, intent(in), optional :: p0 !! first zero-based slice
        integer, intent(in), optional :: p1 !! exclusive final slice
        integer, intent(in), optional :: k_offset !! offset in samples
        character(*), intent(in), optional :: padding !! zeros, edge, even, odd
        logical, intent(in) :: complex_input !! whether original signal is complex
        character(*), intent(in), optional :: detrend !! detrending specification
        procedure(stft_real_detrend_callback), optional :: real_callback !! real slice callback
        procedure(stft_complex_detrend_callback), optional :: complex_callback !! complex slice callback
        type(short_time_fft_result) :: frame
        integer :: a, nb, b

        result%status = short_fft_invalid_input
        a = fft_resolve_axis(axis,2,2)
        if (a == 0 .or. any(shape(x) < 1)) return
        nb = size(x,3-a)
        do b = 1,nb
            if (a == 1) then
                call stft_kernel(self,x(:,b),frame,p0,p1,k_offset,padding,complex_input,detrend, &
                    real_callback,complex_callback)
            else
                call stft_kernel(self,x(b,:),frame,p0,p1,k_offset,padding,complex_input,detrend, &
                    real_callback,complex_callback)
            end if
            if (frame%status /= short_fft_success) then
                if (allocated(result%values)) deallocate(result%values)
                if (allocated(result%times)) deallocate(result%times)
                if (allocated(result%frequencies)) deallocate(result%frequencies)
                return
            end if
            if (b == 1) then
                if (a == 1) then
                    allocate(result%values(size(frame%values,1),nb,size(frame%values,2)))
                else
                    allocate(result%values(nb,size(frame%values,1),size(frame%values,2)))
                end if
                result%times = frame%times
                result%frequencies = frame%frequencies
                result%p0 = frame%p0
            end if
            if (a == 1) then
                result%values(:,b,:) = frame%values
            else
                result%values(b,:,:) = frame%values
            end if
        end do
        result%status = short_fft_success
    end subroutine stft_nd2_kernel

    subroutine stft_nd3_kernel(self, x, result, axis, p0, p1, k_offset, padding, complex_input, detrend, &
            real_callback, complex_callback)
        class(short_time_fft), intent(in) :: self !! configured plan
        complex(dp), intent(in) :: x(:,:,:) !! rank-3 input signals
        type(short_time_fft_result_4d), intent(out) :: result !! rank-4 spectrum
        integer, intent(in), optional :: axis !! signal axis, defaults to 3
        integer, intent(in), optional :: p0 !! first time slice
        integer, intent(in), optional :: p1 !! exclusive final slice
        integer, intent(in), optional :: k_offset !! sample offset
        character(*), intent(in), optional :: padding !! boundary padding
        logical, intent(in) :: complex_input !! original data type
        character(*), intent(in), optional :: detrend !! detrending mode
        procedure(stft_real_detrend_callback), optional :: real_callback !! real slice callback
        procedure(stft_complex_detrend_callback), optional :: complex_callback !! complex slice callback
        type(short_time_fft_result) :: frame
        integer :: a, nb1, nb2, b1, b2

        result%status = short_fft_invalid_input
        a = fft_resolve_axis(axis,3,3)
        if (a == 0 .or. any(shape(x) < 1)) return
        select case (a)
        case (1)
            nb1 = size(x,2)
            nb2 = size(x,3)
        case (2)
            nb1 = size(x,1)
            nb2 = size(x,3)
        case (3)
            nb1 = size(x,1)
            nb2 = size(x,2)
        end select
        do b2 = 1,nb2
            do b1 = 1,nb1
                select case (a)
                case (1)
                    call stft_kernel(self,x(:,b1,b2),frame,p0,p1,k_offset,padding,complex_input,detrend, &
                    real_callback,complex_callback)
                case (2)
                    call stft_kernel(self,x(b1,:,b2),frame,p0,p1,k_offset,padding,complex_input,detrend, &
                    real_callback,complex_callback)
                case (3)
                    call stft_kernel(self,x(b1,b2,:),frame,p0,p1,k_offset,padding,complex_input,detrend, &
                    real_callback,complex_callback)
                end select
                if (frame%status /= short_fft_success) then
                if (allocated(result%values)) deallocate(result%values)
                if (allocated(result%times)) deallocate(result%times)
                if (allocated(result%frequencies)) deallocate(result%frequencies)
                return
            end if
                if (b1 == 1 .and. b2 == 1) then
                    select case (a)
                    case (1)
                        allocate(result%values(size(frame%values,1),nb1,nb2,size(frame%values,2)))
                    case (2)
                        allocate(result%values(nb1,size(frame%values,1),nb2,size(frame%values,2)))
                    case (3)
                        allocate(result%values(nb1,nb2,size(frame%values,1),size(frame%values,2)))
                    end select
                    result%times = frame%times
                    result%frequencies = frame%frequencies
                    result%p0 = frame%p0
                end if
                select case (a)
                case (1)
                    result%values(:,b1,b2,:) = frame%values
                case (2)
                    result%values(b1,:,b2,:) = frame%values
                case (3)
                    result%values(b1,b2,:,:) = frame%values
                end select
            end do
        end do
        result%status = short_fft_success
    end subroutine stft_nd3_kernel

    subroutine istft_nd3_real(self, s, x, status, k0, k1, p0, f_axis, t_axis)
        class(short_time_fft), intent(in) :: self !! transform configuration
        complex(dp), intent(in) :: s(:,:,:) !! STFT coefficients
        real(dp), allocatable, intent(out) :: x(:,:) !! reconstructed signal
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: k0 !! inclusive first sample
        integer, intent(in), optional :: k1 !! exclusive upper sample
        integer, intent(in), optional :: p0 !! first input time slice
        integer, intent(in), optional :: f_axis !! 1-based or negative frequency axis
        integer, intent(in), optional :: t_axis !! 1-based or negative time axis
        complex(dp), allocatable :: z(:,:)

        status = short_fft_invalid_input
        allocate(x(0,0))
        if (self%fft_mode /= 'onesided' .and. self%fft_mode /= 'onesided2X') return
        call inverse_nd3_kernel(self,s,z,status,k0,k1,p0,f_axis,t_axis)
        if (status == short_fft_success) x = real(z,dp)
    end subroutine istft_nd3_real

    subroutine istft_nd3_complex(self, s, x, status, k0, k1, p0, f_axis, t_axis)
        class(short_time_fft), intent(in) :: self !! transform configuration
        complex(dp), intent(in) :: s(:,:,:) !! STFT coefficients
        complex(dp), allocatable, intent(out) :: x(:,:) !! reconstructed signal
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: k0 !! inclusive first sample
        integer, intent(in), optional :: k1 !! exclusive upper sample
        integer, intent(in), optional :: p0 !! first input time slice
        integer, intent(in), optional :: f_axis !! 1-based or negative frequency axis
        integer, intent(in), optional :: t_axis !! 1-based or negative time axis
        complex(dp), allocatable :: z(:,:)

        status = short_fft_invalid_input
        allocate(x(0,0))
        if (self%fft_mode == 'onesided' .or. self%fft_mode == 'onesided2X') return
        call inverse_nd3_kernel(self,s,z,status,k0,k1,p0,f_axis,t_axis)
        if (status == short_fft_success) x = z
    end subroutine istft_nd3_complex

    subroutine istft_nd4_real(self, s, x, status, k0, k1, p0, f_axis, t_axis)
        class(short_time_fft), intent(in) :: self !! transform configuration
        complex(dp), intent(in) :: s(:,:,:,:) !! STFT coefficients
        real(dp), allocatable, intent(out) :: x(:,:,:) !! reconstructed signal
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: k0 !! inclusive first sample
        integer, intent(in), optional :: k1 !! exclusive upper sample
        integer, intent(in), optional :: p0 !! first input time slice
        integer, intent(in), optional :: f_axis !! 1-based or negative frequency axis
        integer, intent(in), optional :: t_axis !! 1-based or negative time axis
        complex(dp), allocatable :: z(:,:,:)

        status = short_fft_invalid_input
        allocate(x(0,0,0))
        if (self%fft_mode /= 'onesided' .and. self%fft_mode /= 'onesided2X') return
        call inverse_nd4_kernel(self,s,z,status,k0,k1,p0,f_axis,t_axis)
        if (status == short_fft_success) x = real(z,dp)
    end subroutine istft_nd4_real

    subroutine istft_nd4_complex(self, s, x, status, k0, k1, p0, f_axis, t_axis)
        class(short_time_fft), intent(in) :: self !! transform configuration
        complex(dp), intent(in) :: s(:,:,:,:) !! STFT coefficients
        complex(dp), allocatable, intent(out) :: x(:,:,:) !! reconstructed signal
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: k0 !! inclusive first sample
        integer, intent(in), optional :: k1 !! exclusive upper sample
        integer, intent(in), optional :: p0 !! first input time slice
        integer, intent(in), optional :: f_axis !! 1-based or negative frequency axis
        integer, intent(in), optional :: t_axis !! 1-based or negative time axis
        complex(dp), allocatable :: z(:,:,:)

        status = short_fft_invalid_input
        allocate(x(0,0,0))
        if (self%fft_mode == 'onesided' .or. self%fft_mode == 'onesided2X') return
        call inverse_nd4_kernel(self,s,z,status,k0,k1,p0,f_axis,t_axis)
        if (status == short_fft_success) x = z
    end subroutine istft_nd4_complex

    subroutine inverse_nd3_kernel(self,s,x,status,k0,k1,p0,f_axis,t_axis)
        class(short_time_fft), intent(in) :: self !! configured plan
        complex(dp), intent(in) :: s(:,:,:) !! multidimensional STFT coefficients
        complex(dp), allocatable, intent(out) :: x(:,:) !! reconstructed samples
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: k0 !! inclusive first sample
        integer, intent(in), optional :: k1 !! exclusive end sample
        integer, intent(in), optional :: p0 !! first time-slice index
        integer, intent(in), optional :: f_axis !! frequency dimension in S
        integer, intent(in), optional :: t_axis !! time dimension in S
        complex(dp), allocatable :: frame(:,:), vec(:)
        integer :: f, t, b_axis, dims(3), out_axes(2), outdims(2)
        integer :: idx(3), outidx(2), b, q, p, k

        allocate(x(0,0))
        status = short_fft_invalid_input
        f = fft_resolve_axis(f_axis,3,2)
        t = fft_resolve_axis(t_axis,3,3)
        if (f == 0 .or. t == 0 .or. f == t) return
        dims = shape(s)
        if (any(dims < 1)) return
        b_axis = 6-f-t
        out_axes = pack([1,2,3],[1,2,3] /= t)
        allocate(frame(dims(f),dims(t)))
        do b = 1,dims(b_axis)
            idx(b_axis) = b
            do p=1,dims(t)
                idx(t) = p
                do q=1,dims(f)
                    idx(f) = q
                    frame(q,p) = s(idx(1),idx(2),idx(3))
                end do
            end do
            call inverse_kernel(self,frame,vec,status,k0,k1,p0)
            if (status /= short_fft_success) return
            if (b == 1) then
                outdims = pack(dims,[1,2,3] /= t)
                outdims(findloc(out_axes,f,dim=1)) = size(vec)
                deallocate(x)
                allocate(x(outdims(1),outdims(2)))
            end if
            outidx(findloc(out_axes,b_axis,dim=1)) = b
            do k = 1,size(vec)
                outidx(findloc(out_axes,f,dim=1)) = k
                x(outidx(1),outidx(2)) = vec(k)
            end do
        end do
        status = short_fft_success
    end subroutine inverse_nd3_kernel

    subroutine inverse_nd4_kernel(self,s,x,status,k0,k1,p0,f_axis,t_axis)
        class(short_time_fft), intent(in) :: self !! configured plan
        complex(dp), intent(in) :: s(:,:,:,:) !! rank-4 STFT coefficients
        complex(dp), allocatable, intent(out) :: x(:,:,:) !! rank-3 reconstructed samples
        integer, intent(out) :: status !! operation status
        integer, intent(in), optional :: k0 !! first output sample
        integer, intent(in), optional :: k1 !! exclusive final sample
        integer, intent(in), optional :: p0 !! first input time slice
        integer, intent(in), optional :: f_axis !! frequency axis
        integer, intent(in), optional :: t_axis !! time axis
        complex(dp), allocatable :: frame(:,:), vec(:)
        integer :: f, t, dims(4), out_axes(3), outdims(3), batch_axes(2)
        integer :: idx(4), outidx(3), b1, b2, q, p, k

        allocate(x(0,0,0))
        status = short_fft_invalid_input
        f = fft_resolve_axis(f_axis,4,3)
        t = fft_resolve_axis(t_axis,4,4)
        if (f == 0 .or. t == 0 .or. f == t) return
        dims = shape(s)
        if (any(dims < 1)) return
        batch_axes = pack([1,2,3,4],[1,2,3,4] /= f .and. [1,2,3,4] /= t)
        out_axes = pack([1,2,3,4],[1,2,3,4] /= t)
        allocate(frame(dims(f),dims(t)))
        do b2=1,dims(batch_axes(2))
            do b1=1,dims(batch_axes(1))
                idx(batch_axes) = [b1,b2]
                do p=1,dims(t)
                    idx(t) = p
                    do q=1,dims(f)
                        idx(f) = q
                        frame(q,p) = s(idx(1),idx(2),idx(3),idx(4))
                    end do
                end do
                call inverse_kernel(self,frame,vec,status,k0,k1,p0)
                if (status /= short_fft_success) return
                if (b1 == 1 .and. b2 == 1) then
                    outdims = pack(dims,[1,2,3,4] /= t)
                    outdims(findloc(out_axes,f,dim=1)) = size(vec)
                    deallocate(x)
                    allocate(x(outdims(1),outdims(2),outdims(3)))
                end if
                outidx(findloc(out_axes,batch_axes(1),dim=1)) = b1
                outidx(findloc(out_axes,batch_axes(2),dim=1)) = b2
                do k=1,size(vec)
                    outidx(findloc(out_axes,f,dim=1)) = k
                    x(outidx(1),outidx(2),outidx(3)) = vec(k)
                end do
            end do
        end do
        status = short_fft_success
    end subroutine inverse_nd4_kernel

    subroutine spectrogram_nd2_real(self,x,result,p0,p1,padding,axis)
        class(short_time_fft), intent(in) :: self !! configured plan
        real(dp), intent(in) :: x(:,:) !! time series
        type(short_time_spectrogram_result_3d), intent(out) :: result !! squared magnitude
        integer, intent(in), optional :: p0 !! starting time slice
        integer, intent(in), optional :: p1 !! exclusive final slice
        character(*), intent(in), optional :: padding !! extension mode
        integer, intent(in), optional :: axis !! input transform axis
        type(short_time_fft_result_3d) :: st

        call self%stft(x,st,axis=axis,p0=p0,p1=p1,padding=padding)
        result%status = st%status
        if (st%status /= short_fft_success) return
        result%power = abs(st%values)**2
        result%times = st%times
        result%frequencies = st%frequencies
        result%p0 = st%p0
    end subroutine spectrogram_nd2_real

    subroutine spectrogram_nd2_complex(self,x,result,p0,p1,padding,axis)
        class(short_time_fft), intent(in) :: self !! configured plan
        complex(dp), intent(in) :: x(:,:) !! time series
        type(short_time_spectrogram_result_3d), intent(out) :: result !! squared magnitude
        integer, intent(in), optional :: p0 !! starting time slice
        integer, intent(in), optional :: p1 !! exclusive final slice
        character(*), intent(in), optional :: padding !! extension mode
        integer, intent(in), optional :: axis !! input transform axis
        type(short_time_fft_result_3d) :: st

        call self%stft(x,st,axis=axis,p0=p0,p1=p1,padding=padding)
        result%status = st%status
        if (st%status /= short_fft_success) return
        result%power = abs(st%values)**2
        result%times = st%times
        result%frequencies = st%frequencies
        result%p0 = st%p0
    end subroutine spectrogram_nd2_complex

    subroutine spectrogram_nd3_real(self,x,result,p0,p1,padding,axis)
        class(short_time_fft), intent(in) :: self !! configured plan
        real(dp), intent(in) :: x(:,:,:) !! time series
        type(short_time_spectrogram_result_4d), intent(out) :: result !! squared magnitude
        integer, intent(in), optional :: p0 !! starting time slice
        integer, intent(in), optional :: p1 !! exclusive final slice
        character(*), intent(in), optional :: padding !! extension mode
        integer, intent(in), optional :: axis !! input transform axis
        type(short_time_fft_result_4d) :: st

        call self%stft(x,st,axis=axis,p0=p0,p1=p1,padding=padding)
        result%status = st%status
        if (st%status /= short_fft_success) return
        result%power = abs(st%values)**2
        result%times = st%times
        result%frequencies = st%frequencies
        result%p0 = st%p0
    end subroutine spectrogram_nd3_real

    subroutine spectrogram_nd3_complex(self,x,result,p0,p1,padding,axis)
        class(short_time_fft), intent(in) :: self !! configured plan
        complex(dp), intent(in) :: x(:,:,:) !! time series
        type(short_time_spectrogram_result_4d), intent(out) :: result !! squared magnitude
        integer, intent(in), optional :: p0 !! starting time slice
        integer, intent(in), optional :: p1 !! exclusive final slice
        character(*), intent(in), optional :: padding !! extension mode
        integer, intent(in), optional :: axis !! input transform axis
        type(short_time_fft_result_4d) :: st

        call self%stft(x,st,axis=axis,p0=p0,p1=p1,padding=padding)
        result%status = st%status
        if (st%status /= short_fft_success) return
        result%power = abs(st%values)**2
        result%times = st%times
        result%frequencies = st%frequencies
        result%p0 = st%p0
    end subroutine spectrogram_nd3_complex

    subroutine cross_spectrogram_nd2_real(self,x,y,values,status,p0,p1,padding,axis)
        class(short_time_fft), intent(in) :: self !! configured plan
        real(dp), intent(in) :: x(:,:) !! first time series batch
        real(dp), intent(in) :: y(:,:) !! second time series batch
        complex(dp), allocatable, intent(out) :: values(:,:,:) !! cross power spectrum
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: p0 !! first time slice
        integer, intent(in), optional :: p1 !! exclusive final time slice
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! time-series axis
        type(short_time_fft_result_3d) :: sx, sy

        allocate(values(0,0,0))
        status = short_fft_invalid_input
        if (any(shape(x) /= shape(y))) return
        call self%stft(x,sx,axis=axis,p0=p0,p1=p1,padding=padding)
        if (sx%status /= short_fft_success) return
        call self%stft(y,sy,axis=axis,p0=p0,p1=p1,padding=padding)
        if (sy%status /= short_fft_success) return
        values = sx%values*conjg(sy%values)
        status = short_fft_success
    end subroutine cross_spectrogram_nd2_real

    subroutine cross_spectrogram_nd2_complex(self,x,y,values,status,p0,p1,padding,axis)
        class(short_time_fft), intent(in) :: self !! configured plan
        complex(dp), intent(in) :: x(:,:) !! first time series batch
        complex(dp), intent(in) :: y(:,:) !! second time series batch
        complex(dp), allocatable, intent(out) :: values(:,:,:) !! cross power spectrum
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: p0 !! first time slice
        integer, intent(in), optional :: p1 !! exclusive final time slice
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! time-series axis
        type(short_time_fft_result_3d) :: sx, sy

        allocate(values(0,0,0))
        status = short_fft_invalid_input
        if (any(shape(x) /= shape(y))) return
        call self%stft(x,sx,axis=axis,p0=p0,p1=p1,padding=padding)
        if (sx%status /= short_fft_success) return
        call self%stft(y,sy,axis=axis,p0=p0,p1=p1,padding=padding)
        if (sy%status /= short_fft_success) return
        values = sx%values*conjg(sy%values)
        status = short_fft_success
    end subroutine cross_spectrogram_nd2_complex

    subroutine cross_spectrogram_nd3_real(self,x,y,values,status,p0,p1,padding,axis)
        class(short_time_fft), intent(in) :: self !! configured plan
        real(dp), intent(in) :: x(:,:,:) !! first time series batch
        real(dp), intent(in) :: y(:,:,:) !! second time series batch
        complex(dp), allocatable, intent(out) :: values(:,:,:,:) !! cross power spectrum
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: p0 !! first time slice
        integer, intent(in), optional :: p1 !! exclusive final time slice
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! time-series axis
        type(short_time_fft_result_4d) :: sx, sy

        allocate(values(0,0,0,0))
        status = short_fft_invalid_input
        if (any(shape(x) /= shape(y))) return
        call self%stft(x,sx,axis=axis,p0=p0,p1=p1,padding=padding)
        if (sx%status /= short_fft_success) return
        call self%stft(y,sy,axis=axis,p0=p0,p1=p1,padding=padding)
        if (sy%status /= short_fft_success) return
        values = sx%values*conjg(sy%values)
        status = short_fft_success
    end subroutine cross_spectrogram_nd3_real

    subroutine cross_spectrogram_nd3_complex(self,x,y,values,status,p0,p1,padding,axis)
        class(short_time_fft), intent(in) :: self !! configured plan
        complex(dp), intent(in) :: x(:,:,:) !! first time series batch
        complex(dp), intent(in) :: y(:,:,:) !! second time series batch
        complex(dp), allocatable, intent(out) :: values(:,:,:,:) !! cross power spectrum
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: p0 !! first time slice
        integer, intent(in), optional :: p1 !! exclusive final time slice
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! time-series axis
        type(short_time_fft_result_4d) :: sx, sy

        allocate(values(0,0,0,0))
        status = short_fft_invalid_input
        if (any(shape(x) /= shape(y))) return
        call self%stft(x,sx,axis=axis,p0=p0,p1=p1,padding=padding)
        if (sx%status /= short_fft_success) return
        call self%stft(y,sy,axis=axis,p0=p0,p1=p1,padding=padding)
        if (sy%status /= short_fft_success) return
        values = sx%values*conjg(sy%values)
        status = short_fft_success
    end subroutine cross_spectrogram_nd3_complex

    subroutine stft_nd4_real(self, x, result, p0, p1, k_offset, padding, axis)
        class(short_time_fft), intent(in) :: self !! configured STFT plan
        real(dp), intent(in) :: x(:,:,:,:) !! input time series
        type(short_time_fft_result_5d), intent(out) :: result !! spectrum
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: p1 !! exclusive upper time slice
        integer, intent(in), optional :: k_offset !! input sample offset
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! time-series axis, 1-based or negative
        call stft_nd4_kernel(self, cmplx(x,0.0_dp,dp), result, axis, p0, p1, &
            k_offset, padding, .false.)
    end subroutine stft_nd4_real

    subroutine stft_detrend_nd4_real(self, x, result, detrend, p0, p1, k_offset, padding, axis)
        class(short_time_fft), intent(in) :: self !! configured STFT plan
        real(dp), intent(in) :: x(:,:,:,:) !! input time series
        type(short_time_fft_result_5d), intent(out) :: result !! spectrum
        character(*), intent(in) :: detrend !! constant or linear segment detrending
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: p1 !! exclusive upper time slice
        integer, intent(in), optional :: k_offset !! input sample offset
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! time-series axis, 1-based or negative
        call stft_nd4_kernel(self, cmplx(x,0.0_dp,dp), result, axis, p0, p1, &
            k_offset, padding, .false., detrend)
    end subroutine stft_detrend_nd4_real

    subroutine stft_nd4_complex(self, x, result, p0, p1, k_offset, padding, axis)
        class(short_time_fft), intent(in) :: self !! configured STFT plan
        complex(dp), intent(in) :: x(:,:,:,:) !! input time series
        type(short_time_fft_result_5d), intent(out) :: result !! spectrum
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: p1 !! exclusive upper time slice
        integer, intent(in), optional :: k_offset !! input sample offset
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! time-series axis, 1-based or negative
        call stft_nd4_kernel(self, x, result, axis, p0, p1, &
            k_offset, padding, .true.)
    end subroutine stft_nd4_complex

    subroutine stft_detrend_nd4_complex(self, x, result, detrend, p0, p1, k_offset, padding, axis)
        class(short_time_fft), intent(in) :: self !! configured STFT plan
        complex(dp), intent(in) :: x(:,:,:,:) !! input time series
        type(short_time_fft_result_5d), intent(out) :: result !! spectrum
        character(*), intent(in) :: detrend !! constant or linear segment detrending
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: p1 !! exclusive upper time slice
        integer, intent(in), optional :: k_offset !! input sample offset
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! time-series axis, 1-based or negative
        call stft_nd4_kernel(self, x, result, axis, p0, p1, &
            k_offset, padding, .true., detrend)
    end subroutine stft_detrend_nd4_complex

    ! Flatten only the batch-index traversal: each one-dimensional signal is
    ! transformed by the same verified STFT kernel as ranks one through three.
    subroutine stft_nd4_kernel(self, x, result, axis, p0, p1, k_offset, padding, complex_input, detrend)
        class(short_time_fft), intent(in) :: self !! configured STFT plan
        complex(dp), intent(in) :: x(:,:,:,:) !! input signals
        type(short_time_fft_result_5d), intent(out) :: result !! output spectrum
        integer, intent(in), optional :: axis !! signal axis
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: p1 !! exclusive upper time slice
        integer, intent(in), optional :: k_offset !! sample offset
        character(*), intent(in), optional :: padding !! extension mode
        logical, intent(in) :: complex_input !! whether input is complex
        character(*), intent(in), optional :: detrend !! detrending mode
        type(short_time_fft_result) :: frame
        complex(dp), allocatable :: trace(:)
        integer :: a, dims(4), dout(5), idx(4), oidx(5)
        integer :: batch_axes(3), b, ib, j, k, q, t, nb

        result%status = short_fft_invalid_input
        a = fft_resolve_axis(axis,4,4)
        if (a == 0) return
        dims = shape(x)
        if (any(dims < 1)) return
        batch_axes = pack([1,2,3,4],[1,2,3,4] /= a)
        nb = product(dims(batch_axes))
        allocate(trace(dims(a)))
        do ib=0,nb-1
            b = ib
            do j=1,3
                idx(batch_axes(j)) = modulo(b,dims(batch_axes(j)))+1
                b = b/dims(batch_axes(j))
            end do
            do k=1,dims(a)
                idx(a) = k
                trace(k) = x(idx(1),idx(2),idx(3),idx(4))
            end do
            call stft_kernel(self,trace,frame,p0,p1,k_offset,padding,complex_input,detrend)
            if (frame%status /= short_fft_success) return
            if (ib == 0) then
                dout(1:4) = dims
                dout(a) = size(frame%values,1)
                dout(5) = size(frame%values,2)
                allocate(result%values(dout(1),dout(2),dout(3),dout(4),dout(5)))
                result%frequencies = frame%frequencies
                result%times = frame%times
                result%p0 = frame%p0
            end if
            do t=1,size(frame%values,2)
                do q=1,size(frame%values,1)
                    oidx(1:4) = idx
                    oidx(a) = q
                    oidx(5) = t
                    result%values(oidx(1),oidx(2),oidx(3),oidx(4),oidx(5)) = frame%values(q,t)
                end do
            end do
        end do
        result%status = short_fft_success
    end subroutine stft_nd4_kernel

    subroutine istft_nd5_real(self, s, x, status, k0, k1, p0, f_axis, t_axis)
        class(short_time_fft), intent(in) :: self !! inverse plan
        complex(dp), intent(in) :: s(:,:,:,:,:) !! STFT coefficients
        real(dp), allocatable, intent(out) :: x(:,:,:,:) !! reconstructed signal
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: k0 !! inclusive first sample
        integer, intent(in), optional :: k1 !! exclusive last sample
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: f_axis !! frequency axis of S
        integer, intent(in), optional :: t_axis !! time axis of S
        complex(dp), allocatable :: z(:,:,:,:)

        allocate(x(0,0,0,0))
        status = short_fft_invalid_input
        if (self%fft_mode /= 'onesided' .and. self%fft_mode /= 'onesided2X') return
        call inverse_nd5_kernel(self,s,z,status,k0,k1,p0,f_axis,t_axis)
        if (status == short_fft_success) x = real(z,dp)
    end subroutine istft_nd5_real

    subroutine istft_nd5_complex(self, s, x, status, k0, k1, p0, f_axis, t_axis)
        class(short_time_fft), intent(in) :: self !! inverse plan
        complex(dp), intent(in) :: s(:,:,:,:,:) !! STFT coefficients
        complex(dp), allocatable, intent(out) :: x(:,:,:,:) !! reconstructed signal
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: k0 !! inclusive first sample
        integer, intent(in), optional :: k1 !! exclusive last sample
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: f_axis !! frequency axis of S
        integer, intent(in), optional :: t_axis !! time axis of S
        complex(dp), allocatable :: z(:,:,:,:)

        allocate(x(0,0,0,0))
        status = short_fft_invalid_input
        if (self%fft_mode == 'onesided' .or. self%fft_mode == 'onesided2X') return
        call inverse_nd5_kernel(self,s,z,status,k0,k1,p0,f_axis,t_axis)
        if (status == short_fft_success) x = z
    end subroutine istft_nd5_complex

    subroutine inverse_nd5_kernel(self, s, x, status, k0, k1, p0, f_axis, t_axis)
        class(short_time_fft), intent(in) :: self !! inverse transform plan
        complex(dp), intent(in) :: s(:,:,:,:,:) !! STFT coefficients
        complex(dp), allocatable, intent(out) :: x(:,:,:,:) !! reconstructed series
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: k0 !! inclusive first sample
        integer, intent(in), optional :: k1 !! exclusive last sample
        integer, intent(in), optional :: p0 !! first time slice
        integer, intent(in), optional :: f_axis !! frequency dimension
        integer, intent(in), optional :: t_axis !! time dimension
        complex(dp), allocatable :: frame(:,:), trace(:)
        integer :: f, t, dims(5), dout(4), batch_axes(3), out_axes(4)
        integer :: idx(5), oidx(4), b, ib, j, q, p, k, nb

        allocate(x(0,0,0,0))
        status = short_fft_invalid_input
        f = fft_resolve_axis(f_axis,5,4)
        t = fft_resolve_axis(t_axis,5,5)
        if (f == 0 .or. t == 0 .or. f == t) return
        dims = shape(s)
        if (any(dims < 1)) return
        batch_axes = pack([1,2,3,4,5],[1,2,3,4,5] /= f .and. [1,2,3,4,5] /= t)
        out_axes = pack([1,2,3,4,5],[1,2,3,4,5] /= t)
        nb = product(dims(batch_axes))
        allocate(frame(dims(f),dims(t)))
        do ib=0,nb-1
            b = ib
            do j=1,3
                idx(batch_axes(j)) = modulo(b,dims(batch_axes(j)))+1
                b = b/dims(batch_axes(j))
            end do
            do p=1,dims(t)
                idx(t) = p
                do q=1,dims(f)
                    idx(f) = q
                    frame(q,p) = s(idx(1),idx(2),idx(3),idx(4),idx(5))
                end do
            end do
            call inverse_kernel(self,frame,trace,status,k0,k1,p0)
            if (status /= short_fft_success) return
            if (ib == 0) then
                dout = pack(dims,[1,2,3,4,5] /= t)
                dout(findloc(out_axes,f,dim=1)) = size(trace)
                deallocate(x)
                allocate(x(dout(1),dout(2),dout(3),dout(4)))
            end if
            do j=1,3
                oidx(findloc(out_axes,batch_axes(j),dim=1)) = idx(batch_axes(j))
            end do
            do k=1,size(trace)
                oidx(findloc(out_axes,f,dim=1)) = k
                x(oidx(1),oidx(2),oidx(3),oidx(4)) = trace(k)
            end do
        end do
        status = short_fft_success
    end subroutine inverse_nd5_kernel

    subroutine spectrogram_nd4_real(self, x, result, p0, p1, padding, axis)
        class(short_time_fft), intent(in) :: self !! transform plan
        real(dp), intent(in) :: x(:,:,:,:) !! signal array
        type(short_time_spectrogram_result_5d), intent(out) :: result !! power spectrum
        integer, intent(in), optional :: p0 !! first time slice
        integer, intent(in), optional :: p1 !! exclusive final slice
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! selected sample dimension
        type(short_time_fft_result_5d) :: spectrum

        call self%stft(x,spectrum,p0=p0,p1=p1,padding=padding,axis=axis)
        result%status = spectrum%status
        if (spectrum%status /= short_fft_success) return
        result%power = abs(spectrum%values)**2
        result%frequencies = spectrum%frequencies
        result%times = spectrum%times
        result%p0 = spectrum%p0
    end subroutine spectrogram_nd4_real

    subroutine spectrogram_nd4_complex(self, x, result, p0, p1, padding, axis)
        class(short_time_fft), intent(in) :: self !! transform plan
        complex(dp), intent(in) :: x(:,:,:,:) !! signal array
        type(short_time_spectrogram_result_5d), intent(out) :: result !! power spectrum
        integer, intent(in), optional :: p0 !! first time slice
        integer, intent(in), optional :: p1 !! exclusive final slice
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! selected sample dimension
        type(short_time_fft_result_5d) :: spectrum

        call self%stft(x,spectrum,p0=p0,p1=p1,padding=padding,axis=axis)
        result%status = spectrum%status
        if (spectrum%status /= short_fft_success) return
        result%power = abs(spectrum%values)**2
        result%frequencies = spectrum%frequencies
        result%times = spectrum%times
        result%p0 = spectrum%p0
    end subroutine spectrogram_nd4_complex

    subroutine cross_spectrogram_nd4_real(self, x, y, values, status, p0, p1, padding, axis)
        class(short_time_fft), intent(in) :: self !! transform plan
        real(dp), intent(in) :: x(:,:,:,:) !! first array of signals
        real(dp), intent(in) :: y(:,:,:,:) !! second array of signals
        complex(dp), allocatable, intent(out) :: values(:,:,:,:,:) !! cross spectrum
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: p0 !! first time slice
        integer, intent(in), optional :: p1 !! exclusive final slice
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! selected sample dimension
        type(short_time_fft_result_5d) :: sx, sy

        allocate(values(0,0,0,0,0))
        status = short_fft_invalid_input
        if (any(shape(x) /= shape(y))) return
        call self%stft(x,sx,p0=p0,p1=p1,padding=padding,axis=axis)
        if (sx%status /= short_fft_success) return
        call self%stft(y,sy,p0=p0,p1=p1,padding=padding,axis=axis)
        if (sy%status /= short_fft_success) return
        values = sx%values*conjg(sy%values)
        status = short_fft_success
    end subroutine cross_spectrogram_nd4_real

    subroutine cross_spectrogram_nd4_complex(self, x, y, values, status, p0, p1, padding, axis)
        class(short_time_fft), intent(in) :: self !! transform plan
        complex(dp), intent(in) :: x(:,:,:,:) !! first array of signals
        complex(dp), intent(in) :: y(:,:,:,:) !! second array of signals
        complex(dp), allocatable, intent(out) :: values(:,:,:,:,:) !! cross spectrum
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: p0 !! first time slice
        integer, intent(in), optional :: p1 !! exclusive final slice
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! selected sample dimension
        type(short_time_fft_result_5d) :: sx, sy

        allocate(values(0,0,0,0,0))
        status = short_fft_invalid_input
        if (any(shape(x) /= shape(y))) return
        call self%stft(x,sx,p0=p0,p1=p1,padding=padding,axis=axis)
        if (sx%status /= short_fft_success) return
        call self%stft(y,sy,p0=p0,p1=p1,padding=padding,axis=axis)
        if (sy%status /= short_fft_success) return
        values = sx%values*conjg(sy%values)
        status = short_fft_success
    end subroutine cross_spectrogram_nd4_complex

    subroutine stft_nd5_real(self, x, result, p0, p1, k_offset, padding, axis)
        class(short_time_fft), intent(in) :: self !! configured STFT plan
        real(dp), intent(in) :: x(:,:,:,:,:) !! input time series
        type(short_time_fft_result_6d), intent(out) :: result !! spectrum
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: p1 !! exclusive upper time slice
        integer, intent(in), optional :: k_offset !! input sample offset
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! time-series axis, 1-based or negative
        call stft_nd5_kernel(self, cmplx(x,0.0_dp,dp), result, axis, p0, p1, &
            k_offset, padding, .false.)
    end subroutine stft_nd5_real

    subroutine stft_detrend_nd5_real(self, x, result, detrend, p0, p1, k_offset, padding, axis)
        class(short_time_fft), intent(in) :: self !! configured STFT plan
        real(dp), intent(in) :: x(:,:,:,:,:) !! input time series
        type(short_time_fft_result_6d), intent(out) :: result !! spectrum
        character(*), intent(in) :: detrend !! constant or linear segment detrending
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: p1 !! exclusive upper time slice
        integer, intent(in), optional :: k_offset !! input sample offset
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! time-series axis, 1-based or negative
        call stft_nd5_kernel(self, cmplx(x,0.0_dp,dp), result, axis, p0, p1, &
            k_offset, padding, .false., detrend)
    end subroutine stft_detrend_nd5_real

    subroutine stft_nd5_complex(self, x, result, p0, p1, k_offset, padding, axis)
        class(short_time_fft), intent(in) :: self !! configured STFT plan
        complex(dp), intent(in) :: x(:,:,:,:,:) !! input time series
        type(short_time_fft_result_6d), intent(out) :: result !! spectrum
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: p1 !! exclusive upper time slice
        integer, intent(in), optional :: k_offset !! input sample offset
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! time-series axis, 1-based or negative
        call stft_nd5_kernel(self, x, result, axis, p0, p1, &
            k_offset, padding, .true.)
    end subroutine stft_nd5_complex

    subroutine stft_detrend_nd5_complex(self, x, result, detrend, p0, p1, k_offset, padding, axis)
        class(short_time_fft), intent(in) :: self !! configured STFT plan
        complex(dp), intent(in) :: x(:,:,:,:,:) !! input time series
        type(short_time_fft_result_6d), intent(out) :: result !! spectrum
        character(*), intent(in) :: detrend !! constant or linear segment detrending
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: p1 !! exclusive upper time slice
        integer, intent(in), optional :: k_offset !! input sample offset
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! time-series axis, 1-based or negative
        call stft_nd5_kernel(self, x, result, axis, p0, p1, &
            k_offset, padding, .true., detrend)
    end subroutine stft_detrend_nd5_complex

    ! Flatten only the batch-index traversal: each one-dimensional signal is
    ! transformed by the same verified STFT kernel as ranks one through three.
    subroutine stft_nd5_kernel(self, x, result, axis, p0, p1, k_offset, padding, complex_input, detrend)
        class(short_time_fft), intent(in) :: self !! configured STFT plan
        complex(dp), intent(in) :: x(:,:,:,:,:) !! input signals
        type(short_time_fft_result_6d), intent(out) :: result !! output spectrum
        integer, intent(in), optional :: axis !! signal axis
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: p1 !! exclusive upper time slice
        integer, intent(in), optional :: k_offset !! sample offset
        character(*), intent(in), optional :: padding !! extension mode
        logical, intent(in) :: complex_input !! whether input is complex
        character(*), intent(in), optional :: detrend !! detrending mode
        type(short_time_fft_result) :: frame
        complex(dp), allocatable :: trace(:)
        integer :: a, dims(5), dout(6), idx(5), oidx(6)
        integer :: batch_axes(4), b, ib, j, k, q, t, nb

        result%status = short_fft_invalid_input
        a = fft_resolve_axis(axis,5,5)
        if (a == 0) return
        dims = shape(x)
        if (any(dims < 1)) return
        batch_axes = pack([1,2,3,4,5],[1,2,3,4,5] /= a)
        nb = product(dims(batch_axes))
        allocate(trace(dims(a)))
        do ib=0,nb-1
            b = ib
            do j=1,4
                idx(batch_axes(j)) = modulo(b,dims(batch_axes(j)))+1
                b = b/dims(batch_axes(j))
            end do
            do k=1,dims(a)
                idx(a) = k
                trace(k) = x(idx(1),idx(2),idx(3),idx(4),idx(5))
            end do
            call stft_kernel(self,trace,frame,p0,p1,k_offset,padding,complex_input,detrend)
            if (frame%status /= short_fft_success) return
            if (ib == 0) then
                dout(1:5) = dims
                dout(a) = size(frame%values,1)
                dout(6) = size(frame%values,2)
                allocate(result%values(dout(1),dout(2),dout(3),dout(4),dout(5),dout(6)))
                result%frequencies = frame%frequencies
                result%times = frame%times
                result%p0 = frame%p0
            end if
            do t=1,size(frame%values,2)
                do q=1,size(frame%values,1)
                    oidx(1:5) = idx
                    oidx(a) = q
                    oidx(6) = t
                    result%values(oidx(1),oidx(2),oidx(3),oidx(4),oidx(5),oidx(6)) = frame%values(q,t)
                end do
            end do
        end do
        result%status = short_fft_success
    end subroutine stft_nd5_kernel

    subroutine istft_nd6_real(self, s, x, status, k0, k1, p0, f_axis, t_axis)
        class(short_time_fft), intent(in) :: self !! inverse plan
        complex(dp), intent(in) :: s(:,:,:,:,:,:) !! STFT coefficients
        real(dp), allocatable, intent(out) :: x(:,:,:,:,:) !! reconstructed signal
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: k0 !! inclusive first sample
        integer, intent(in), optional :: k1 !! exclusive last sample
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: f_axis !! frequency axis of S
        integer, intent(in), optional :: t_axis !! time axis of S
        complex(dp), allocatable :: z(:,:,:,:,:)

        allocate(x(0,0,0,0,0))
        status = short_fft_invalid_input
        if (self%fft_mode /= 'onesided' .and. self%fft_mode /= 'onesided2X') return
        call inverse_nd6_kernel(self,s,z,status,k0,k1,p0,f_axis,t_axis)
        if (status == short_fft_success) x = real(z,dp)
    end subroutine istft_nd6_real

    subroutine istft_nd6_complex(self, s, x, status, k0, k1, p0, f_axis, t_axis)
        class(short_time_fft), intent(in) :: self !! inverse plan
        complex(dp), intent(in) :: s(:,:,:,:,:,:) !! STFT coefficients
        complex(dp), allocatable, intent(out) :: x(:,:,:,:,:) !! reconstructed signal
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: k0 !! inclusive first sample
        integer, intent(in), optional :: k1 !! exclusive last sample
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: f_axis !! frequency axis of S
        integer, intent(in), optional :: t_axis !! time axis of S
        complex(dp), allocatable :: z(:,:,:,:,:)

        allocate(x(0,0,0,0,0))
        status = short_fft_invalid_input
        if (self%fft_mode == 'onesided' .or. self%fft_mode == 'onesided2X') return
        call inverse_nd6_kernel(self,s,z,status,k0,k1,p0,f_axis,t_axis)
        if (status == short_fft_success) x = z
    end subroutine istft_nd6_complex

    subroutine inverse_nd6_kernel(self, s, x, status, k0, k1, p0, f_axis, t_axis)
        class(short_time_fft), intent(in) :: self !! inverse transform plan
        complex(dp), intent(in) :: s(:,:,:,:,:,:) !! STFT coefficients
        complex(dp), allocatable, intent(out) :: x(:,:,:,:,:) !! reconstructed series
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: k0 !! inclusive first sample
        integer, intent(in), optional :: k1 !! exclusive last sample
        integer, intent(in), optional :: p0 !! first time slice
        integer, intent(in), optional :: f_axis !! frequency dimension
        integer, intent(in), optional :: t_axis !! time dimension
        complex(dp), allocatable :: frame(:,:), trace(:)
        integer :: f, t, dims(6), dout(5), batch_axes(4), out_axes(5)
        integer :: idx(6), oidx(5), b, ib, j, q, p, k, nb

        allocate(x(0,0,0,0,0))
        status = short_fft_invalid_input
        f = fft_resolve_axis(f_axis,6,5)
        t = fft_resolve_axis(t_axis,6,6)
        if (f == 0 .or. t == 0 .or. f == t) return
        dims = shape(s)
        if (any(dims < 1)) return
        batch_axes = pack([1,2,3,4,5,6],[1,2,3,4,5,6] /= f .and. [1,2,3,4,5,6] /= t)
        out_axes = pack([1,2,3,4,5,6],[1,2,3,4,5,6] /= t)
        nb = product(dims(batch_axes))
        allocate(frame(dims(f),dims(t)))
        do ib=0,nb-1
            b = ib
            do j=1,4
                idx(batch_axes(j)) = modulo(b,dims(batch_axes(j)))+1
                b = b/dims(batch_axes(j))
            end do
            do p=1,dims(t)
                idx(t) = p
                do q=1,dims(f)
                    idx(f) = q
                    frame(q,p) = s(idx(1),idx(2),idx(3),idx(4),idx(5),idx(6))
                end do
            end do
            call inverse_kernel(self,frame,trace,status,k0,k1,p0)
            if (status /= short_fft_success) return
            if (ib == 0) then
                dout = pack(dims,[1,2,3,4,5,6] /= t)
                dout(findloc(out_axes,f,dim=1)) = size(trace)
                deallocate(x)
                allocate(x(dout(1),dout(2),dout(3),dout(4),dout(5)))
            end if
            do j=1,4
                oidx(findloc(out_axes,batch_axes(j),dim=1)) = idx(batch_axes(j))
            end do
            do k=1,size(trace)
                oidx(findloc(out_axes,f,dim=1)) = k
                x(oidx(1),oidx(2),oidx(3),oidx(4),oidx(5)) = trace(k)
            end do
        end do
        status = short_fft_success
    end subroutine inverse_nd6_kernel

    subroutine spectrogram_nd5_real(self, x, result, p0, p1, padding, axis)
        class(short_time_fft), intent(in) :: self !! transform plan
        real(dp), intent(in) :: x(:,:,:,:,:) !! signal array
        type(short_time_spectrogram_result_6d), intent(out) :: result !! power spectrum
        integer, intent(in), optional :: p0 !! first time slice
        integer, intent(in), optional :: p1 !! exclusive final slice
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! selected sample dimension
        type(short_time_fft_result_6d) :: spectrum

        call self%stft(x,spectrum,p0=p0,p1=p1,padding=padding,axis=axis)
        result%status = spectrum%status
        if (spectrum%status /= short_fft_success) return
        result%power = abs(spectrum%values)**2
        result%frequencies = spectrum%frequencies
        result%times = spectrum%times
        result%p0 = spectrum%p0
    end subroutine spectrogram_nd5_real

    subroutine spectrogram_nd5_complex(self, x, result, p0, p1, padding, axis)
        class(short_time_fft), intent(in) :: self !! transform plan
        complex(dp), intent(in) :: x(:,:,:,:,:) !! signal array
        type(short_time_spectrogram_result_6d), intent(out) :: result !! power spectrum
        integer, intent(in), optional :: p0 !! first time slice
        integer, intent(in), optional :: p1 !! exclusive final slice
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! selected sample dimension
        type(short_time_fft_result_6d) :: spectrum

        call self%stft(x,spectrum,p0=p0,p1=p1,padding=padding,axis=axis)
        result%status = spectrum%status
        if (spectrum%status /= short_fft_success) return
        result%power = abs(spectrum%values)**2
        result%frequencies = spectrum%frequencies
        result%times = spectrum%times
        result%p0 = spectrum%p0
    end subroutine spectrogram_nd5_complex

    subroutine cross_spectrogram_nd5_real(self, x, y, values, status, p0, p1, padding, axis)
        class(short_time_fft), intent(in) :: self !! transform plan
        real(dp), intent(in) :: x(:,:,:,:,:) !! first array of signals
        real(dp), intent(in) :: y(:,:,:,:,:) !! second array of signals
        complex(dp), allocatable, intent(out) :: values(:,:,:,:,:,:) !! cross spectrum
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: p0 !! first time slice
        integer, intent(in), optional :: p1 !! exclusive final slice
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! selected sample dimension
        type(short_time_fft_result_6d) :: sx, sy

        allocate(values(0,0,0,0,0,0))
        status = short_fft_invalid_input
        if (any(shape(x) /= shape(y))) return
        call self%stft(x,sx,p0=p0,p1=p1,padding=padding,axis=axis)
        if (sx%status /= short_fft_success) return
        call self%stft(y,sy,p0=p0,p1=p1,padding=padding,axis=axis)
        if (sy%status /= short_fft_success) return
        values = sx%values*conjg(sy%values)
        status = short_fft_success
    end subroutine cross_spectrogram_nd5_real

    subroutine cross_spectrogram_nd5_complex(self, x, y, values, status, p0, p1, padding, axis)
        class(short_time_fft), intent(in) :: self !! transform plan
        complex(dp), intent(in) :: x(:,:,:,:,:) !! first array of signals
        complex(dp), intent(in) :: y(:,:,:,:,:) !! second array of signals
        complex(dp), allocatable, intent(out) :: values(:,:,:,:,:,:) !! cross spectrum
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: p0 !! first time slice
        integer, intent(in), optional :: p1 !! exclusive final slice
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! selected sample dimension
        type(short_time_fft_result_6d) :: sx, sy

        allocate(values(0,0,0,0,0,0))
        status = short_fft_invalid_input
        if (any(shape(x) /= shape(y))) return
        call self%stft(x,sx,p0=p0,p1=p1,padding=padding,axis=axis)
        if (sx%status /= short_fft_success) return
        call self%stft(y,sy,p0=p0,p1=p1,padding=padding,axis=axis)
        if (sy%status /= short_fft_success) return
        values = sx%values*conjg(sy%values)
        status = short_fft_success
    end subroutine cross_spectrogram_nd5_complex

    subroutine stft_nd6_real(self, x, result, p0, p1, k_offset, padding, axis)
        class(short_time_fft), intent(in) :: self !! configured STFT plan
        real(dp), intent(in) :: x(:,:,:,:,:,:) !! input time series
        type(short_time_fft_result_7d), intent(out) :: result !! spectrum
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: p1 !! exclusive upper time slice
        integer, intent(in), optional :: k_offset !! input sample offset
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! time-series axis, 1-based or negative
        call stft_nd6_kernel(self, cmplx(x,0.0_dp,dp), result, axis, p0, p1, &
            k_offset, padding, .false.)
    end subroutine stft_nd6_real

    subroutine stft_detrend_nd6_real(self, x, result, detrend, p0, p1, k_offset, padding, axis)
        class(short_time_fft), intent(in) :: self !! configured STFT plan
        real(dp), intent(in) :: x(:,:,:,:,:,:) !! input time series
        type(short_time_fft_result_7d), intent(out) :: result !! spectrum
        character(*), intent(in) :: detrend !! constant or linear segment detrending
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: p1 !! exclusive upper time slice
        integer, intent(in), optional :: k_offset !! input sample offset
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! time-series axis, 1-based or negative
        call stft_nd6_kernel(self, cmplx(x,0.0_dp,dp), result, axis, p0, p1, &
            k_offset, padding, .false., detrend)
    end subroutine stft_detrend_nd6_real

    subroutine stft_nd6_complex(self, x, result, p0, p1, k_offset, padding, axis)
        class(short_time_fft), intent(in) :: self !! configured STFT plan
        complex(dp), intent(in) :: x(:,:,:,:,:,:) !! input time series
        type(short_time_fft_result_7d), intent(out) :: result !! spectrum
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: p1 !! exclusive upper time slice
        integer, intent(in), optional :: k_offset !! input sample offset
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! time-series axis, 1-based or negative
        call stft_nd6_kernel(self, x, result, axis, p0, p1, &
            k_offset, padding, .true.)
    end subroutine stft_nd6_complex

    subroutine stft_detrend_nd6_complex(self, x, result, detrend, p0, p1, k_offset, padding, axis)
        class(short_time_fft), intent(in) :: self !! configured STFT plan
        complex(dp), intent(in) :: x(:,:,:,:,:,:) !! input time series
        type(short_time_fft_result_7d), intent(out) :: result !! spectrum
        character(*), intent(in) :: detrend !! constant or linear segment detrending
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: p1 !! exclusive upper time slice
        integer, intent(in), optional :: k_offset !! input sample offset
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! time-series axis, 1-based or negative
        call stft_nd6_kernel(self, x, result, axis, p0, p1, &
            k_offset, padding, .true., detrend)
    end subroutine stft_detrend_nd6_complex

    ! Flatten only the batch-index traversal: each one-dimensional signal is
    ! transformed by the same verified STFT kernel as ranks one through three.
    subroutine stft_nd6_kernel(self, x, result, axis, p0, p1, k_offset, padding, complex_input, detrend)
        class(short_time_fft), intent(in) :: self !! configured STFT plan
        complex(dp), intent(in) :: x(:,:,:,:,:,:) !! input signals
        type(short_time_fft_result_7d), intent(out) :: result !! output spectrum
        integer, intent(in), optional :: axis !! signal axis
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: p1 !! exclusive upper time slice
        integer, intent(in), optional :: k_offset !! sample offset
        character(*), intent(in), optional :: padding !! extension mode
        logical, intent(in) :: complex_input !! whether input is complex
        character(*), intent(in), optional :: detrend !! detrending mode
        type(short_time_fft_result) :: frame
        complex(dp), allocatable :: trace(:)
        integer :: a, dims(6), dout(7), idx(6), oidx(7)
        integer :: batch_axes(5), b, ib, j, k, q, t, nb

        result%status = short_fft_invalid_input
        a = fft_resolve_axis(axis,6,6)
        if (a == 0) return
        dims = shape(x)
        if (any(dims < 1)) return
        batch_axes = pack([1,2,3,4,5,6],[1,2,3,4,5,6] /= a)
        nb = product(dims(batch_axes))
        allocate(trace(dims(a)))
        do ib=0,nb-1
            b = ib
            do j=1,5
                idx(batch_axes(j)) = modulo(b,dims(batch_axes(j)))+1
                b = b/dims(batch_axes(j))
            end do
            do k=1,dims(a)
                idx(a) = k
                trace(k) = x(idx(1),idx(2),idx(3),idx(4),idx(5),idx(6))
            end do
            call stft_kernel(self,trace,frame,p0,p1,k_offset,padding,complex_input,detrend)
            if (frame%status /= short_fft_success) return
            if (ib == 0) then
                dout(1:6) = dims
                dout(a) = size(frame%values,1)
                dout(7) = size(frame%values,2)
                allocate(result%values(dout(1),dout(2),dout(3),dout(4),dout(5),dout(6),dout(7)))
                result%frequencies = frame%frequencies
                result%times = frame%times
                result%p0 = frame%p0
            end if
            do t=1,size(frame%values,2)
                do q=1,size(frame%values,1)
                    oidx(1:6) = idx
                    oidx(a) = q
                    oidx(7) = t
                    result%values(oidx(1),oidx(2),oidx(3),oidx(4),oidx(5),oidx(6),oidx(7)) = frame%values(q,t)
                end do
            end do
        end do
        result%status = short_fft_success
    end subroutine stft_nd6_kernel

    subroutine istft_nd7_real(self, s, x, status, k0, k1, p0, f_axis, t_axis)
        class(short_time_fft), intent(in) :: self !! inverse plan
        complex(dp), intent(in) :: s(:,:,:,:,:,:,:) !! STFT coefficients
        real(dp), allocatable, intent(out) :: x(:,:,:,:,:,:) !! reconstructed signal
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: k0 !! inclusive first sample
        integer, intent(in), optional :: k1 !! exclusive last sample
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: f_axis !! frequency axis of S
        integer, intent(in), optional :: t_axis !! time axis of S
        complex(dp), allocatable :: z(:,:,:,:,:,:)

        allocate(x(0,0,0,0,0,0))
        status = short_fft_invalid_input
        if (self%fft_mode /= 'onesided' .and. self%fft_mode /= 'onesided2X') return
        call inverse_nd7_kernel(self,s,z,status,k0,k1,p0,f_axis,t_axis)
        if (status == short_fft_success) x = real(z,dp)
    end subroutine istft_nd7_real

    subroutine istft_nd7_complex(self, s, x, status, k0, k1, p0, f_axis, t_axis)
        class(short_time_fft), intent(in) :: self !! inverse plan
        complex(dp), intent(in) :: s(:,:,:,:,:,:,:) !! STFT coefficients
        complex(dp), allocatable, intent(out) :: x(:,:,:,:,:,:) !! reconstructed signal
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: k0 !! inclusive first sample
        integer, intent(in), optional :: k1 !! exclusive last sample
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: f_axis !! frequency axis of S
        integer, intent(in), optional :: t_axis !! time axis of S
        complex(dp), allocatable :: z(:,:,:,:,:,:)

        allocate(x(0,0,0,0,0,0))
        status = short_fft_invalid_input
        if (self%fft_mode == 'onesided' .or. self%fft_mode == 'onesided2X') return
        call inverse_nd7_kernel(self,s,z,status,k0,k1,p0,f_axis,t_axis)
        if (status == short_fft_success) x = z
    end subroutine istft_nd7_complex

    subroutine inverse_nd7_kernel(self, s, x, status, k0, k1, p0, f_axis, t_axis)
        class(short_time_fft), intent(in) :: self !! inverse transform plan
        complex(dp), intent(in) :: s(:,:,:,:,:,:,:) !! STFT coefficients
        complex(dp), allocatable, intent(out) :: x(:,:,:,:,:,:) !! reconstructed series
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: k0 !! inclusive first sample
        integer, intent(in), optional :: k1 !! exclusive last sample
        integer, intent(in), optional :: p0 !! first time slice
        integer, intent(in), optional :: f_axis !! frequency dimension
        integer, intent(in), optional :: t_axis !! time dimension
        complex(dp), allocatable :: frame(:,:), trace(:)
        integer :: f, t, dims(7), dout(6), batch_axes(5), out_axes(6)
        integer :: idx(7), oidx(6), b, ib, j, q, p, k, nb

        allocate(x(0,0,0,0,0,0))
        status = short_fft_invalid_input
        f = fft_resolve_axis(f_axis,7,6)
        t = fft_resolve_axis(t_axis,7,7)
        if (f == 0 .or. t == 0 .or. f == t) return
        dims = shape(s)
        if (any(dims < 1)) return
        batch_axes = pack([1,2,3,4,5,6,7],[1,2,3,4,5,6,7] /= f .and. [1,2,3,4,5,6,7] /= t)
        out_axes = pack([1,2,3,4,5,6,7],[1,2,3,4,5,6,7] /= t)
        nb = product(dims(batch_axes))
        allocate(frame(dims(f),dims(t)))
        do ib=0,nb-1
            b = ib
            do j=1,5
                idx(batch_axes(j)) = modulo(b,dims(batch_axes(j)))+1
                b = b/dims(batch_axes(j))
            end do
            do p=1,dims(t)
                idx(t) = p
                do q=1,dims(f)
                    idx(f) = q
                    frame(q,p) = s(idx(1),idx(2),idx(3),idx(4),idx(5),idx(6),idx(7))
                end do
            end do
            call inverse_kernel(self,frame,trace,status,k0,k1,p0)
            if (status /= short_fft_success) return
            if (ib == 0) then
                dout = pack(dims,[1,2,3,4,5,6,7] /= t)
                dout(findloc(out_axes,f,dim=1)) = size(trace)
                deallocate(x)
                allocate(x(dout(1),dout(2),dout(3),dout(4),dout(5),dout(6)))
            end if
            do j=1,5
                oidx(findloc(out_axes,batch_axes(j),dim=1)) = idx(batch_axes(j))
            end do
            do k=1,size(trace)
                oidx(findloc(out_axes,f,dim=1)) = k
                x(oidx(1),oidx(2),oidx(3),oidx(4),oidx(5),oidx(6)) = trace(k)
            end do
        end do
        status = short_fft_success
    end subroutine inverse_nd7_kernel

    subroutine spectrogram_nd6_real(self, x, result, p0, p1, padding, axis)
        class(short_time_fft), intent(in) :: self !! transform plan
        real(dp), intent(in) :: x(:,:,:,:,:,:) !! signal array
        type(short_time_spectrogram_result_7d), intent(out) :: result !! power spectrum
        integer, intent(in), optional :: p0 !! first time slice
        integer, intent(in), optional :: p1 !! exclusive final slice
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! selected sample dimension
        type(short_time_fft_result_7d) :: spectrum

        call self%stft(x,spectrum,p0=p0,p1=p1,padding=padding,axis=axis)
        result%status = spectrum%status
        if (spectrum%status /= short_fft_success) return
        result%power = abs(spectrum%values)**2
        result%frequencies = spectrum%frequencies
        result%times = spectrum%times
        result%p0 = spectrum%p0
    end subroutine spectrogram_nd6_real

    subroutine spectrogram_nd6_complex(self, x, result, p0, p1, padding, axis)
        class(short_time_fft), intent(in) :: self !! transform plan
        complex(dp), intent(in) :: x(:,:,:,:,:,:) !! signal array
        type(short_time_spectrogram_result_7d), intent(out) :: result !! power spectrum
        integer, intent(in), optional :: p0 !! first time slice
        integer, intent(in), optional :: p1 !! exclusive final slice
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! selected sample dimension
        type(short_time_fft_result_7d) :: spectrum

        call self%stft(x,spectrum,p0=p0,p1=p1,padding=padding,axis=axis)
        result%status = spectrum%status
        if (spectrum%status /= short_fft_success) return
        result%power = abs(spectrum%values)**2
        result%frequencies = spectrum%frequencies
        result%times = spectrum%times
        result%p0 = spectrum%p0
    end subroutine spectrogram_nd6_complex

    subroutine cross_spectrogram_nd6_real(self, x, y, values, status, p0, p1, padding, axis)
        class(short_time_fft), intent(in) :: self !! transform plan
        real(dp), intent(in) :: x(:,:,:,:,:,:) !! first array of signals
        real(dp), intent(in) :: y(:,:,:,:,:,:) !! second array of signals
        complex(dp), allocatable, intent(out) :: values(:,:,:,:,:,:,:) !! cross spectrum
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: p0 !! first time slice
        integer, intent(in), optional :: p1 !! exclusive final slice
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! selected sample dimension
        type(short_time_fft_result_7d) :: sx, sy

        allocate(values(0,0,0,0,0,0,0))
        status = short_fft_invalid_input
        if (any(shape(x) /= shape(y))) return
        call self%stft(x,sx,p0=p0,p1=p1,padding=padding,axis=axis)
        if (sx%status /= short_fft_success) return
        call self%stft(y,sy,p0=p0,p1=p1,padding=padding,axis=axis)
        if (sy%status /= short_fft_success) return
        values = sx%values*conjg(sy%values)
        status = short_fft_success
    end subroutine cross_spectrogram_nd6_real

    subroutine cross_spectrogram_nd6_complex(self, x, y, values, status, p0, p1, padding, axis)
        class(short_time_fft), intent(in) :: self !! transform plan
        complex(dp), intent(in) :: x(:,:,:,:,:,:) !! first array of signals
        complex(dp), intent(in) :: y(:,:,:,:,:,:) !! second array of signals
        complex(dp), allocatable, intent(out) :: values(:,:,:,:,:,:,:) !! cross spectrum
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: p0 !! first time slice
        integer, intent(in), optional :: p1 !! exclusive final slice
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! selected sample dimension
        type(short_time_fft_result_7d) :: sx, sy

        allocate(values(0,0,0,0,0,0,0))
        status = short_fft_invalid_input
        if (any(shape(x) /= shape(y))) return
        call self%stft(x,sx,p0=p0,p1=p1,padding=padding,axis=axis)
        if (sx%status /= short_fft_success) return
        call self%stft(y,sy,p0=p0,p1=p1,padding=padding,axis=axis)
        if (sy%status /= short_fft_success) return
        values = sx%values*conjg(sy%values)
        status = short_fft_success
    end subroutine cross_spectrogram_nd6_complex

    subroutine stft_nd7_real(self, x, result, p0, p1, k_offset, padding, axis)
        class(short_time_fft), intent(in) :: self !! configured STFT plan
        real(dp), intent(in) :: x(:,:,:,:,:,:,:) !! input time series
        type(short_time_fft_result_8d), intent(out) :: result !! spectrum
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: p1 !! exclusive upper time slice
        integer, intent(in), optional :: k_offset !! input sample offset
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! time-series axis, 1-based or negative
        call stft_nd7_kernel(self, cmplx(x,0.0_dp,dp), result, axis, p0, p1, &
            k_offset, padding, .false.)
    end subroutine stft_nd7_real

    subroutine stft_detrend_nd7_real(self, x, result, detrend, p0, p1, k_offset, padding, axis)
        class(short_time_fft), intent(in) :: self !! configured STFT plan
        real(dp), intent(in) :: x(:,:,:,:,:,:,:) !! input time series
        type(short_time_fft_result_8d), intent(out) :: result !! spectrum
        character(*), intent(in) :: detrend !! constant or linear segment detrending
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: p1 !! exclusive upper time slice
        integer, intent(in), optional :: k_offset !! input sample offset
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! time-series axis, 1-based or negative
        call stft_nd7_kernel(self, cmplx(x,0.0_dp,dp), result, axis, p0, p1, &
            k_offset, padding, .false., detrend)
    end subroutine stft_detrend_nd7_real

    subroutine stft_nd7_complex(self, x, result, p0, p1, k_offset, padding, axis)
        class(short_time_fft), intent(in) :: self !! configured STFT plan
        complex(dp), intent(in) :: x(:,:,:,:,:,:,:) !! input time series
        type(short_time_fft_result_8d), intent(out) :: result !! spectrum
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: p1 !! exclusive upper time slice
        integer, intent(in), optional :: k_offset !! input sample offset
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! time-series axis, 1-based or negative
        call stft_nd7_kernel(self, x, result, axis, p0, p1, &
            k_offset, padding, .true.)
    end subroutine stft_nd7_complex

    subroutine stft_detrend_nd7_complex(self, x, result, detrend, p0, p1, k_offset, padding, axis)
        class(short_time_fft), intent(in) :: self !! configured STFT plan
        complex(dp), intent(in) :: x(:,:,:,:,:,:,:) !! input time series
        type(short_time_fft_result_8d), intent(out) :: result !! spectrum
        character(*), intent(in) :: detrend !! constant or linear segment detrending
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: p1 !! exclusive upper time slice
        integer, intent(in), optional :: k_offset !! input sample offset
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! time-series axis, 1-based or negative
        call stft_nd7_kernel(self, x, result, axis, p0, p1, &
            k_offset, padding, .true., detrend)
    end subroutine stft_detrend_nd7_complex

    ! Flatten only the batch-index traversal: each one-dimensional signal is
    ! transformed by the same verified STFT kernel as ranks one through three.
    subroutine stft_nd7_kernel(self, x, result, axis, p0, p1, k_offset, padding, complex_input, detrend)
        class(short_time_fft), intent(in) :: self !! configured STFT plan
        complex(dp), intent(in) :: x(:,:,:,:,:,:,:) !! input signals
        type(short_time_fft_result_8d), intent(out) :: result !! output spectrum
        integer, intent(in), optional :: axis !! signal axis
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: p1 !! exclusive upper time slice
        integer, intent(in), optional :: k_offset !! sample offset
        character(*), intent(in), optional :: padding !! extension mode
        logical, intent(in) :: complex_input !! whether input is complex
        character(*), intent(in), optional :: detrend !! detrending mode
        type(short_time_fft_result) :: frame
        complex(dp), allocatable :: trace(:)
        integer :: a, dims(7), dout(8), idx(7), oidx(8)
        integer :: batch_axes(6), b, ib, j, k, q, t, nb

        result%status = short_fft_invalid_input
        a = fft_resolve_axis(axis,7,7)
        if (a == 0) return
        dims = shape(x)
        if (any(dims < 1)) return
        batch_axes = pack([1,2,3,4,5,6,7],[1,2,3,4,5,6,7] /= a)
        nb = product(dims(batch_axes))
        allocate(trace(dims(a)))
        do ib=0,nb-1
            b = ib
            do j=1,6
                idx(batch_axes(j)) = modulo(b,dims(batch_axes(j)))+1
                b = b/dims(batch_axes(j))
            end do
            do k=1,dims(a)
                idx(a) = k
                trace(k) = x(idx(1),idx(2),idx(3),idx(4),idx(5),idx(6),idx(7))
            end do
            call stft_kernel(self,trace,frame,p0,p1,k_offset,padding,complex_input,detrend)
            if (frame%status /= short_fft_success) return
            if (ib == 0) then
                dout(1:7) = dims
                dout(a) = size(frame%values,1)
                dout(8) = size(frame%values,2)
                allocate(result%values(dout(1),dout(2),dout(3),dout(4),dout(5),dout(6),dout(7),dout(8)))
                result%frequencies = frame%frequencies
                result%times = frame%times
                result%p0 = frame%p0
            end if
            do t=1,size(frame%values,2)
                do q=1,size(frame%values,1)
                    oidx(1:7) = idx
                    oidx(a) = q
                    oidx(8) = t
                    result%values(oidx(1),oidx(2),oidx(3),oidx(4),oidx(5),oidx(6),oidx(7),oidx(8)) = frame%values(q,t)
                end do
            end do
        end do
        result%status = short_fft_success
    end subroutine stft_nd7_kernel

    subroutine istft_nd8_real(self, s, x, status, k0, k1, p0, f_axis, t_axis)
        class(short_time_fft), intent(in) :: self !! inverse plan
        complex(dp), intent(in) :: s(:,:,:,:,:,:,:,:) !! STFT coefficients
        real(dp), allocatable, intent(out) :: x(:,:,:,:,:,:,:) !! reconstructed signal
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: k0 !! inclusive first sample
        integer, intent(in), optional :: k1 !! exclusive last sample
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: f_axis !! frequency axis of S
        integer, intent(in), optional :: t_axis !! time axis of S
        complex(dp), allocatable :: z(:,:,:,:,:,:,:)

        allocate(x(0,0,0,0,0,0,0))
        status = short_fft_invalid_input
        if (self%fft_mode /= 'onesided' .and. self%fft_mode /= 'onesided2X') return
        call inverse_nd8_kernel(self,s,z,status,k0,k1,p0,f_axis,t_axis)
        if (status == short_fft_success) x = real(z,dp)
    end subroutine istft_nd8_real

    subroutine istft_nd8_complex(self, s, x, status, k0, k1, p0, f_axis, t_axis)
        class(short_time_fft), intent(in) :: self !! inverse plan
        complex(dp), intent(in) :: s(:,:,:,:,:,:,:,:) !! STFT coefficients
        complex(dp), allocatable, intent(out) :: x(:,:,:,:,:,:,:) !! reconstructed signal
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: k0 !! inclusive first sample
        integer, intent(in), optional :: k1 !! exclusive last sample
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: f_axis !! frequency axis of S
        integer, intent(in), optional :: t_axis !! time axis of S
        complex(dp), allocatable :: z(:,:,:,:,:,:,:)

        allocate(x(0,0,0,0,0,0,0))
        status = short_fft_invalid_input
        if (self%fft_mode == 'onesided' .or. self%fft_mode == 'onesided2X') return
        call inverse_nd8_kernel(self,s,z,status,k0,k1,p0,f_axis,t_axis)
        if (status == short_fft_success) x = z
    end subroutine istft_nd8_complex

    subroutine inverse_nd8_kernel(self, s, x, status, k0, k1, p0, f_axis, t_axis)
        class(short_time_fft), intent(in) :: self !! inverse transform plan
        complex(dp), intent(in) :: s(:,:,:,:,:,:,:,:) !! STFT coefficients
        complex(dp), allocatable, intent(out) :: x(:,:,:,:,:,:,:) !! reconstructed series
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: k0 !! inclusive first sample
        integer, intent(in), optional :: k1 !! exclusive last sample
        integer, intent(in), optional :: p0 !! first time slice
        integer, intent(in), optional :: f_axis !! frequency dimension
        integer, intent(in), optional :: t_axis !! time dimension
        complex(dp), allocatable :: frame(:,:), trace(:)
        integer :: f, t, dims(8), dout(7), batch_axes(6), out_axes(7)
        integer :: idx(8), oidx(7), b, ib, j, q, p, k, nb

        allocate(x(0,0,0,0,0,0,0))
        status = short_fft_invalid_input
        f = fft_resolve_axis(f_axis,8,7)
        t = fft_resolve_axis(t_axis,8,8)
        if (f == 0 .or. t == 0 .or. f == t) return
        dims = shape(s)
        if (any(dims < 1)) return
        batch_axes = pack([1,2,3,4,5,6,7,8],[1,2,3,4,5,6,7,8] /= f .and. [1,2,3,4,5,6,7,8] /= t)
        out_axes = pack([1,2,3,4,5,6,7,8],[1,2,3,4,5,6,7,8] /= t)
        nb = product(dims(batch_axes))
        allocate(frame(dims(f),dims(t)))
        do ib=0,nb-1
            b = ib
            do j=1,6
                idx(batch_axes(j)) = modulo(b,dims(batch_axes(j)))+1
                b = b/dims(batch_axes(j))
            end do
            do p=1,dims(t)
                idx(t) = p
                do q=1,dims(f)
                    idx(f) = q
                    frame(q,p) = s(idx(1),idx(2),idx(3),idx(4),idx(5),idx(6),idx(7),idx(8))
                end do
            end do
            call inverse_kernel(self,frame,trace,status,k0,k1,p0)
            if (status /= short_fft_success) return
            if (ib == 0) then
                dout = pack(dims,[1,2,3,4,5,6,7,8] /= t)
                dout(findloc(out_axes,f,dim=1)) = size(trace)
                deallocate(x)
                allocate(x(dout(1),dout(2),dout(3),dout(4),dout(5),dout(6),dout(7)))
            end if
            do j=1,6
                oidx(findloc(out_axes,batch_axes(j),dim=1)) = idx(batch_axes(j))
            end do
            do k=1,size(trace)
                oidx(findloc(out_axes,f,dim=1)) = k
                x(oidx(1),oidx(2),oidx(3),oidx(4),oidx(5),oidx(6),oidx(7)) = trace(k)
            end do
        end do
        status = short_fft_success
    end subroutine inverse_nd8_kernel

    subroutine spectrogram_nd7_real(self, x, result, p0, p1, padding, axis)
        class(short_time_fft), intent(in) :: self !! transform plan
        real(dp), intent(in) :: x(:,:,:,:,:,:,:) !! signal array
        type(short_time_spectrogram_result_8d), intent(out) :: result !! power spectrum
        integer, intent(in), optional :: p0 !! first time slice
        integer, intent(in), optional :: p1 !! exclusive final slice
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! selected sample dimension
        type(short_time_fft_result_8d) :: spectrum

        call self%stft(x,spectrum,p0=p0,p1=p1,padding=padding,axis=axis)
        result%status = spectrum%status
        if (spectrum%status /= short_fft_success) return
        result%power = abs(spectrum%values)**2
        result%frequencies = spectrum%frequencies
        result%times = spectrum%times
        result%p0 = spectrum%p0
    end subroutine spectrogram_nd7_real

    subroutine spectrogram_nd7_complex(self, x, result, p0, p1, padding, axis)
        class(short_time_fft), intent(in) :: self !! transform plan
        complex(dp), intent(in) :: x(:,:,:,:,:,:,:) !! signal array
        type(short_time_spectrogram_result_8d), intent(out) :: result !! power spectrum
        integer, intent(in), optional :: p0 !! first time slice
        integer, intent(in), optional :: p1 !! exclusive final slice
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! selected sample dimension
        type(short_time_fft_result_8d) :: spectrum

        call self%stft(x,spectrum,p0=p0,p1=p1,padding=padding,axis=axis)
        result%status = spectrum%status
        if (spectrum%status /= short_fft_success) return
        result%power = abs(spectrum%values)**2
        result%frequencies = spectrum%frequencies
        result%times = spectrum%times
        result%p0 = spectrum%p0
    end subroutine spectrogram_nd7_complex

    subroutine cross_spectrogram_nd7_real(self, x, y, values, status, p0, p1, padding, axis)
        class(short_time_fft), intent(in) :: self !! transform plan
        real(dp), intent(in) :: x(:,:,:,:,:,:,:) !! first array of signals
        real(dp), intent(in) :: y(:,:,:,:,:,:,:) !! second array of signals
        complex(dp), allocatable, intent(out) :: values(:,:,:,:,:,:,:,:) !! cross spectrum
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: p0 !! first time slice
        integer, intent(in), optional :: p1 !! exclusive final slice
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! selected sample dimension
        type(short_time_fft_result_8d) :: sx, sy

        allocate(values(0,0,0,0,0,0,0,0))
        status = short_fft_invalid_input
        if (any(shape(x) /= shape(y))) return
        call self%stft(x,sx,p0=p0,p1=p1,padding=padding,axis=axis)
        if (sx%status /= short_fft_success) return
        call self%stft(y,sy,p0=p0,p1=p1,padding=padding,axis=axis)
        if (sy%status /= short_fft_success) return
        values = sx%values*conjg(sy%values)
        status = short_fft_success
    end subroutine cross_spectrogram_nd7_real

    subroutine cross_spectrogram_nd7_complex(self, x, y, values, status, p0, p1, padding, axis)
        class(short_time_fft), intent(in) :: self !! transform plan
        complex(dp), intent(in) :: x(:,:,:,:,:,:,:) !! first array of signals
        complex(dp), intent(in) :: y(:,:,:,:,:,:,:) !! second array of signals
        complex(dp), allocatable, intent(out) :: values(:,:,:,:,:,:,:,:) !! cross spectrum
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: p0 !! first time slice
        integer, intent(in), optional :: p1 !! exclusive final slice
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! selected sample dimension
        type(short_time_fft_result_8d) :: sx, sy

        allocate(values(0,0,0,0,0,0,0,0))
        status = short_fft_invalid_input
        if (any(shape(x) /= shape(y))) return
        call self%stft(x,sx,p0=p0,p1=p1,padding=padding,axis=axis)
        if (sx%status /= short_fft_success) return
        call self%stft(y,sy,p0=p0,p1=p1,padding=padding,axis=axis)
        if (sy%status /= short_fft_success) return
        values = sx%values*conjg(sy%values)
        status = short_fft_success
    end subroutine cross_spectrogram_nd7_complex

    subroutine stft_nd8_real(self, x, result, p0, p1, k_offset, padding, axis)
        class(short_time_fft), intent(in) :: self !! configured STFT plan
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:) !! input time series
        type(short_time_fft_result_9d), intent(out) :: result !! spectrum
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: p1 !! exclusive upper time slice
        integer, intent(in), optional :: k_offset !! input sample offset
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! time-series axis, 1-based or negative
        call stft_nd8_kernel(self, cmplx(x,0.0_dp,dp), result, axis, p0, p1, &
            k_offset, padding, .false.)
    end subroutine stft_nd8_real

    subroutine stft_detrend_nd8_real(self, x, result, detrend, p0, p1, k_offset, padding, axis)
        class(short_time_fft), intent(in) :: self !! configured STFT plan
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:) !! input time series
        type(short_time_fft_result_9d), intent(out) :: result !! spectrum
        character(*), intent(in) :: detrend !! constant or linear segment detrending
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: p1 !! exclusive upper time slice
        integer, intent(in), optional :: k_offset !! input sample offset
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! time-series axis, 1-based or negative
        call stft_nd8_kernel(self, cmplx(x,0.0_dp,dp), result, axis, p0, p1, &
            k_offset, padding, .false., detrend)
    end subroutine stft_detrend_nd8_real

    subroutine stft_nd8_complex(self, x, result, p0, p1, k_offset, padding, axis)
        class(short_time_fft), intent(in) :: self !! configured STFT plan
        complex(dp), intent(in) :: x(:,:,:,:,:,:,:,:) !! input time series
        type(short_time_fft_result_9d), intent(out) :: result !! spectrum
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: p1 !! exclusive upper time slice
        integer, intent(in), optional :: k_offset !! input sample offset
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! time-series axis, 1-based or negative
        call stft_nd8_kernel(self, x, result, axis, p0, p1, &
            k_offset, padding, .true.)
    end subroutine stft_nd8_complex

    subroutine stft_detrend_nd8_complex(self, x, result, detrend, p0, p1, k_offset, padding, axis)
        class(short_time_fft), intent(in) :: self !! configured STFT plan
        complex(dp), intent(in) :: x(:,:,:,:,:,:,:,:) !! input time series
        type(short_time_fft_result_9d), intent(out) :: result !! spectrum
        character(*), intent(in) :: detrend !! constant or linear segment detrending
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: p1 !! exclusive upper time slice
        integer, intent(in), optional :: k_offset !! input sample offset
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! time-series axis, 1-based or negative
        call stft_nd8_kernel(self, x, result, axis, p0, p1, &
            k_offset, padding, .true., detrend)
    end subroutine stft_detrend_nd8_complex

    ! Flatten only the batch-index traversal: each one-dimensional signal is
    ! transformed by the same verified STFT kernel as ranks one through three.
    subroutine stft_nd8_kernel(self, x, result, axis, p0, p1, k_offset, padding, complex_input, detrend)
        class(short_time_fft), intent(in) :: self !! configured STFT plan
        complex(dp), intent(in) :: x(:,:,:,:,:,:,:,:) !! input signals
        type(short_time_fft_result_9d), intent(out) :: result !! output spectrum
        integer, intent(in), optional :: axis !! signal axis
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: p1 !! exclusive upper time slice
        integer, intent(in), optional :: k_offset !! sample offset
        character(*), intent(in), optional :: padding !! extension mode
        logical, intent(in) :: complex_input !! whether input is complex
        character(*), intent(in), optional :: detrend !! detrending mode
        type(short_time_fft_result) :: frame
        complex(dp), allocatable :: trace(:)
        integer :: a, dims(8), dout(9), idx(8), oidx(9)
        integer :: batch_axes(7), b, ib, j, k, q, t, nb

        result%status = short_fft_invalid_input
        a = fft_resolve_axis(axis,8,8)
        if (a == 0) return
        dims = shape(x)
        if (any(dims < 1)) return
        batch_axes = pack([1,2,3,4,5,6,7,8],[1,2,3,4,5,6,7,8] /= a)
        nb = product(dims(batch_axes))
        allocate(trace(dims(a)))
        do ib=0,nb-1
            b = ib
            do j=1,7
                idx(batch_axes(j)) = modulo(b,dims(batch_axes(j)))+1
                b = b/dims(batch_axes(j))
            end do
            do k=1,dims(a)
                idx(a) = k
                trace(k) = x(idx(1),idx(2),idx(3),idx(4),idx(5),idx(6),idx(7),idx(8))
            end do
            call stft_kernel(self,trace,frame,p0,p1,k_offset,padding,complex_input,detrend)
            if (frame%status /= short_fft_success) return
            if (ib == 0) then
                dout(1:8) = dims
                dout(a) = size(frame%values,1)
                dout(9) = size(frame%values,2)
                allocate(result%values(dout(1),dout(2),dout(3),dout(4),dout(5),dout(6),dout(7),dout(8),dout(9)))
                result%frequencies = frame%frequencies
                result%times = frame%times
                result%p0 = frame%p0
            end if
            do t=1,size(frame%values,2)
                do q=1,size(frame%values,1)
                    oidx(1:8) = idx
                    oidx(a) = q
                    oidx(9) = t
                    result%values(oidx(1),oidx(2),oidx(3),oidx(4),oidx(5),oidx(6),oidx(7),oidx(8),oidx(9)) = frame%values(q,t)
                end do
            end do
        end do
        result%status = short_fft_success
    end subroutine stft_nd8_kernel

    subroutine istft_nd9_real(self, s, x, status, k0, k1, p0, f_axis, t_axis)
        class(short_time_fft), intent(in) :: self !! inverse plan
        complex(dp), intent(in) :: s(:,:,:,:,:,:,:,:,:) !! STFT coefficients
        real(dp), allocatable, intent(out) :: x(:,:,:,:,:,:,:,:) !! reconstructed signal
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: k0 !! inclusive first sample
        integer, intent(in), optional :: k1 !! exclusive last sample
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: f_axis !! frequency axis of S
        integer, intent(in), optional :: t_axis !! time axis of S
        complex(dp), allocatable :: z(:,:,:,:,:,:,:,:)

        allocate(x(0,0,0,0,0,0,0,0))
        status = short_fft_invalid_input
        if (self%fft_mode /= 'onesided' .and. self%fft_mode /= 'onesided2X') return
        call inverse_nd9_kernel(self,s,z,status,k0,k1,p0,f_axis,t_axis)
        if (status == short_fft_success) x = real(z,dp)
    end subroutine istft_nd9_real

    subroutine istft_nd9_complex(self, s, x, status, k0, k1, p0, f_axis, t_axis)
        class(short_time_fft), intent(in) :: self !! inverse plan
        complex(dp), intent(in) :: s(:,:,:,:,:,:,:,:,:) !! STFT coefficients
        complex(dp), allocatable, intent(out) :: x(:,:,:,:,:,:,:,:) !! reconstructed signal
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: k0 !! inclusive first sample
        integer, intent(in), optional :: k1 !! exclusive last sample
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: f_axis !! frequency axis of S
        integer, intent(in), optional :: t_axis !! time axis of S
        complex(dp), allocatable :: z(:,:,:,:,:,:,:,:)

        allocate(x(0,0,0,0,0,0,0,0))
        status = short_fft_invalid_input
        if (self%fft_mode == 'onesided' .or. self%fft_mode == 'onesided2X') return
        call inverse_nd9_kernel(self,s,z,status,k0,k1,p0,f_axis,t_axis)
        if (status == short_fft_success) x = z
    end subroutine istft_nd9_complex

    subroutine inverse_nd9_kernel(self, s, x, status, k0, k1, p0, f_axis, t_axis)
        class(short_time_fft), intent(in) :: self !! inverse transform plan
        complex(dp), intent(in) :: s(:,:,:,:,:,:,:,:,:) !! STFT coefficients
        complex(dp), allocatable, intent(out) :: x(:,:,:,:,:,:,:,:) !! reconstructed series
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: k0 !! inclusive first sample
        integer, intent(in), optional :: k1 !! exclusive last sample
        integer, intent(in), optional :: p0 !! first time slice
        integer, intent(in), optional :: f_axis !! frequency dimension
        integer, intent(in), optional :: t_axis !! time dimension
        complex(dp), allocatable :: frame(:,:), trace(:)
        integer :: f, t, dims(9), dout(8), batch_axes(7), out_axes(8)
        integer :: idx(9), oidx(8), b, ib, j, q, p, k, nb

        allocate(x(0,0,0,0,0,0,0,0))
        status = short_fft_invalid_input
        f = fft_resolve_axis(f_axis,9,8)
        t = fft_resolve_axis(t_axis,9,9)
        if (f == 0 .or. t == 0 .or. f == t) return
        dims = shape(s)
        if (any(dims < 1)) return
        batch_axes = pack([1,2,3,4,5,6,7,8,9],[1,2,3,4,5,6,7,8,9] /= f .and. [1,2,3,4,5,6,7,8,9] /= t)
        out_axes = pack([1,2,3,4,5,6,7,8,9],[1,2,3,4,5,6,7,8,9] /= t)
        nb = product(dims(batch_axes))
        allocate(frame(dims(f),dims(t)))
        do ib=0,nb-1
            b = ib
            do j=1,7
                idx(batch_axes(j)) = modulo(b,dims(batch_axes(j)))+1
                b = b/dims(batch_axes(j))
            end do
            do p=1,dims(t)
                idx(t) = p
                do q=1,dims(f)
                    idx(f) = q
                    frame(q,p) = s(idx(1),idx(2),idx(3),idx(4),idx(5),idx(6),idx(7),idx(8),idx(9))
                end do
            end do
            call inverse_kernel(self,frame,trace,status,k0,k1,p0)
            if (status /= short_fft_success) return
            if (ib == 0) then
                dout = pack(dims,[1,2,3,4,5,6,7,8,9] /= t)
                dout(findloc(out_axes,f,dim=1)) = size(trace)
                deallocate(x)
                allocate(x(dout(1),dout(2),dout(3),dout(4),dout(5),dout(6),dout(7),dout(8)))
            end if
            do j=1,7
                oidx(findloc(out_axes,batch_axes(j),dim=1)) = idx(batch_axes(j))
            end do
            do k=1,size(trace)
                oidx(findloc(out_axes,f,dim=1)) = k
                x(oidx(1),oidx(2),oidx(3),oidx(4),oidx(5),oidx(6),oidx(7),oidx(8)) = trace(k)
            end do
        end do
        status = short_fft_success
    end subroutine inverse_nd9_kernel

    subroutine spectrogram_nd8_real(self, x, result, p0, p1, padding, axis)
        class(short_time_fft), intent(in) :: self !! transform plan
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:) !! signal array
        type(short_time_spectrogram_result_9d), intent(out) :: result !! power spectrum
        integer, intent(in), optional :: p0 !! first time slice
        integer, intent(in), optional :: p1 !! exclusive final slice
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! selected sample dimension
        type(short_time_fft_result_9d) :: spectrum

        call self%stft(x,spectrum,p0=p0,p1=p1,padding=padding,axis=axis)
        result%status = spectrum%status
        if (spectrum%status /= short_fft_success) return
        result%power = abs(spectrum%values)**2
        result%frequencies = spectrum%frequencies
        result%times = spectrum%times
        result%p0 = spectrum%p0
    end subroutine spectrogram_nd8_real

    subroutine spectrogram_nd8_complex(self, x, result, p0, p1, padding, axis)
        class(short_time_fft), intent(in) :: self !! transform plan
        complex(dp), intent(in) :: x(:,:,:,:,:,:,:,:) !! signal array
        type(short_time_spectrogram_result_9d), intent(out) :: result !! power spectrum
        integer, intent(in), optional :: p0 !! first time slice
        integer, intent(in), optional :: p1 !! exclusive final slice
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! selected sample dimension
        type(short_time_fft_result_9d) :: spectrum

        call self%stft(x,spectrum,p0=p0,p1=p1,padding=padding,axis=axis)
        result%status = spectrum%status
        if (spectrum%status /= short_fft_success) return
        result%power = abs(spectrum%values)**2
        result%frequencies = spectrum%frequencies
        result%times = spectrum%times
        result%p0 = spectrum%p0
    end subroutine spectrogram_nd8_complex

    subroutine cross_spectrogram_nd8_real(self, x, y, values, status, p0, p1, padding, axis)
        class(short_time_fft), intent(in) :: self !! transform plan
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:) !! first array of signals
        real(dp), intent(in) :: y(:,:,:,:,:,:,:,:) !! second array of signals
        complex(dp), allocatable, intent(out) :: values(:,:,:,:,:,:,:,:,:) !! cross spectrum
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: p0 !! first time slice
        integer, intent(in), optional :: p1 !! exclusive final slice
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! selected sample dimension
        type(short_time_fft_result_9d) :: sx, sy

        allocate(values(0,0,0,0,0,0,0,0,0))
        status = short_fft_invalid_input
        if (any(shape(x) /= shape(y))) return
        call self%stft(x,sx,p0=p0,p1=p1,padding=padding,axis=axis)
        if (sx%status /= short_fft_success) return
        call self%stft(y,sy,p0=p0,p1=p1,padding=padding,axis=axis)
        if (sy%status /= short_fft_success) return
        values = sx%values*conjg(sy%values)
        status = short_fft_success
    end subroutine cross_spectrogram_nd8_real

    subroutine cross_spectrogram_nd8_complex(self, x, y, values, status, p0, p1, padding, axis)
        class(short_time_fft), intent(in) :: self !! transform plan
        complex(dp), intent(in) :: x(:,:,:,:,:,:,:,:) !! first array of signals
        complex(dp), intent(in) :: y(:,:,:,:,:,:,:,:) !! second array of signals
        complex(dp), allocatable, intent(out) :: values(:,:,:,:,:,:,:,:,:) !! cross spectrum
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: p0 !! first time slice
        integer, intent(in), optional :: p1 !! exclusive final slice
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! selected sample dimension
        type(short_time_fft_result_9d) :: sx, sy

        allocate(values(0,0,0,0,0,0,0,0,0))
        status = short_fft_invalid_input
        if (any(shape(x) /= shape(y))) return
        call self%stft(x,sx,p0=p0,p1=p1,padding=padding,axis=axis)
        if (sx%status /= short_fft_success) return
        call self%stft(y,sy,p0=p0,p1=p1,padding=padding,axis=axis)
        if (sy%status /= short_fft_success) return
        values = sx%values*conjg(sy%values)
        status = short_fft_success
    end subroutine cross_spectrogram_nd8_complex

    subroutine stft_nd9_real(self, x, result, p0, p1, k_offset, padding, axis)
        class(short_time_fft), intent(in) :: self !! configured STFT plan
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:) !! input time series
        type(short_time_fft_result_10d), intent(out) :: result !! spectrum
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: p1 !! exclusive upper time slice
        integer, intent(in), optional :: k_offset !! input sample offset
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! time-series axis, 1-based or negative
        call stft_nd9_kernel(self, cmplx(x,0.0_dp,dp), result, axis, p0, p1, &
            k_offset, padding, .false.)
    end subroutine stft_nd9_real

    subroutine stft_detrend_nd9_real(self, x, result, detrend, p0, p1, k_offset, padding, axis)
        class(short_time_fft), intent(in) :: self !! configured STFT plan
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:) !! input time series
        type(short_time_fft_result_10d), intent(out) :: result !! spectrum
        character(*), intent(in) :: detrend !! constant or linear segment detrending
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: p1 !! exclusive upper time slice
        integer, intent(in), optional :: k_offset !! input sample offset
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! time-series axis, 1-based or negative
        call stft_nd9_kernel(self, cmplx(x,0.0_dp,dp), result, axis, p0, p1, &
            k_offset, padding, .false., detrend)
    end subroutine stft_detrend_nd9_real

    subroutine stft_nd9_complex(self, x, result, p0, p1, k_offset, padding, axis)
        class(short_time_fft), intent(in) :: self !! configured STFT plan
        complex(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:) !! input time series
        type(short_time_fft_result_10d), intent(out) :: result !! spectrum
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: p1 !! exclusive upper time slice
        integer, intent(in), optional :: k_offset !! input sample offset
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! time-series axis, 1-based or negative
        call stft_nd9_kernel(self, x, result, axis, p0, p1, &
            k_offset, padding, .true.)
    end subroutine stft_nd9_complex

    subroutine stft_detrend_nd9_complex(self, x, result, detrend, p0, p1, k_offset, padding, axis)
        class(short_time_fft), intent(in) :: self !! configured STFT plan
        complex(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:) !! input time series
        type(short_time_fft_result_10d), intent(out) :: result !! spectrum
        character(*), intent(in) :: detrend !! constant or linear segment detrending
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: p1 !! exclusive upper time slice
        integer, intent(in), optional :: k_offset !! input sample offset
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! time-series axis, 1-based or negative
        call stft_nd9_kernel(self, x, result, axis, p0, p1, &
            k_offset, padding, .true., detrend)
    end subroutine stft_detrend_nd9_complex

    ! Flatten only the batch-index traversal: each one-dimensional signal is
    ! transformed by the same verified STFT kernel as ranks one through three.
    subroutine stft_nd9_kernel(self, x, result, axis, p0, p1, k_offset, padding, complex_input, detrend)
        class(short_time_fft), intent(in) :: self !! configured STFT plan
        complex(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:) !! input signals
        type(short_time_fft_result_10d), intent(out) :: result !! output spectrum
        integer, intent(in), optional :: axis !! signal axis
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: p1 !! exclusive upper time slice
        integer, intent(in), optional :: k_offset !! sample offset
        character(*), intent(in), optional :: padding !! extension mode
        logical, intent(in) :: complex_input !! whether input is complex
        character(*), intent(in), optional :: detrend !! detrending mode
        type(short_time_fft_result) :: frame
        complex(dp), allocatable :: trace(:)
        integer :: a, dims(9), dout(10), idx(9), oidx(10)
        integer :: batch_axes(8), b, ib, j, k, q, t, nb

        result%status = short_fft_invalid_input
        a = fft_resolve_axis(axis,9,9)
        if (a == 0) return
        dims = shape(x)
        if (any(dims < 1)) return
        batch_axes = pack([1,2,3,4,5,6,7,8,9],[1,2,3,4,5,6,7,8,9] /= a)
        nb = product(dims(batch_axes))
        allocate(trace(dims(a)))
        do ib=0,nb-1
            b = ib
            do j=1,8
                idx(batch_axes(j)) = modulo(b,dims(batch_axes(j)))+1
                b = b/dims(batch_axes(j))
            end do
            do k=1,dims(a)
                idx(a) = k
                trace(k) = x(idx(1),idx(2),idx(3),idx(4),idx(5),idx(6),idx(7),idx(8),idx(9))
            end do
            call stft_kernel(self,trace,frame,p0,p1,k_offset,padding,complex_input,detrend)
            if (frame%status /= short_fft_success) return
            if (ib == 0) then
                dout(1:9) = dims
                dout(a) = size(frame%values,1)
                dout(10) = size(frame%values,2)
                allocate(result%values(dout(1),dout(2),dout(3),dout(4),dout(5),dout(6),dout(7),dout(8),dout(9),dout(10)))
                result%frequencies = frame%frequencies
                result%times = frame%times
                result%p0 = frame%p0
            end if
            do t=1,size(frame%values,2)
                do q=1,size(frame%values,1)
                    oidx(1:9) = idx
                    oidx(a) = q
                    oidx(10) = t
                    result%values(oidx(1),oidx(2),oidx(3),oidx(4),oidx(5),oidx(6),oidx(7),oidx(8),oidx(9), &
            oidx(10)) = frame%values(q,t)
                end do
            end do
        end do
        result%status = short_fft_success
    end subroutine stft_nd9_kernel

    subroutine istft_nd10_real(self, s, x, status, k0, k1, p0, f_axis, t_axis)
        class(short_time_fft), intent(in) :: self !! inverse plan
        complex(dp), intent(in) :: s(:,:,:,:,:,:,:,:,:,:) !! STFT coefficients
        real(dp), allocatable, intent(out) :: x(:,:,:,:,:,:,:,:,:) !! reconstructed signal
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: k0 !! inclusive first sample
        integer, intent(in), optional :: k1 !! exclusive last sample
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: f_axis !! frequency axis of S
        integer, intent(in), optional :: t_axis !! time axis of S
        complex(dp), allocatable :: z(:,:,:,:,:,:,:,:,:)

        allocate(x(0,0,0,0,0,0,0,0,0))
        status = short_fft_invalid_input
        if (self%fft_mode /= 'onesided' .and. self%fft_mode /= 'onesided2X') return
        call inverse_nd10_kernel(self,s,z,status,k0,k1,p0,f_axis,t_axis)
        if (status == short_fft_success) x = real(z,dp)
    end subroutine istft_nd10_real

    subroutine istft_nd10_complex(self, s, x, status, k0, k1, p0, f_axis, t_axis)
        class(short_time_fft), intent(in) :: self !! inverse plan
        complex(dp), intent(in) :: s(:,:,:,:,:,:,:,:,:,:) !! STFT coefficients
        complex(dp), allocatable, intent(out) :: x(:,:,:,:,:,:,:,:,:) !! reconstructed signal
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: k0 !! inclusive first sample
        integer, intent(in), optional :: k1 !! exclusive last sample
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: f_axis !! frequency axis of S
        integer, intent(in), optional :: t_axis !! time axis of S
        complex(dp), allocatable :: z(:,:,:,:,:,:,:,:,:)

        allocate(x(0,0,0,0,0,0,0,0,0))
        status = short_fft_invalid_input
        if (self%fft_mode == 'onesided' .or. self%fft_mode == 'onesided2X') return
        call inverse_nd10_kernel(self,s,z,status,k0,k1,p0,f_axis,t_axis)
        if (status == short_fft_success) x = z
    end subroutine istft_nd10_complex

    subroutine inverse_nd10_kernel(self, s, x, status, k0, k1, p0, f_axis, t_axis)
        class(short_time_fft), intent(in) :: self !! inverse transform plan
        complex(dp), intent(in) :: s(:,:,:,:,:,:,:,:,:,:) !! STFT coefficients
        complex(dp), allocatable, intent(out) :: x(:,:,:,:,:,:,:,:,:) !! reconstructed series
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: k0 !! inclusive first sample
        integer, intent(in), optional :: k1 !! exclusive last sample
        integer, intent(in), optional :: p0 !! first time slice
        integer, intent(in), optional :: f_axis !! frequency dimension
        integer, intent(in), optional :: t_axis !! time dimension
        complex(dp), allocatable :: frame(:,:), trace(:)
        integer :: f, t, dims(10), dout(9), batch_axes(8), out_axes(9)
        integer :: idx(10), oidx(9), b, ib, j, q, p, k, nb

        allocate(x(0,0,0,0,0,0,0,0,0))
        status = short_fft_invalid_input
        f = fft_resolve_axis(f_axis,10,9)
        t = fft_resolve_axis(t_axis,10,10)
        if (f == 0 .or. t == 0 .or. f == t) return
        dims = shape(s)
        if (any(dims < 1)) return
        batch_axes = pack([1,2,3,4,5,6,7,8,9,10],[1,2,3,4,5,6,7,8,9,10] /= f .and. [1,2,3,4,5,6,7,8,9,10] /= t)
        out_axes = pack([1,2,3,4,5,6,7,8,9,10],[1,2,3,4,5,6,7,8,9,10] /= t)
        nb = product(dims(batch_axes))
        allocate(frame(dims(f),dims(t)))
        do ib=0,nb-1
            b = ib
            do j=1,8
                idx(batch_axes(j)) = modulo(b,dims(batch_axes(j)))+1
                b = b/dims(batch_axes(j))
            end do
            do p=1,dims(t)
                idx(t) = p
                do q=1,dims(f)
                    idx(f) = q
                    frame(q,p) = s(idx(1),idx(2),idx(3),idx(4),idx(5),idx(6),idx(7),idx(8),idx(9),idx(10))
                end do
            end do
            call inverse_kernel(self,frame,trace,status,k0,k1,p0)
            if (status /= short_fft_success) return
            if (ib == 0) then
                dout = pack(dims,[1,2,3,4,5,6,7,8,9,10] /= t)
                dout(findloc(out_axes,f,dim=1)) = size(trace)
                deallocate(x)
                allocate(x(dout(1),dout(2),dout(3),dout(4),dout(5),dout(6),dout(7),dout(8),dout(9)))
            end if
            do j=1,8
                oidx(findloc(out_axes,batch_axes(j),dim=1)) = idx(batch_axes(j))
            end do
            do k=1,size(trace)
                oidx(findloc(out_axes,f,dim=1)) = k
                x(oidx(1),oidx(2),oidx(3),oidx(4),oidx(5),oidx(6),oidx(7),oidx(8),oidx(9)) = trace(k)
            end do
        end do
        status = short_fft_success
    end subroutine inverse_nd10_kernel

    subroutine spectrogram_nd9_real(self, x, result, p0, p1, padding, axis)
        class(short_time_fft), intent(in) :: self !! transform plan
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:) !! signal array
        type(short_time_spectrogram_result_10d), intent(out) :: result !! power spectrum
        integer, intent(in), optional :: p0 !! first time slice
        integer, intent(in), optional :: p1 !! exclusive final slice
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! selected sample dimension
        type(short_time_fft_result_10d) :: spectrum

        call self%stft(x,spectrum,p0=p0,p1=p1,padding=padding,axis=axis)
        result%status = spectrum%status
        if (spectrum%status /= short_fft_success) return
        result%power = abs(spectrum%values)**2
        result%frequencies = spectrum%frequencies
        result%times = spectrum%times
        result%p0 = spectrum%p0
    end subroutine spectrogram_nd9_real

    subroutine spectrogram_nd9_complex(self, x, result, p0, p1, padding, axis)
        class(short_time_fft), intent(in) :: self !! transform plan
        complex(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:) !! signal array
        type(short_time_spectrogram_result_10d), intent(out) :: result !! power spectrum
        integer, intent(in), optional :: p0 !! first time slice
        integer, intent(in), optional :: p1 !! exclusive final slice
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! selected sample dimension
        type(short_time_fft_result_10d) :: spectrum

        call self%stft(x,spectrum,p0=p0,p1=p1,padding=padding,axis=axis)
        result%status = spectrum%status
        if (spectrum%status /= short_fft_success) return
        result%power = abs(spectrum%values)**2
        result%frequencies = spectrum%frequencies
        result%times = spectrum%times
        result%p0 = spectrum%p0
    end subroutine spectrogram_nd9_complex

    subroutine cross_spectrogram_nd9_real(self, x, y, values, status, p0, p1, padding, axis)
        class(short_time_fft), intent(in) :: self !! transform plan
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:) !! first array of signals
        real(dp), intent(in) :: y(:,:,:,:,:,:,:,:,:) !! second array of signals
        complex(dp), allocatable, intent(out) :: values(:,:,:,:,:,:,:,:,:,:) !! cross spectrum
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: p0 !! first time slice
        integer, intent(in), optional :: p1 !! exclusive final slice
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! selected sample dimension
        type(short_time_fft_result_10d) :: sx, sy

        allocate(values(0,0,0,0,0,0,0,0,0,0))
        status = short_fft_invalid_input
        if (any(shape(x) /= shape(y))) return
        call self%stft(x,sx,p0=p0,p1=p1,padding=padding,axis=axis)
        if (sx%status /= short_fft_success) return
        call self%stft(y,sy,p0=p0,p1=p1,padding=padding,axis=axis)
        if (sy%status /= short_fft_success) return
        values = sx%values*conjg(sy%values)
        status = short_fft_success
    end subroutine cross_spectrogram_nd9_real

    subroutine cross_spectrogram_nd9_complex(self, x, y, values, status, p0, p1, padding, axis)
        class(short_time_fft), intent(in) :: self !! transform plan
        complex(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:) !! first array of signals
        complex(dp), intent(in) :: y(:,:,:,:,:,:,:,:,:) !! second array of signals
        complex(dp), allocatable, intent(out) :: values(:,:,:,:,:,:,:,:,:,:) !! cross spectrum
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: p0 !! first time slice
        integer, intent(in), optional :: p1 !! exclusive final slice
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! selected sample dimension
        type(short_time_fft_result_10d) :: sx, sy

        allocate(values(0,0,0,0,0,0,0,0,0,0))
        status = short_fft_invalid_input
        if (any(shape(x) /= shape(y))) return
        call self%stft(x,sx,p0=p0,p1=p1,padding=padding,axis=axis)
        if (sx%status /= short_fft_success) return
        call self%stft(y,sy,p0=p0,p1=p1,padding=padding,axis=axis)
        if (sy%status /= short_fft_success) return
        values = sx%values*conjg(sy%values)
        status = short_fft_success
    end subroutine cross_spectrogram_nd9_complex

    subroutine stft_nd10_real(self, x, result, p0, p1, k_offset, padding, axis)
        class(short_time_fft), intent(in) :: self !! configured STFT plan
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:) !! input time series
        type(short_time_fft_result_11d), intent(out) :: result !! spectrum
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: p1 !! exclusive upper time slice
        integer, intent(in), optional :: k_offset !! input sample offset
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! time-series axis, 1-based or negative
        call stft_nd10_kernel(self, cmplx(x,0.0_dp,dp), result, axis, p0, p1, &
            k_offset, padding, .false.)
    end subroutine stft_nd10_real

    subroutine stft_detrend_nd10_real(self, x, result, detrend, p0, p1, k_offset, padding, axis)
        class(short_time_fft), intent(in) :: self !! configured STFT plan
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:) !! input time series
        type(short_time_fft_result_11d), intent(out) :: result !! spectrum
        character(*), intent(in) :: detrend !! constant or linear segment detrending
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: p1 !! exclusive upper time slice
        integer, intent(in), optional :: k_offset !! input sample offset
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! time-series axis, 1-based or negative
        call stft_nd10_kernel(self, cmplx(x,0.0_dp,dp), result, axis, p0, p1, &
            k_offset, padding, .false., detrend)
    end subroutine stft_detrend_nd10_real

    subroutine stft_nd10_complex(self, x, result, p0, p1, k_offset, padding, axis)
        class(short_time_fft), intent(in) :: self !! configured STFT plan
        complex(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:) !! input time series
        type(short_time_fft_result_11d), intent(out) :: result !! spectrum
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: p1 !! exclusive upper time slice
        integer, intent(in), optional :: k_offset !! input sample offset
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! time-series axis, 1-based or negative
        call stft_nd10_kernel(self, x, result, axis, p0, p1, &
            k_offset, padding, .true.)
    end subroutine stft_nd10_complex

    subroutine stft_detrend_nd10_complex(self, x, result, detrend, p0, p1, k_offset, padding, axis)
        class(short_time_fft), intent(in) :: self !! configured STFT plan
        complex(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:) !! input time series
        type(short_time_fft_result_11d), intent(out) :: result !! spectrum
        character(*), intent(in) :: detrend !! constant or linear segment detrending
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: p1 !! exclusive upper time slice
        integer, intent(in), optional :: k_offset !! input sample offset
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! time-series axis, 1-based or negative
        call stft_nd10_kernel(self, x, result, axis, p0, p1, &
            k_offset, padding, .true., detrend)
    end subroutine stft_detrend_nd10_complex

    ! Flatten only the batch-index traversal: each one-dimensional signal is
    ! transformed by the same verified STFT kernel as ranks one through three.
    subroutine stft_nd10_kernel(self, x, result, axis, p0, p1, k_offset, padding, complex_input, detrend)
        class(short_time_fft), intent(in) :: self !! configured STFT plan
        complex(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:) !! input signals
        type(short_time_fft_result_11d), intent(out) :: result !! output spectrum
        integer, intent(in), optional :: axis !! signal axis
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: p1 !! exclusive upper time slice
        integer, intent(in), optional :: k_offset !! sample offset
        character(*), intent(in), optional :: padding !! extension mode
        logical, intent(in) :: complex_input !! whether input is complex
        character(*), intent(in), optional :: detrend !! detrending mode
        type(short_time_fft_result) :: frame
        complex(dp), allocatable :: trace(:)
        integer :: a, dims(10), dout(11), idx(10), oidx(11)
        integer :: batch_axes(9), b, ib, j, k, q, t, nb

        result%status = short_fft_invalid_input
        a = fft_resolve_axis(axis,10,10)
        if (a == 0) return
        dims = shape(x)
        if (any(dims < 1)) return
        batch_axes = pack([1,2,3,4,5,6,7,8,9,10],[1,2,3,4,5,6,7,8,9,10] /= a)
        nb = product(dims(batch_axes))
        allocate(trace(dims(a)))
        do ib=0,nb-1
            b = ib
            do j=1,9
                idx(batch_axes(j)) = modulo(b,dims(batch_axes(j)))+1
                b = b/dims(batch_axes(j))
            end do
            do k=1,dims(a)
                idx(a) = k
                trace(k) = x(idx(1),idx(2),idx(3),idx(4),idx(5),idx(6),idx(7),idx(8),idx(9),idx(10))
            end do
            call stft_kernel(self,trace,frame,p0,p1,k_offset,padding,complex_input,detrend)
            if (frame%status /= short_fft_success) return
            if (ib == 0) then
                dout(1:10) = dims
                dout(a) = size(frame%values,1)
                dout(11) = size(frame%values,2)
                allocate(result%values(dout(1),dout(2),dout(3),dout(4),dout(5),dout(6),dout(7),dout(8),dout(9),dout(10),dout(11)))
                result%frequencies = frame%frequencies
                result%times = frame%times
                result%p0 = frame%p0
            end if
            do t=1,size(frame%values,2)
                do q=1,size(frame%values,1)
                    oidx(1:10) = idx
                    oidx(a) = q
                    oidx(11) = t
                    result%values(oidx(1),oidx(2),oidx(3),oidx(4),oidx(5),oidx(6),oidx(7),oidx(8),oidx(9), &
            oidx(10),oidx(11)) = frame%values(q,t)
                end do
            end do
        end do
        result%status = short_fft_success
    end subroutine stft_nd10_kernel

    subroutine istft_nd11_real(self, s, x, status, k0, k1, p0, f_axis, t_axis)
        class(short_time_fft), intent(in) :: self !! inverse plan
        complex(dp), intent(in) :: s(:,:,:,:,:,:,:,:,:,:,:) !! STFT coefficients
        real(dp), allocatable, intent(out) :: x(:,:,:,:,:,:,:,:,:,:) !! reconstructed signal
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: k0 !! inclusive first sample
        integer, intent(in), optional :: k1 !! exclusive last sample
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: f_axis !! frequency axis of S
        integer, intent(in), optional :: t_axis !! time axis of S
        complex(dp), allocatable :: z(:,:,:,:,:,:,:,:,:,:)

        allocate(x(0,0,0,0,0,0,0,0,0,0))
        status = short_fft_invalid_input
        if (self%fft_mode /= 'onesided' .and. self%fft_mode /= 'onesided2X') return
        call inverse_nd11_kernel(self,s,z,status,k0,k1,p0,f_axis,t_axis)
        if (status == short_fft_success) x = real(z,dp)
    end subroutine istft_nd11_real

    subroutine istft_nd11_complex(self, s, x, status, k0, k1, p0, f_axis, t_axis)
        class(short_time_fft), intent(in) :: self !! inverse plan
        complex(dp), intent(in) :: s(:,:,:,:,:,:,:,:,:,:,:) !! STFT coefficients
        complex(dp), allocatable, intent(out) :: x(:,:,:,:,:,:,:,:,:,:) !! reconstructed signal
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: k0 !! inclusive first sample
        integer, intent(in), optional :: k1 !! exclusive last sample
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: f_axis !! frequency axis of S
        integer, intent(in), optional :: t_axis !! time axis of S
        complex(dp), allocatable :: z(:,:,:,:,:,:,:,:,:,:)

        allocate(x(0,0,0,0,0,0,0,0,0,0))
        status = short_fft_invalid_input
        if (self%fft_mode == 'onesided' .or. self%fft_mode == 'onesided2X') return
        call inverse_nd11_kernel(self,s,z,status,k0,k1,p0,f_axis,t_axis)
        if (status == short_fft_success) x = z
    end subroutine istft_nd11_complex

    subroutine inverse_nd11_kernel(self, s, x, status, k0, k1, p0, f_axis, t_axis)
        class(short_time_fft), intent(in) :: self !! inverse transform plan
        complex(dp), intent(in) :: s(:,:,:,:,:,:,:,:,:,:,:) !! STFT coefficients
        complex(dp), allocatable, intent(out) :: x(:,:,:,:,:,:,:,:,:,:) !! reconstructed series
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: k0 !! inclusive first sample
        integer, intent(in), optional :: k1 !! exclusive last sample
        integer, intent(in), optional :: p0 !! first time slice
        integer, intent(in), optional :: f_axis !! frequency dimension
        integer, intent(in), optional :: t_axis !! time dimension
        complex(dp), allocatable :: frame(:,:), trace(:)
        integer :: f, t, dims(11), dout(10), batch_axes(9), out_axes(10)
        integer :: idx(11), oidx(10), b, ib, j, q, p, k, nb

        allocate(x(0,0,0,0,0,0,0,0,0,0))
        status = short_fft_invalid_input
        f = fft_resolve_axis(f_axis,11,10)
        t = fft_resolve_axis(t_axis,11,11)
        if (f == 0 .or. t == 0 .or. f == t) return
        dims = shape(s)
        if (any(dims < 1)) return
        batch_axes = pack([1,2,3,4,5,6,7,8,9,10,11],[1,2,3,4,5,6,7,8,9,10,11] /= f .and. [1,2,3,4,5,6,7,8,9,10,11] /= t)
        out_axes = pack([1,2,3,4,5,6,7,8,9,10,11],[1,2,3,4,5,6,7,8,9,10,11] /= t)
        nb = product(dims(batch_axes))
        allocate(frame(dims(f),dims(t)))
        do ib=0,nb-1
            b = ib
            do j=1,9
                idx(batch_axes(j)) = modulo(b,dims(batch_axes(j)))+1
                b = b/dims(batch_axes(j))
            end do
            do p=1,dims(t)
                idx(t) = p
                do q=1,dims(f)
                    idx(f) = q
                    frame(q,p) = s(idx(1),idx(2),idx(3),idx(4),idx(5),idx(6),idx(7),idx(8),idx(9),idx(10),idx(11))
                end do
            end do
            call inverse_kernel(self,frame,trace,status,k0,k1,p0)
            if (status /= short_fft_success) return
            if (ib == 0) then
                dout = pack(dims,[1,2,3,4,5,6,7,8,9,10,11] /= t)
                dout(findloc(out_axes,f,dim=1)) = size(trace)
                deallocate(x)
                allocate(x(dout(1),dout(2),dout(3),dout(4),dout(5),dout(6),dout(7),dout(8),dout(9),dout(10)))
            end if
            do j=1,9
                oidx(findloc(out_axes,batch_axes(j),dim=1)) = idx(batch_axes(j))
            end do
            do k=1,size(trace)
                oidx(findloc(out_axes,f,dim=1)) = k
                x(oidx(1),oidx(2),oidx(3),oidx(4),oidx(5),oidx(6),oidx(7),oidx(8),oidx(9),oidx(10)) = trace(k)
            end do
        end do
        status = short_fft_success
    end subroutine inverse_nd11_kernel

    subroutine spectrogram_nd10_real(self, x, result, p0, p1, padding, axis)
        class(short_time_fft), intent(in) :: self !! transform plan
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:) !! signal array
        type(short_time_spectrogram_result_11d), intent(out) :: result !! power spectrum
        integer, intent(in), optional :: p0 !! first time slice
        integer, intent(in), optional :: p1 !! exclusive final slice
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! selected sample dimension
        type(short_time_fft_result_11d) :: spectrum

        call self%stft(x,spectrum,p0=p0,p1=p1,padding=padding,axis=axis)
        result%status = spectrum%status
        if (spectrum%status /= short_fft_success) return
        result%power = abs(spectrum%values)**2
        result%frequencies = spectrum%frequencies
        result%times = spectrum%times
        result%p0 = spectrum%p0
    end subroutine spectrogram_nd10_real

    subroutine spectrogram_nd10_complex(self, x, result, p0, p1, padding, axis)
        class(short_time_fft), intent(in) :: self !! transform plan
        complex(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:) !! signal array
        type(short_time_spectrogram_result_11d), intent(out) :: result !! power spectrum
        integer, intent(in), optional :: p0 !! first time slice
        integer, intent(in), optional :: p1 !! exclusive final slice
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! selected sample dimension
        type(short_time_fft_result_11d) :: spectrum

        call self%stft(x,spectrum,p0=p0,p1=p1,padding=padding,axis=axis)
        result%status = spectrum%status
        if (spectrum%status /= short_fft_success) return
        result%power = abs(spectrum%values)**2
        result%frequencies = spectrum%frequencies
        result%times = spectrum%times
        result%p0 = spectrum%p0
    end subroutine spectrogram_nd10_complex

    subroutine cross_spectrogram_nd10_real(self, x, y, values, status, p0, p1, padding, axis)
        class(short_time_fft), intent(in) :: self !! transform plan
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:) !! first array of signals
        real(dp), intent(in) :: y(:,:,:,:,:,:,:,:,:,:) !! second array of signals
        complex(dp), allocatable, intent(out) :: values(:,:,:,:,:,:,:,:,:,:,:) !! cross spectrum
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: p0 !! first time slice
        integer, intent(in), optional :: p1 !! exclusive final slice
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! selected sample dimension
        type(short_time_fft_result_11d) :: sx, sy

        allocate(values(0,0,0,0,0,0,0,0,0,0,0))
        status = short_fft_invalid_input
        if (any(shape(x) /= shape(y))) return
        call self%stft(x,sx,p0=p0,p1=p1,padding=padding,axis=axis)
        if (sx%status /= short_fft_success) return
        call self%stft(y,sy,p0=p0,p1=p1,padding=padding,axis=axis)
        if (sy%status /= short_fft_success) return
        values = sx%values*conjg(sy%values)
        status = short_fft_success
    end subroutine cross_spectrogram_nd10_real

    subroutine cross_spectrogram_nd10_complex(self, x, y, values, status, p0, p1, padding, axis)
        class(short_time_fft), intent(in) :: self !! transform plan
        complex(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:) !! first array of signals
        complex(dp), intent(in) :: y(:,:,:,:,:,:,:,:,:,:) !! second array of signals
        complex(dp), allocatable, intent(out) :: values(:,:,:,:,:,:,:,:,:,:,:) !! cross spectrum
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: p0 !! first time slice
        integer, intent(in), optional :: p1 !! exclusive final slice
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! selected sample dimension
        type(short_time_fft_result_11d) :: sx, sy

        allocate(values(0,0,0,0,0,0,0,0,0,0,0))
        status = short_fft_invalid_input
        if (any(shape(x) /= shape(y))) return
        call self%stft(x,sx,p0=p0,p1=p1,padding=padding,axis=axis)
        if (sx%status /= short_fft_success) return
        call self%stft(y,sy,p0=p0,p1=p1,padding=padding,axis=axis)
        if (sy%status /= short_fft_success) return
        values = sx%values*conjg(sy%values)
        status = short_fft_success
    end subroutine cross_spectrogram_nd10_complex

    subroutine stft_nd11_real(self, x, result, p0, p1, k_offset, padding, axis)
        class(short_time_fft), intent(in) :: self !! configured STFT plan
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:) !! input time series
        type(short_time_fft_result_12d), intent(out) :: result !! spectrum
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: p1 !! exclusive upper time slice
        integer, intent(in), optional :: k_offset !! input sample offset
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! time-series axis, 1-based or negative
        call stft_nd11_kernel(self, cmplx(x,0.0_dp,dp), result, axis, p0, p1, &
            k_offset, padding, .false.)
    end subroutine stft_nd11_real

    subroutine stft_detrend_nd11_real(self, x, result, detrend, p0, p1, k_offset, padding, axis)
        class(short_time_fft), intent(in) :: self !! configured STFT plan
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:) !! input time series
        type(short_time_fft_result_12d), intent(out) :: result !! spectrum
        character(*), intent(in) :: detrend !! constant or linear segment detrending
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: p1 !! exclusive upper time slice
        integer, intent(in), optional :: k_offset !! input sample offset
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! time-series axis, 1-based or negative
        call stft_nd11_kernel(self, cmplx(x,0.0_dp,dp), result, axis, p0, p1, &
            k_offset, padding, .false., detrend)
    end subroutine stft_detrend_nd11_real

    subroutine stft_nd11_complex(self, x, result, p0, p1, k_offset, padding, axis)
        class(short_time_fft), intent(in) :: self !! configured STFT plan
        complex(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:) !! input time series
        type(short_time_fft_result_12d), intent(out) :: result !! spectrum
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: p1 !! exclusive upper time slice
        integer, intent(in), optional :: k_offset !! input sample offset
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! time-series axis, 1-based or negative
        call stft_nd11_kernel(self, x, result, axis, p0, p1, &
            k_offset, padding, .true.)
    end subroutine stft_nd11_complex

    subroutine stft_detrend_nd11_complex(self, x, result, detrend, p0, p1, k_offset, padding, axis)
        class(short_time_fft), intent(in) :: self !! configured STFT plan
        complex(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:) !! input time series
        type(short_time_fft_result_12d), intent(out) :: result !! spectrum
        character(*), intent(in) :: detrend !! constant or linear segment detrending
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: p1 !! exclusive upper time slice
        integer, intent(in), optional :: k_offset !! input sample offset
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! time-series axis, 1-based or negative
        call stft_nd11_kernel(self, x, result, axis, p0, p1, &
            k_offset, padding, .true., detrend)
    end subroutine stft_detrend_nd11_complex

    ! Flatten only the batch-index traversal: each one-dimensional signal is
    ! transformed by the same verified STFT kernel as ranks one through three.
    subroutine stft_nd11_kernel(self, x, result, axis, p0, p1, k_offset, padding, complex_input, detrend)
        class(short_time_fft), intent(in) :: self !! configured STFT plan
        complex(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:) !! input signals
        type(short_time_fft_result_12d), intent(out) :: result !! output spectrum
        integer, intent(in), optional :: axis !! signal axis
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: p1 !! exclusive upper time slice
        integer, intent(in), optional :: k_offset !! sample offset
        character(*), intent(in), optional :: padding !! extension mode
        logical, intent(in) :: complex_input !! whether input is complex
        character(*), intent(in), optional :: detrend !! detrending mode
        type(short_time_fft_result) :: frame
        complex(dp), allocatable :: trace(:)
        integer :: a, dims(11), dout(12), idx(11), oidx(12)
        integer :: batch_axes(10), b, ib, j, k, q, t, nb

        result%status = short_fft_invalid_input
        a = fft_resolve_axis(axis,11,11)
        if (a == 0) return
        dims = shape(x)
        if (any(dims < 1)) return
        batch_axes = pack([1,2,3,4,5,6,7,8,9,10,11],[1,2,3,4,5,6,7,8,9,10,11] /= a)
        nb = product(dims(batch_axes))
        allocate(trace(dims(a)))
        do ib=0,nb-1
            b = ib
            do j=1,10
                idx(batch_axes(j)) = modulo(b,dims(batch_axes(j)))+1
                b = b/dims(batch_axes(j))
            end do
            do k=1,dims(a)
                idx(a) = k
                trace(k) = x(idx(1),idx(2),idx(3),idx(4),idx(5),idx(6),idx(7),idx(8),idx(9),idx(10),idx(11))
            end do
            call stft_kernel(self,trace,frame,p0,p1,k_offset,padding,complex_input,detrend)
            if (frame%status /= short_fft_success) return
            if (ib == 0) then
                dout(1:11) = dims
                dout(a) = size(frame%values,1)
                dout(12) = size(frame%values,2)
                allocate(result%values(dout(1),dout(2),dout(3),dout(4),dout(5),dout(6),dout(7),dout(8), &
            dout(9),dout(10),dout(11),dout(12)))
                result%frequencies = frame%frequencies
                result%times = frame%times
                result%p0 = frame%p0
            end if
            do t=1,size(frame%values,2)
                do q=1,size(frame%values,1)
                    oidx(1:11) = idx
                    oidx(a) = q
                    oidx(12) = t
                    result%values(oidx(1),oidx(2),oidx(3),oidx(4),oidx(5),oidx(6),oidx(7),oidx(8),oidx(9), &
            oidx(10),oidx(11),oidx(12)) = frame%values(q,t)
                end do
            end do
        end do
        result%status = short_fft_success
    end subroutine stft_nd11_kernel

    subroutine istft_nd12_real(self, s, x, status, k0, k1, p0, f_axis, t_axis)
        class(short_time_fft), intent(in) :: self !! inverse plan
        complex(dp), intent(in) :: s(:,:,:,:,:,:,:,:,:,:,:,:) !! STFT coefficients
        real(dp), allocatable, intent(out) :: x(:,:,:,:,:,:,:,:,:,:,:) !! reconstructed signal
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: k0 !! inclusive first sample
        integer, intent(in), optional :: k1 !! exclusive last sample
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: f_axis !! frequency axis of S
        integer, intent(in), optional :: t_axis !! time axis of S
        complex(dp), allocatable :: z(:,:,:,:,:,:,:,:,:,:,:)

        allocate(x(0,0,0,0,0,0,0,0,0,0,0))
        status = short_fft_invalid_input
        if (self%fft_mode /= 'onesided' .and. self%fft_mode /= 'onesided2X') return
        call inverse_nd12_kernel(self,s,z,status,k0,k1,p0,f_axis,t_axis)
        if (status == short_fft_success) x = real(z,dp)
    end subroutine istft_nd12_real

    subroutine istft_nd12_complex(self, s, x, status, k0, k1, p0, f_axis, t_axis)
        class(short_time_fft), intent(in) :: self !! inverse plan
        complex(dp), intent(in) :: s(:,:,:,:,:,:,:,:,:,:,:,:) !! STFT coefficients
        complex(dp), allocatable, intent(out) :: x(:,:,:,:,:,:,:,:,:,:,:) !! reconstructed signal
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: k0 !! inclusive first sample
        integer, intent(in), optional :: k1 !! exclusive last sample
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: f_axis !! frequency axis of S
        integer, intent(in), optional :: t_axis !! time axis of S
        complex(dp), allocatable :: z(:,:,:,:,:,:,:,:,:,:,:)

        allocate(x(0,0,0,0,0,0,0,0,0,0,0))
        status = short_fft_invalid_input
        if (self%fft_mode == 'onesided' .or. self%fft_mode == 'onesided2X') return
        call inverse_nd12_kernel(self,s,z,status,k0,k1,p0,f_axis,t_axis)
        if (status == short_fft_success) x = z
    end subroutine istft_nd12_complex

    subroutine inverse_nd12_kernel(self, s, x, status, k0, k1, p0, f_axis, t_axis)
        class(short_time_fft), intent(in) :: self !! inverse transform plan
        complex(dp), intent(in) :: s(:,:,:,:,:,:,:,:,:,:,:,:) !! STFT coefficients
        complex(dp), allocatable, intent(out) :: x(:,:,:,:,:,:,:,:,:,:,:) !! reconstructed series
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: k0 !! inclusive first sample
        integer, intent(in), optional :: k1 !! exclusive last sample
        integer, intent(in), optional :: p0 !! first time slice
        integer, intent(in), optional :: f_axis !! frequency dimension
        integer, intent(in), optional :: t_axis !! time dimension
        complex(dp), allocatable :: frame(:,:), trace(:)
        integer :: f, t, dims(12), dout(11), batch_axes(10), out_axes(11)
        integer :: idx(12), oidx(11), b, ib, j, q, p, k, nb

        allocate(x(0,0,0,0,0,0,0,0,0,0,0))
        status = short_fft_invalid_input
        f = fft_resolve_axis(f_axis,12,11)
        t = fft_resolve_axis(t_axis,12,12)
        if (f == 0 .or. t == 0 .or. f == t) return
        dims = shape(s)
        if (any(dims < 1)) return
        batch_axes = pack([1,2,3,4,5,6,7,8,9,10,11,12],[1,2,3,4,5,6,7,8,9,10,11,12] /= f .and. [1,2,3,4,5,6,7,8,9,10,11,12] /= t)
        out_axes = pack([1,2,3,4,5,6,7,8,9,10,11,12],[1,2,3,4,5,6,7,8,9,10,11,12] /= t)
        nb = product(dims(batch_axes))
        allocate(frame(dims(f),dims(t)))
        do ib=0,nb-1
            b = ib
            do j=1,10
                idx(batch_axes(j)) = modulo(b,dims(batch_axes(j)))+1
                b = b/dims(batch_axes(j))
            end do
            do p=1,dims(t)
                idx(t) = p
                do q=1,dims(f)
                    idx(f) = q
                    frame(q,p) = s(idx(1),idx(2),idx(3),idx(4),idx(5),idx(6),idx(7),idx(8),idx(9),idx(10),idx(11),idx(12))
                end do
            end do
            call inverse_kernel(self,frame,trace,status,k0,k1,p0)
            if (status /= short_fft_success) return
            if (ib == 0) then
                dout = pack(dims,[1,2,3,4,5,6,7,8,9,10,11,12] /= t)
                dout(findloc(out_axes,f,dim=1)) = size(trace)
                deallocate(x)
                allocate(x(dout(1),dout(2),dout(3),dout(4),dout(5),dout(6),dout(7),dout(8),dout(9),dout(10),dout(11)))
            end if
            do j=1,10
                oidx(findloc(out_axes,batch_axes(j),dim=1)) = idx(batch_axes(j))
            end do
            do k=1,size(trace)
                oidx(findloc(out_axes,f,dim=1)) = k
                x(oidx(1),oidx(2),oidx(3),oidx(4),oidx(5),oidx(6),oidx(7),oidx(8),oidx(9),oidx(10),oidx(11)) = trace(k)
            end do
        end do
        status = short_fft_success
    end subroutine inverse_nd12_kernel

    subroutine spectrogram_nd11_real(self, x, result, p0, p1, padding, axis)
        class(short_time_fft), intent(in) :: self !! transform plan
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:) !! signal array
        type(short_time_spectrogram_result_12d), intent(out) :: result !! power spectrum
        integer, intent(in), optional :: p0 !! first time slice
        integer, intent(in), optional :: p1 !! exclusive final slice
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! selected sample dimension
        type(short_time_fft_result_12d) :: spectrum

        call self%stft(x,spectrum,p0=p0,p1=p1,padding=padding,axis=axis)
        result%status = spectrum%status
        if (spectrum%status /= short_fft_success) return
        result%power = abs(spectrum%values)**2
        result%frequencies = spectrum%frequencies
        result%times = spectrum%times
        result%p0 = spectrum%p0
    end subroutine spectrogram_nd11_real

    subroutine spectrogram_nd11_complex(self, x, result, p0, p1, padding, axis)
        class(short_time_fft), intent(in) :: self !! transform plan
        complex(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:) !! signal array
        type(short_time_spectrogram_result_12d), intent(out) :: result !! power spectrum
        integer, intent(in), optional :: p0 !! first time slice
        integer, intent(in), optional :: p1 !! exclusive final slice
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! selected sample dimension
        type(short_time_fft_result_12d) :: spectrum

        call self%stft(x,spectrum,p0=p0,p1=p1,padding=padding,axis=axis)
        result%status = spectrum%status
        if (spectrum%status /= short_fft_success) return
        result%power = abs(spectrum%values)**2
        result%frequencies = spectrum%frequencies
        result%times = spectrum%times
        result%p0 = spectrum%p0
    end subroutine spectrogram_nd11_complex

    subroutine cross_spectrogram_nd11_real(self, x, y, values, status, p0, p1, padding, axis)
        class(short_time_fft), intent(in) :: self !! transform plan
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:) !! first array of signals
        real(dp), intent(in) :: y(:,:,:,:,:,:,:,:,:,:,:) !! second array of signals
        complex(dp), allocatable, intent(out) :: values(:,:,:,:,:,:,:,:,:,:,:,:) !! cross spectrum
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: p0 !! first time slice
        integer, intent(in), optional :: p1 !! exclusive final slice
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! selected sample dimension
        type(short_time_fft_result_12d) :: sx, sy

        allocate(values(0,0,0,0,0,0,0,0,0,0,0,0))
        status = short_fft_invalid_input
        if (any(shape(x) /= shape(y))) return
        call self%stft(x,sx,p0=p0,p1=p1,padding=padding,axis=axis)
        if (sx%status /= short_fft_success) return
        call self%stft(y,sy,p0=p0,p1=p1,padding=padding,axis=axis)
        if (sy%status /= short_fft_success) return
        values = sx%values*conjg(sy%values)
        status = short_fft_success
    end subroutine cross_spectrogram_nd11_real

    subroutine cross_spectrogram_nd11_complex(self, x, y, values, status, p0, p1, padding, axis)
        class(short_time_fft), intent(in) :: self !! transform plan
        complex(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:) !! first array of signals
        complex(dp), intent(in) :: y(:,:,:,:,:,:,:,:,:,:,:) !! second array of signals
        complex(dp), allocatable, intent(out) :: values(:,:,:,:,:,:,:,:,:,:,:,:) !! cross spectrum
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: p0 !! first time slice
        integer, intent(in), optional :: p1 !! exclusive final slice
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! selected sample dimension
        type(short_time_fft_result_12d) :: sx, sy

        allocate(values(0,0,0,0,0,0,0,0,0,0,0,0))
        status = short_fft_invalid_input
        if (any(shape(x) /= shape(y))) return
        call self%stft(x,sx,p0=p0,p1=p1,padding=padding,axis=axis)
        if (sx%status /= short_fft_success) return
        call self%stft(y,sy,p0=p0,p1=p1,padding=padding,axis=axis)
        if (sy%status /= short_fft_success) return
        values = sx%values*conjg(sy%values)
        status = short_fft_success
    end subroutine cross_spectrogram_nd11_complex

    subroutine stft_nd12_real(self, x, result, p0, p1, k_offset, padding, axis)
        class(short_time_fft), intent(in) :: self !! configured STFT plan
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:) !! input time series
        type(short_time_fft_result_13d), intent(out) :: result !! spectrum
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: p1 !! exclusive upper time slice
        integer, intent(in), optional :: k_offset !! input sample offset
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! time-series axis, 1-based or negative
        call stft_nd12_kernel(self, cmplx(x,0.0_dp,dp), result, axis, p0, p1, &
            k_offset, padding, .false.)
    end subroutine stft_nd12_real

    subroutine stft_detrend_nd12_real(self, x, result, detrend, p0, p1, k_offset, padding, axis)
        class(short_time_fft), intent(in) :: self !! configured STFT plan
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:) !! input time series
        type(short_time_fft_result_13d), intent(out) :: result !! spectrum
        character(*), intent(in) :: detrend !! constant or linear segment detrending
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: p1 !! exclusive upper time slice
        integer, intent(in), optional :: k_offset !! input sample offset
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! time-series axis, 1-based or negative
        call stft_nd12_kernel(self, cmplx(x,0.0_dp,dp), result, axis, p0, p1, &
            k_offset, padding, .false., detrend)
    end subroutine stft_detrend_nd12_real

    subroutine stft_nd12_complex(self, x, result, p0, p1, k_offset, padding, axis)
        class(short_time_fft), intent(in) :: self !! configured STFT plan
        complex(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:) !! input time series
        type(short_time_fft_result_13d), intent(out) :: result !! spectrum
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: p1 !! exclusive upper time slice
        integer, intent(in), optional :: k_offset !! input sample offset
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! time-series axis, 1-based or negative
        call stft_nd12_kernel(self, x, result, axis, p0, p1, &
            k_offset, padding, .true.)
    end subroutine stft_nd12_complex

    subroutine stft_detrend_nd12_complex(self, x, result, detrend, p0, p1, k_offset, padding, axis)
        class(short_time_fft), intent(in) :: self !! configured STFT plan
        complex(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:) !! input time series
        type(short_time_fft_result_13d), intent(out) :: result !! spectrum
        character(*), intent(in) :: detrend !! constant or linear segment detrending
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: p1 !! exclusive upper time slice
        integer, intent(in), optional :: k_offset !! input sample offset
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! time-series axis, 1-based or negative
        call stft_nd12_kernel(self, x, result, axis, p0, p1, &
            k_offset, padding, .true., detrend)
    end subroutine stft_detrend_nd12_complex

    ! Flatten only the batch-index traversal: each one-dimensional signal is
    ! transformed by the same verified STFT kernel as ranks one through three.
    subroutine stft_nd12_kernel(self, x, result, axis, p0, p1, k_offset, padding, complex_input, detrend)
        class(short_time_fft), intent(in) :: self !! configured STFT plan
        complex(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:) !! input signals
        type(short_time_fft_result_13d), intent(out) :: result !! output spectrum
        integer, intent(in), optional :: axis !! signal axis
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: p1 !! exclusive upper time slice
        integer, intent(in), optional :: k_offset !! sample offset
        character(*), intent(in), optional :: padding !! extension mode
        logical, intent(in) :: complex_input !! whether input is complex
        character(*), intent(in), optional :: detrend !! detrending mode
        type(short_time_fft_result) :: frame
        complex(dp), allocatable :: trace(:)
        integer :: a, dims(12), dout(13), idx(12), oidx(13)
        integer :: batch_axes(11), b, ib, j, k, q, t, nb

        result%status = short_fft_invalid_input
        a = fft_resolve_axis(axis,12,12)
        if (a == 0) return
        dims = shape(x)
        if (any(dims < 1)) return
        batch_axes = pack([1,2,3,4,5,6,7,8,9,10,11,12],[1,2,3,4,5,6,7,8,9,10,11,12] /= a)
        nb = product(dims(batch_axes))
        allocate(trace(dims(a)))
        do ib=0,nb-1
            b = ib
            do j=1,11
                idx(batch_axes(j)) = modulo(b,dims(batch_axes(j)))+1
                b = b/dims(batch_axes(j))
            end do
            do k=1,dims(a)
                idx(a) = k
                trace(k) = x(idx(1),idx(2),idx(3),idx(4),idx(5),idx(6),idx(7),idx(8),idx(9),idx(10),idx(11),idx(12))
            end do
            call stft_kernel(self,trace,frame,p0,p1,k_offset,padding,complex_input,detrend)
            if (frame%status /= short_fft_success) return
            if (ib == 0) then
                dout(1:12) = dims
                dout(a) = size(frame%values,1)
                dout(13) = size(frame%values,2)
                allocate(result%values(dout(1),dout(2),dout(3),dout(4),dout(5),dout(6),dout(7),dout(8), &
            dout(9),dout(10),dout(11),dout(12),dout(13)))
                result%frequencies = frame%frequencies
                result%times = frame%times
                result%p0 = frame%p0
            end if
            do t=1,size(frame%values,2)
                do q=1,size(frame%values,1)
                    oidx(1:12) = idx
                    oidx(a) = q
                    oidx(13) = t
                    result%values(oidx(1),oidx(2),oidx(3),oidx(4),oidx(5),oidx(6),oidx(7),oidx(8),oidx(9), &
            oidx(10),oidx(11),oidx(12),oidx(13)) = frame%values(q,t)
                end do
            end do
        end do
        result%status = short_fft_success
    end subroutine stft_nd12_kernel

    subroutine istft_nd13_real(self, s, x, status, k0, k1, p0, f_axis, t_axis)
        class(short_time_fft), intent(in) :: self !! inverse plan
        complex(dp), intent(in) :: s(:,:,:,:,:,:,:,:,:,:,:,:,:) !! STFT coefficients
        real(dp), allocatable, intent(out) :: x(:,:,:,:,:,:,:,:,:,:,:,:) !! reconstructed signal
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: k0 !! inclusive first sample
        integer, intent(in), optional :: k1 !! exclusive last sample
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: f_axis !! frequency axis of S
        integer, intent(in), optional :: t_axis !! time axis of S
        complex(dp), allocatable :: z(:,:,:,:,:,:,:,:,:,:,:,:)

        allocate(x(0,0,0,0,0,0,0,0,0,0,0,0))
        status = short_fft_invalid_input
        if (self%fft_mode /= 'onesided' .and. self%fft_mode /= 'onesided2X') return
        call inverse_nd13_kernel(self,s,z,status,k0,k1,p0,f_axis,t_axis)
        if (status == short_fft_success) x = real(z,dp)
    end subroutine istft_nd13_real

    subroutine istft_nd13_complex(self, s, x, status, k0, k1, p0, f_axis, t_axis)
        class(short_time_fft), intent(in) :: self !! inverse plan
        complex(dp), intent(in) :: s(:,:,:,:,:,:,:,:,:,:,:,:,:) !! STFT coefficients
        complex(dp), allocatable, intent(out) :: x(:,:,:,:,:,:,:,:,:,:,:,:) !! reconstructed signal
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: k0 !! inclusive first sample
        integer, intent(in), optional :: k1 !! exclusive last sample
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: f_axis !! frequency axis of S
        integer, intent(in), optional :: t_axis !! time axis of S
        complex(dp), allocatable :: z(:,:,:,:,:,:,:,:,:,:,:,:)

        allocate(x(0,0,0,0,0,0,0,0,0,0,0,0))
        status = short_fft_invalid_input
        if (self%fft_mode == 'onesided' .or. self%fft_mode == 'onesided2X') return
        call inverse_nd13_kernel(self,s,z,status,k0,k1,p0,f_axis,t_axis)
        if (status == short_fft_success) x = z
    end subroutine istft_nd13_complex

    subroutine inverse_nd13_kernel(self, s, x, status, k0, k1, p0, f_axis, t_axis)
        class(short_time_fft), intent(in) :: self !! inverse transform plan
        complex(dp), intent(in) :: s(:,:,:,:,:,:,:,:,:,:,:,:,:) !! STFT coefficients
        complex(dp), allocatable, intent(out) :: x(:,:,:,:,:,:,:,:,:,:,:,:) !! reconstructed series
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: k0 !! inclusive first sample
        integer, intent(in), optional :: k1 !! exclusive last sample
        integer, intent(in), optional :: p0 !! first time slice
        integer, intent(in), optional :: f_axis !! frequency dimension
        integer, intent(in), optional :: t_axis !! time dimension
        complex(dp), allocatable :: frame(:,:), trace(:)
        integer :: f, t, dims(13), dout(12), batch_axes(11), out_axes(12)
        integer :: idx(13), oidx(12), b, ib, j, q, p, k, nb

        allocate(x(0,0,0,0,0,0,0,0,0,0,0,0))
        status = short_fft_invalid_input
        f = fft_resolve_axis(f_axis,13,12)
        t = fft_resolve_axis(t_axis,13,13)
        if (f == 0 .or. t == 0 .or. f == t) return
        dims = shape(s)
        if (any(dims < 1)) return
        batch_axes = pack([1,2,3,4,5,6,7,8,9,10,11,12,13],[1,2,3,4,5,6,7,8,9,10,11,12,13] /= f .and. [1,2,3, &
            4,5,6,7,8,9,10,11,12,13] /= t)
        out_axes = pack([1,2,3,4,5,6,7,8,9,10,11,12,13],[1,2,3,4,5,6,7,8,9,10,11,12,13] /= t)
        nb = product(dims(batch_axes))
        allocate(frame(dims(f),dims(t)))
        do ib=0,nb-1
            b = ib
            do j=1,11
                idx(batch_axes(j)) = modulo(b,dims(batch_axes(j)))+1
                b = b/dims(batch_axes(j))
            end do
            do p=1,dims(t)
                idx(t) = p
                do q=1,dims(f)
                    idx(f) = q
                    frame(q,p) = s(idx(1),idx(2),idx(3),idx(4),idx(5),idx(6),idx(7),idx(8),idx(9),idx(10),idx(11),idx(12),idx(13))
                end do
            end do
            call inverse_kernel(self,frame,trace,status,k0,k1,p0)
            if (status /= short_fft_success) return
            if (ib == 0) then
                dout = pack(dims,[1,2,3,4,5,6,7,8,9,10,11,12,13] /= t)
                dout(findloc(out_axes,f,dim=1)) = size(trace)
                deallocate(x)
                allocate(x(dout(1),dout(2),dout(3),dout(4),dout(5),dout(6),dout(7),dout(8),dout(9),dout(10),dout(11),dout(12)))
            end if
            do j=1,11
                oidx(findloc(out_axes,batch_axes(j),dim=1)) = idx(batch_axes(j))
            end do
            do k=1,size(trace)
                oidx(findloc(out_axes,f,dim=1)) = k
                x(oidx(1),oidx(2),oidx(3),oidx(4),oidx(5),oidx(6),oidx(7),oidx(8),oidx(9),oidx(10),oidx(11),oidx(12)) = trace(k)
            end do
        end do
        status = short_fft_success
    end subroutine inverse_nd13_kernel

    subroutine spectrogram_nd12_real(self, x, result, p0, p1, padding, axis)
        class(short_time_fft), intent(in) :: self !! transform plan
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:) !! signal array
        type(short_time_spectrogram_result_13d), intent(out) :: result !! power spectrum
        integer, intent(in), optional :: p0 !! first time slice
        integer, intent(in), optional :: p1 !! exclusive final slice
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! selected sample dimension
        type(short_time_fft_result_13d) :: spectrum

        call self%stft(x,spectrum,p0=p0,p1=p1,padding=padding,axis=axis)
        result%status = spectrum%status
        if (spectrum%status /= short_fft_success) return
        result%power = abs(spectrum%values)**2
        result%frequencies = spectrum%frequencies
        result%times = spectrum%times
        result%p0 = spectrum%p0
    end subroutine spectrogram_nd12_real

    subroutine spectrogram_nd12_complex(self, x, result, p0, p1, padding, axis)
        class(short_time_fft), intent(in) :: self !! transform plan
        complex(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:) !! signal array
        type(short_time_spectrogram_result_13d), intent(out) :: result !! power spectrum
        integer, intent(in), optional :: p0 !! first time slice
        integer, intent(in), optional :: p1 !! exclusive final slice
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! selected sample dimension
        type(short_time_fft_result_13d) :: spectrum

        call self%stft(x,spectrum,p0=p0,p1=p1,padding=padding,axis=axis)
        result%status = spectrum%status
        if (spectrum%status /= short_fft_success) return
        result%power = abs(spectrum%values)**2
        result%frequencies = spectrum%frequencies
        result%times = spectrum%times
        result%p0 = spectrum%p0
    end subroutine spectrogram_nd12_complex

    subroutine cross_spectrogram_nd12_real(self, x, y, values, status, p0, p1, padding, axis)
        class(short_time_fft), intent(in) :: self !! transform plan
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:) !! first array of signals
        real(dp), intent(in) :: y(:,:,:,:,:,:,:,:,:,:,:,:) !! second array of signals
        complex(dp), allocatable, intent(out) :: values(:,:,:,:,:,:,:,:,:,:,:,:,:) !! cross spectrum
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: p0 !! first time slice
        integer, intent(in), optional :: p1 !! exclusive final slice
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! selected sample dimension
        type(short_time_fft_result_13d) :: sx, sy

        allocate(values(0,0,0,0,0,0,0,0,0,0,0,0,0))
        status = short_fft_invalid_input
        if (any(shape(x) /= shape(y))) return
        call self%stft(x,sx,p0=p0,p1=p1,padding=padding,axis=axis)
        if (sx%status /= short_fft_success) return
        call self%stft(y,sy,p0=p0,p1=p1,padding=padding,axis=axis)
        if (sy%status /= short_fft_success) return
        values = sx%values*conjg(sy%values)
        status = short_fft_success
    end subroutine cross_spectrogram_nd12_real

    subroutine cross_spectrogram_nd12_complex(self, x, y, values, status, p0, p1, padding, axis)
        class(short_time_fft), intent(in) :: self !! transform plan
        complex(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:) !! first array of signals
        complex(dp), intent(in) :: y(:,:,:,:,:,:,:,:,:,:,:,:) !! second array of signals
        complex(dp), allocatable, intent(out) :: values(:,:,:,:,:,:,:,:,:,:,:,:,:) !! cross spectrum
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: p0 !! first time slice
        integer, intent(in), optional :: p1 !! exclusive final slice
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! selected sample dimension
        type(short_time_fft_result_13d) :: sx, sy

        allocate(values(0,0,0,0,0,0,0,0,0,0,0,0,0))
        status = short_fft_invalid_input
        if (any(shape(x) /= shape(y))) return
        call self%stft(x,sx,p0=p0,p1=p1,padding=padding,axis=axis)
        if (sx%status /= short_fft_success) return
        call self%stft(y,sy,p0=p0,p1=p1,padding=padding,axis=axis)
        if (sy%status /= short_fft_success) return
        values = sx%values*conjg(sy%values)
        status = short_fft_success
    end subroutine cross_spectrogram_nd12_complex

    subroutine stft_nd13_real(self, x, result, p0, p1, k_offset, padding, axis)
        class(short_time_fft), intent(in) :: self !! configured STFT plan
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:,:) !! input time series
        type(short_time_fft_result_14d), intent(out) :: result !! spectrum
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: p1 !! exclusive upper time slice
        integer, intent(in), optional :: k_offset !! input sample offset
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! time-series axis, 1-based or negative
        call stft_nd13_kernel(self, cmplx(x,0.0_dp,dp), result, axis, p0, p1, &
            k_offset, padding, .false.)
    end subroutine stft_nd13_real

    subroutine stft_detrend_nd13_real(self, x, result, detrend, p0, p1, k_offset, padding, axis)
        class(short_time_fft), intent(in) :: self !! configured STFT plan
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:,:) !! input time series
        type(short_time_fft_result_14d), intent(out) :: result !! spectrum
        character(*), intent(in) :: detrend !! constant or linear segment detrending
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: p1 !! exclusive upper time slice
        integer, intent(in), optional :: k_offset !! input sample offset
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! time-series axis, 1-based or negative
        call stft_nd13_kernel(self, cmplx(x,0.0_dp,dp), result, axis, p0, p1, &
            k_offset, padding, .false., detrend)
    end subroutine stft_detrend_nd13_real

    subroutine stft_nd13_complex(self, x, result, p0, p1, k_offset, padding, axis)
        class(short_time_fft), intent(in) :: self !! configured STFT plan
        complex(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:,:) !! input time series
        type(short_time_fft_result_14d), intent(out) :: result !! spectrum
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: p1 !! exclusive upper time slice
        integer, intent(in), optional :: k_offset !! input sample offset
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! time-series axis, 1-based or negative
        call stft_nd13_kernel(self, x, result, axis, p0, p1, &
            k_offset, padding, .true.)
    end subroutine stft_nd13_complex

    subroutine stft_detrend_nd13_complex(self, x, result, detrend, p0, p1, k_offset, padding, axis)
        class(short_time_fft), intent(in) :: self !! configured STFT plan
        complex(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:,:) !! input time series
        type(short_time_fft_result_14d), intent(out) :: result !! spectrum
        character(*), intent(in) :: detrend !! constant or linear segment detrending
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: p1 !! exclusive upper time slice
        integer, intent(in), optional :: k_offset !! input sample offset
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! time-series axis, 1-based or negative
        call stft_nd13_kernel(self, x, result, axis, p0, p1, &
            k_offset, padding, .true., detrend)
    end subroutine stft_detrend_nd13_complex

    ! Flatten only the batch-index traversal: each one-dimensional signal is
    ! transformed by the same verified STFT kernel as ranks one through three.
    subroutine stft_nd13_kernel(self, x, result, axis, p0, p1, k_offset, padding, complex_input, detrend)
        class(short_time_fft), intent(in) :: self !! configured STFT plan
        complex(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:,:) !! input signals
        type(short_time_fft_result_14d), intent(out) :: result !! output spectrum
        integer, intent(in), optional :: axis !! signal axis
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: p1 !! exclusive upper time slice
        integer, intent(in), optional :: k_offset !! sample offset
        character(*), intent(in), optional :: padding !! extension mode
        logical, intent(in) :: complex_input !! whether input is complex
        character(*), intent(in), optional :: detrend !! detrending mode
        type(short_time_fft_result) :: frame
        complex(dp), allocatable :: trace(:)
        integer :: a, dims(13), dout(14), idx(13), oidx(14)
        integer :: batch_axes(12), b, ib, j, k, q, t, nb

        result%status = short_fft_invalid_input
        a = fft_resolve_axis(axis,13,13)
        if (a == 0) return
        dims = shape(x)
        if (any(dims < 1)) return
        batch_axes = pack([1,2,3,4,5,6,7,8,9,10,11,12,13],[1,2,3,4,5,6,7,8,9,10,11,12,13] /= a)
        nb = product(dims(batch_axes))
        allocate(trace(dims(a)))
        do ib=0,nb-1
            b = ib
            do j=1,12
                idx(batch_axes(j)) = modulo(b,dims(batch_axes(j)))+1
                b = b/dims(batch_axes(j))
            end do
            do k=1,dims(a)
                idx(a) = k
                trace(k) = x(idx(1),idx(2),idx(3),idx(4),idx(5),idx(6),idx(7),idx(8),idx(9),idx(10),idx(11),idx(12),idx(13))
            end do
            call stft_kernel(self,trace,frame,p0,p1,k_offset,padding,complex_input,detrend)
            if (frame%status /= short_fft_success) return
            if (ib == 0) then
                dout(1:13) = dims
                dout(a) = size(frame%values,1)
                dout(14) = size(frame%values,2)
                allocate(result%values(dout(1),dout(2),dout(3),dout(4),dout(5),dout(6),dout(7),dout(8), &
            dout(9),dout(10),dout(11),dout(12),dout(13),dout(14)))
                result%frequencies = frame%frequencies
                result%times = frame%times
                result%p0 = frame%p0
            end if
            do t=1,size(frame%values,2)
                do q=1,size(frame%values,1)
                    oidx(1:13) = idx
                    oidx(a) = q
                    oidx(14) = t
                    result%values(oidx(1),oidx(2),oidx(3),oidx(4),oidx(5),oidx(6),oidx(7),oidx(8),oidx(9), &
            oidx(10),oidx(11),oidx(12),oidx(13),oidx(14)) = frame%values(q,t)
                end do
            end do
        end do
        result%status = short_fft_success
    end subroutine stft_nd13_kernel

    subroutine istft_nd14_real(self, s, x, status, k0, k1, p0, f_axis, t_axis)
        class(short_time_fft), intent(in) :: self !! inverse plan
        complex(dp), intent(in) :: s(:,:,:,:,:,:,:,:,:,:,:,:,:,:) !! STFT coefficients
        real(dp), allocatable, intent(out) :: x(:,:,:,:,:,:,:,:,:,:,:,:,:) !! reconstructed signal
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: k0 !! inclusive first sample
        integer, intent(in), optional :: k1 !! exclusive last sample
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: f_axis !! frequency axis of S
        integer, intent(in), optional :: t_axis !! time axis of S
        complex(dp), allocatable :: z(:,:,:,:,:,:,:,:,:,:,:,:,:)

        allocate(x(0,0,0,0,0,0,0,0,0,0,0,0,0))
        status = short_fft_invalid_input
        if (self%fft_mode /= 'onesided' .and. self%fft_mode /= 'onesided2X') return
        call inverse_nd14_kernel(self,s,z,status,k0,k1,p0,f_axis,t_axis)
        if (status == short_fft_success) x = real(z,dp)
    end subroutine istft_nd14_real

    subroutine istft_nd14_complex(self, s, x, status, k0, k1, p0, f_axis, t_axis)
        class(short_time_fft), intent(in) :: self !! inverse plan
        complex(dp), intent(in) :: s(:,:,:,:,:,:,:,:,:,:,:,:,:,:) !! STFT coefficients
        complex(dp), allocatable, intent(out) :: x(:,:,:,:,:,:,:,:,:,:,:,:,:) !! reconstructed signal
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: k0 !! inclusive first sample
        integer, intent(in), optional :: k1 !! exclusive last sample
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: f_axis !! frequency axis of S
        integer, intent(in), optional :: t_axis !! time axis of S
        complex(dp), allocatable :: z(:,:,:,:,:,:,:,:,:,:,:,:,:)

        allocate(x(0,0,0,0,0,0,0,0,0,0,0,0,0))
        status = short_fft_invalid_input
        if (self%fft_mode == 'onesided' .or. self%fft_mode == 'onesided2X') return
        call inverse_nd14_kernel(self,s,z,status,k0,k1,p0,f_axis,t_axis)
        if (status == short_fft_success) x = z
    end subroutine istft_nd14_complex

    subroutine inverse_nd14_kernel(self, s, x, status, k0, k1, p0, f_axis, t_axis)
        class(short_time_fft), intent(in) :: self !! inverse transform plan
        complex(dp), intent(in) :: s(:,:,:,:,:,:,:,:,:,:,:,:,:,:) !! STFT coefficients
        complex(dp), allocatable, intent(out) :: x(:,:,:,:,:,:,:,:,:,:,:,:,:) !! reconstructed series
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: k0 !! inclusive first sample
        integer, intent(in), optional :: k1 !! exclusive last sample
        integer, intent(in), optional :: p0 !! first time slice
        integer, intent(in), optional :: f_axis !! frequency dimension
        integer, intent(in), optional :: t_axis !! time dimension
        complex(dp), allocatable :: frame(:,:), trace(:)
        integer :: f, t, dims(14), dout(13), batch_axes(12), out_axes(13)
        integer :: idx(14), oidx(13), b, ib, j, q, p, k, nb

        allocate(x(0,0,0,0,0,0,0,0,0,0,0,0,0))
        status = short_fft_invalid_input
        f = fft_resolve_axis(f_axis,14,13)
        t = fft_resolve_axis(t_axis,14,14)
        if (f == 0 .or. t == 0 .or. f == t) return
        dims = shape(s)
        if (any(dims < 1)) return
        batch_axes = pack([1,2,3,4,5,6,7,8,9,10,11,12,13,14],[1,2,3,4,5,6,7,8,9,10,11,12,13, &
            14] /= f .and. [1,2,3,4,5,6,7,8,9,10,11,12,13,14] /= t)
        out_axes = pack([1,2,3,4,5,6,7,8,9,10,11,12,13,14],[1,2,3,4,5,6,7,8,9,10,11,12,13,14] /= t)
        nb = product(dims(batch_axes))
        allocate(frame(dims(f),dims(t)))
        do ib=0,nb-1
            b = ib
            do j=1,12
                idx(batch_axes(j)) = modulo(b,dims(batch_axes(j)))+1
                b = b/dims(batch_axes(j))
            end do
            do p=1,dims(t)
                idx(t) = p
                do q=1,dims(f)
                    idx(f) = q
                    frame(q,p) = s(idx(1),idx(2),idx(3),idx(4),idx(5),idx(6),idx(7),idx(8),idx(9),idx(10), &
            idx(11),idx(12),idx(13),idx(14))
                end do
            end do
            call inverse_kernel(self,frame,trace,status,k0,k1,p0)
            if (status /= short_fft_success) return
            if (ib == 0) then
                dout = pack(dims,[1,2,3,4,5,6,7,8,9,10,11,12,13,14] /= t)
                dout(findloc(out_axes,f,dim=1)) = size(trace)
                deallocate(x)
                allocate(x(dout(1),dout(2),dout(3),dout(4),dout(5),dout(6),dout(7),dout(8),dout(9),dout(10), &
            dout(11),dout(12),dout(13)))
            end if
            do j=1,12
                oidx(findloc(out_axes,batch_axes(j),dim=1)) = idx(batch_axes(j))
            end do
            do k=1,size(trace)
                oidx(findloc(out_axes,f,dim=1)) = k
                x(oidx(1),oidx(2),oidx(3),oidx(4),oidx(5),oidx(6),oidx(7),oidx(8),oidx(9),oidx(10),oidx(11), &
            oidx(12),oidx(13)) = trace(k)
            end do
        end do
        status = short_fft_success
    end subroutine inverse_nd14_kernel

    subroutine spectrogram_nd13_real(self, x, result, p0, p1, padding, axis)
        class(short_time_fft), intent(in) :: self !! transform plan
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:,:) !! signal array
        type(short_time_spectrogram_result_14d), intent(out) :: result !! power spectrum
        integer, intent(in), optional :: p0 !! first time slice
        integer, intent(in), optional :: p1 !! exclusive final slice
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! selected sample dimension
        type(short_time_fft_result_14d) :: spectrum

        call self%stft(x,spectrum,p0=p0,p1=p1,padding=padding,axis=axis)
        result%status = spectrum%status
        if (spectrum%status /= short_fft_success) return
        result%power = abs(spectrum%values)**2
        result%frequencies = spectrum%frequencies
        result%times = spectrum%times
        result%p0 = spectrum%p0
    end subroutine spectrogram_nd13_real

    subroutine spectrogram_nd13_complex(self, x, result, p0, p1, padding, axis)
        class(short_time_fft), intent(in) :: self !! transform plan
        complex(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:,:) !! signal array
        type(short_time_spectrogram_result_14d), intent(out) :: result !! power spectrum
        integer, intent(in), optional :: p0 !! first time slice
        integer, intent(in), optional :: p1 !! exclusive final slice
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! selected sample dimension
        type(short_time_fft_result_14d) :: spectrum

        call self%stft(x,spectrum,p0=p0,p1=p1,padding=padding,axis=axis)
        result%status = spectrum%status
        if (spectrum%status /= short_fft_success) return
        result%power = abs(spectrum%values)**2
        result%frequencies = spectrum%frequencies
        result%times = spectrum%times
        result%p0 = spectrum%p0
    end subroutine spectrogram_nd13_complex

    subroutine cross_spectrogram_nd13_real(self, x, y, values, status, p0, p1, padding, axis)
        class(short_time_fft), intent(in) :: self !! transform plan
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:,:) !! first array of signals
        real(dp), intent(in) :: y(:,:,:,:,:,:,:,:,:,:,:,:,:) !! second array of signals
        complex(dp), allocatable, intent(out) :: values(:,:,:,:,:,:,:,:,:,:,:,:,:,:) !! cross spectrum
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: p0 !! first time slice
        integer, intent(in), optional :: p1 !! exclusive final slice
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! selected sample dimension
        type(short_time_fft_result_14d) :: sx, sy

        allocate(values(0,0,0,0,0,0,0,0,0,0,0,0,0,0))
        status = short_fft_invalid_input
        if (any(shape(x) /= shape(y))) return
        call self%stft(x,sx,p0=p0,p1=p1,padding=padding,axis=axis)
        if (sx%status /= short_fft_success) return
        call self%stft(y,sy,p0=p0,p1=p1,padding=padding,axis=axis)
        if (sy%status /= short_fft_success) return
        values = sx%values*conjg(sy%values)
        status = short_fft_success
    end subroutine cross_spectrogram_nd13_real

    subroutine cross_spectrogram_nd13_complex(self, x, y, values, status, p0, p1, padding, axis)
        class(short_time_fft), intent(in) :: self !! transform plan
        complex(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:,:) !! first array of signals
        complex(dp), intent(in) :: y(:,:,:,:,:,:,:,:,:,:,:,:,:) !! second array of signals
        complex(dp), allocatable, intent(out) :: values(:,:,:,:,:,:,:,:,:,:,:,:,:,:) !! cross spectrum
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: p0 !! first time slice
        integer, intent(in), optional :: p1 !! exclusive final slice
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! selected sample dimension
        type(short_time_fft_result_14d) :: sx, sy

        allocate(values(0,0,0,0,0,0,0,0,0,0,0,0,0,0))
        status = short_fft_invalid_input
        if (any(shape(x) /= shape(y))) return
        call self%stft(x,sx,p0=p0,p1=p1,padding=padding,axis=axis)
        if (sx%status /= short_fft_success) return
        call self%stft(y,sy,p0=p0,p1=p1,padding=padding,axis=axis)
        if (sy%status /= short_fft_success) return
        values = sx%values*conjg(sy%values)
        status = short_fft_success
    end subroutine cross_spectrogram_nd13_complex

    subroutine stft_nd14_real(self, x, result, p0, p1, k_offset, padding, axis)
        class(short_time_fft), intent(in) :: self !! configured STFT plan
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:,:,:) !! input time series
        type(short_time_fft_result_15d), intent(out) :: result !! spectrum
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: p1 !! exclusive upper time slice
        integer, intent(in), optional :: k_offset !! input sample offset
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! time-series axis, 1-based or negative
        call stft_nd14_kernel(self, cmplx(x,0.0_dp,dp), result, axis, p0, p1, &
            k_offset, padding, .false.)
    end subroutine stft_nd14_real

    subroutine stft_detrend_nd14_real(self, x, result, detrend, p0, p1, k_offset, padding, axis)
        class(short_time_fft), intent(in) :: self !! configured STFT plan
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:,:,:) !! input time series
        type(short_time_fft_result_15d), intent(out) :: result !! spectrum
        character(*), intent(in) :: detrend !! constant or linear segment detrending
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: p1 !! exclusive upper time slice
        integer, intent(in), optional :: k_offset !! input sample offset
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! time-series axis, 1-based or negative
        call stft_nd14_kernel(self, cmplx(x,0.0_dp,dp), result, axis, p0, p1, &
            k_offset, padding, .false., detrend)
    end subroutine stft_detrend_nd14_real

    subroutine stft_nd14_complex(self, x, result, p0, p1, k_offset, padding, axis)
        class(short_time_fft), intent(in) :: self !! configured STFT plan
        complex(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:,:,:) !! input time series
        type(short_time_fft_result_15d), intent(out) :: result !! spectrum
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: p1 !! exclusive upper time slice
        integer, intent(in), optional :: k_offset !! input sample offset
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! time-series axis, 1-based or negative
        call stft_nd14_kernel(self, x, result, axis, p0, p1, &
            k_offset, padding, .true.)
    end subroutine stft_nd14_complex

    subroutine stft_detrend_nd14_complex(self, x, result, detrend, p0, p1, k_offset, padding, axis)
        class(short_time_fft), intent(in) :: self !! configured STFT plan
        complex(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:,:,:) !! input time series
        type(short_time_fft_result_15d), intent(out) :: result !! spectrum
        character(*), intent(in) :: detrend !! constant or linear segment detrending
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: p1 !! exclusive upper time slice
        integer, intent(in), optional :: k_offset !! input sample offset
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! time-series axis, 1-based or negative
        call stft_nd14_kernel(self, x, result, axis, p0, p1, &
            k_offset, padding, .true., detrend)
    end subroutine stft_detrend_nd14_complex

    ! Flatten only the batch-index traversal: each one-dimensional signal is
    ! transformed by the same verified STFT kernel as ranks one through three.
    subroutine stft_nd14_kernel(self, x, result, axis, p0, p1, k_offset, padding, complex_input, detrend)
        class(short_time_fft), intent(in) :: self !! configured STFT plan
        complex(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:,:,:) !! input signals
        type(short_time_fft_result_15d), intent(out) :: result !! output spectrum
        integer, intent(in), optional :: axis !! signal axis
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: p1 !! exclusive upper time slice
        integer, intent(in), optional :: k_offset !! sample offset
        character(*), intent(in), optional :: padding !! extension mode
        logical, intent(in) :: complex_input !! whether input is complex
        character(*), intent(in), optional :: detrend !! detrending mode
        type(short_time_fft_result) :: frame
        complex(dp), allocatable :: trace(:)
        integer :: a, dims(14), dout(15), idx(14), oidx(15)
        integer :: batch_axes(13), b, ib, j, k, q, t, nb

        result%status = short_fft_invalid_input
        a = fft_resolve_axis(axis,14,14)
        if (a == 0) return
        dims = shape(x)
        if (any(dims < 1)) return
        batch_axes = pack([1,2,3,4,5,6,7,8,9,10,11,12,13,14],[1,2,3,4,5,6,7,8,9,10,11,12,13,14] /= a)
        nb = product(dims(batch_axes))
        allocate(trace(dims(a)))
        do ib=0,nb-1
            b = ib
            do j=1,13
                idx(batch_axes(j)) = modulo(b,dims(batch_axes(j)))+1
                b = b/dims(batch_axes(j))
            end do
            do k=1,dims(a)
                idx(a) = k
                trace(k) = x(idx(1),idx(2),idx(3),idx(4),idx(5),idx(6),idx(7),idx(8),idx(9),idx(10),idx(11),idx(12),idx(13),idx(14))
            end do
            call stft_kernel(self,trace,frame,p0,p1,k_offset,padding,complex_input,detrend)
            if (frame%status /= short_fft_success) return
            if (ib == 0) then
                dout(1:14) = dims
                dout(a) = size(frame%values,1)
                dout(15) = size(frame%values,2)
                allocate(result%values(dout(1),dout(2),dout(3),dout(4),dout(5),dout(6),dout(7),dout(8), &
            dout(9),dout(10),dout(11),dout(12),dout(13),dout(14),dout(15)))
                result%frequencies = frame%frequencies
                result%times = frame%times
                result%p0 = frame%p0
            end if
            do t=1,size(frame%values,2)
                do q=1,size(frame%values,1)
                    oidx(1:14) = idx
                    oidx(a) = q
                    oidx(15) = t
                    result%values(oidx(1),oidx(2),oidx(3),oidx(4),oidx(5),oidx(6),oidx(7),oidx(8),oidx(9), &
            oidx(10),oidx(11),oidx(12),oidx(13),oidx(14),oidx(15)) = frame%values(q,t)
                end do
            end do
        end do
        result%status = short_fft_success
    end subroutine stft_nd14_kernel

    subroutine istft_nd15_real(self, s, x, status, k0, k1, p0, f_axis, t_axis)
        class(short_time_fft), intent(in) :: self !! inverse plan
        complex(dp), intent(in) :: s(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:) !! STFT coefficients
        real(dp), allocatable, intent(out) :: x(:,:,:,:,:,:,:,:,:,:,:,:,:,:) !! reconstructed signal
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: k0 !! inclusive first sample
        integer, intent(in), optional :: k1 !! exclusive last sample
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: f_axis !! frequency axis of S
        integer, intent(in), optional :: t_axis !! time axis of S
        complex(dp), allocatable :: z(:,:,:,:,:,:,:,:,:,:,:,:,:,:)

        allocate(x(0,0,0,0,0,0,0,0,0,0,0,0,0,0))
        status = short_fft_invalid_input
        if (self%fft_mode /= 'onesided' .and. self%fft_mode /= 'onesided2X') return
        call inverse_nd15_kernel(self,s,z,status,k0,k1,p0,f_axis,t_axis)
        if (status == short_fft_success) x = real(z,dp)
    end subroutine istft_nd15_real

    subroutine istft_nd15_complex(self, s, x, status, k0, k1, p0, f_axis, t_axis)
        class(short_time_fft), intent(in) :: self !! inverse plan
        complex(dp), intent(in) :: s(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:) !! STFT coefficients
        complex(dp), allocatable, intent(out) :: x(:,:,:,:,:,:,:,:,:,:,:,:,:,:) !! reconstructed signal
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: k0 !! inclusive first sample
        integer, intent(in), optional :: k1 !! exclusive last sample
        integer, intent(in), optional :: p0 !! first zero-based time slice
        integer, intent(in), optional :: f_axis !! frequency axis of S
        integer, intent(in), optional :: t_axis !! time axis of S
        complex(dp), allocatable :: z(:,:,:,:,:,:,:,:,:,:,:,:,:,:)

        allocate(x(0,0,0,0,0,0,0,0,0,0,0,0,0,0))
        status = short_fft_invalid_input
        if (self%fft_mode == 'onesided' .or. self%fft_mode == 'onesided2X') return
        call inverse_nd15_kernel(self,s,z,status,k0,k1,p0,f_axis,t_axis)
        if (status == short_fft_success) x = z
    end subroutine istft_nd15_complex

    subroutine inverse_nd15_kernel(self, s, x, status, k0, k1, p0, f_axis, t_axis)
        class(short_time_fft), intent(in) :: self !! inverse transform plan
        complex(dp), intent(in) :: s(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:) !! STFT coefficients
        complex(dp), allocatable, intent(out) :: x(:,:,:,:,:,:,:,:,:,:,:,:,:,:) !! reconstructed series
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: k0 !! inclusive first sample
        integer, intent(in), optional :: k1 !! exclusive last sample
        integer, intent(in), optional :: p0 !! first time slice
        integer, intent(in), optional :: f_axis !! frequency dimension
        integer, intent(in), optional :: t_axis !! time dimension
        complex(dp), allocatable :: frame(:,:), trace(:)
        integer :: f, t, dims(15), dout(14), batch_axes(13), out_axes(14)
        integer :: idx(15), oidx(14), b, ib, j, q, p, k, nb

        allocate(x(0,0,0,0,0,0,0,0,0,0,0,0,0,0))
        status = short_fft_invalid_input
        f = fft_resolve_axis(f_axis,15,14)
        t = fft_resolve_axis(t_axis,15,15)
        if (f == 0 .or. t == 0 .or. f == t) return
        dims = shape(s)
        if (any(dims < 1)) return
        batch_axes = pack([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15],[1,2,3,4,5,6,7,8,9,10,11,12,13,14, &
            15] /= f .and. [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15] /= t)
        out_axes = pack([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15],[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15] /= t)
        nb = product(dims(batch_axes))
        allocate(frame(dims(f),dims(t)))
        do ib=0,nb-1
            b = ib
            do j=1,13
                idx(batch_axes(j)) = modulo(b,dims(batch_axes(j)))+1
                b = b/dims(batch_axes(j))
            end do
            do p=1,dims(t)
                idx(t) = p
                do q=1,dims(f)
                    idx(f) = q
                    frame(q,p) = s(idx(1),idx(2),idx(3),idx(4),idx(5),idx(6),idx(7),idx(8),idx(9),idx(10), &
            idx(11),idx(12),idx(13),idx(14),idx(15))
                end do
            end do
            call inverse_kernel(self,frame,trace,status,k0,k1,p0)
            if (status /= short_fft_success) return
            if (ib == 0) then
                dout = pack(dims,[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15] /= t)
                dout(findloc(out_axes,f,dim=1)) = size(trace)
                deallocate(x)
                allocate(x(dout(1),dout(2),dout(3),dout(4),dout(5),dout(6),dout(7),dout(8),dout(9),dout(10), &
            dout(11),dout(12),dout(13),dout(14)))
            end if
            do j=1,13
                oidx(findloc(out_axes,batch_axes(j),dim=1)) = idx(batch_axes(j))
            end do
            do k=1,size(trace)
                oidx(findloc(out_axes,f,dim=1)) = k
                x(oidx(1),oidx(2),oidx(3),oidx(4),oidx(5),oidx(6),oidx(7),oidx(8),oidx(9),oidx(10),oidx(11), &
            oidx(12),oidx(13),oidx(14)) = trace(k)
            end do
        end do
        status = short_fft_success
    end subroutine inverse_nd15_kernel

    subroutine spectrogram_nd14_real(self, x, result, p0, p1, padding, axis)
        class(short_time_fft), intent(in) :: self !! transform plan
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:,:,:) !! signal array
        type(short_time_spectrogram_result_15d), intent(out) :: result !! power spectrum
        integer, intent(in), optional :: p0 !! first time slice
        integer, intent(in), optional :: p1 !! exclusive final slice
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! selected sample dimension
        type(short_time_fft_result_15d) :: spectrum

        call self%stft(x,spectrum,p0=p0,p1=p1,padding=padding,axis=axis)
        result%status = spectrum%status
        if (spectrum%status /= short_fft_success) return
        result%power = abs(spectrum%values)**2
        result%frequencies = spectrum%frequencies
        result%times = spectrum%times
        result%p0 = spectrum%p0
    end subroutine spectrogram_nd14_real

    subroutine spectrogram_nd14_complex(self, x, result, p0, p1, padding, axis)
        class(short_time_fft), intent(in) :: self !! transform plan
        complex(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:,:,:) !! signal array
        type(short_time_spectrogram_result_15d), intent(out) :: result !! power spectrum
        integer, intent(in), optional :: p0 !! first time slice
        integer, intent(in), optional :: p1 !! exclusive final slice
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! selected sample dimension
        type(short_time_fft_result_15d) :: spectrum

        call self%stft(x,spectrum,p0=p0,p1=p1,padding=padding,axis=axis)
        result%status = spectrum%status
        if (spectrum%status /= short_fft_success) return
        result%power = abs(spectrum%values)**2
        result%frequencies = spectrum%frequencies
        result%times = spectrum%times
        result%p0 = spectrum%p0
    end subroutine spectrogram_nd14_complex

    subroutine cross_spectrogram_nd14_real(self, x, y, values, status, p0, p1, padding, axis)
        class(short_time_fft), intent(in) :: self !! transform plan
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:,:,:) !! first array of signals
        real(dp), intent(in) :: y(:,:,:,:,:,:,:,:,:,:,:,:,:,:) !! second array of signals
        complex(dp), allocatable, intent(out) :: values(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:) !! cross spectrum
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: p0 !! first time slice
        integer, intent(in), optional :: p1 !! exclusive final slice
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! selected sample dimension
        type(short_time_fft_result_15d) :: sx, sy

        allocate(values(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))
        status = short_fft_invalid_input
        if (any(shape(x) /= shape(y))) return
        call self%stft(x,sx,p0=p0,p1=p1,padding=padding,axis=axis)
        if (sx%status /= short_fft_success) return
        call self%stft(y,sy,p0=p0,p1=p1,padding=padding,axis=axis)
        if (sy%status /= short_fft_success) return
        values = sx%values*conjg(sy%values)
        status = short_fft_success
    end subroutine cross_spectrogram_nd14_real

    subroutine cross_spectrogram_nd14_complex(self, x, y, values, status, p0, p1, padding, axis)
        class(short_time_fft), intent(in) :: self !! transform plan
        complex(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:,:,:) !! first array of signals
        complex(dp), intent(in) :: y(:,:,:,:,:,:,:,:,:,:,:,:,:,:) !! second array of signals
        complex(dp), allocatable, intent(out) :: values(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:) !! cross spectrum
        integer, intent(out) :: status !! zero on success
        integer, intent(in), optional :: p0 !! first time slice
        integer, intent(in), optional :: p1 !! exclusive final slice
        character(*), intent(in), optional :: padding !! boundary extension
        integer, intent(in), optional :: axis !! selected sample dimension
        type(short_time_fft_result_15d) :: sx, sy

        allocate(values(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))
        status = short_fft_invalid_input
        if (any(shape(x) /= shape(y))) return
        call self%stft(x,sx,p0=p0,p1=p1,padding=padding,axis=axis)
        if (sx%status /= short_fft_success) return
        call self%stft(y,sy,p0=p0,p1=p1,padding=padding,axis=axis)
        if (sy%status /= short_fft_success) return
        values = sx%values*conjg(sy%values)
        status = short_fft_success
    end subroutine cross_spectrogram_nd14_complex

end module scifort_signal_short_time_fft
