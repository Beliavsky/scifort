! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Direct discrete linear convolution and cross-correlation in one dimension.
module scifort_signal_convolve
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_signal_fftconvolve, only : fftconvolve
    implicit none
    private

    integer, parameter, public :: convolution_status_success = 0
    integer, parameter, public :: convolution_status_invalid_input = 1

    interface convolve
        module procedure convolve_real
        module procedure convolve_complex
        module procedure convolve_real_complex, convolve_complex_real
    end interface convolve

    interface correlate
        module procedure correlate_real
        module procedure correlate_complex
        module procedure correlate_real_complex, correlate_complex_real
    end interface correlate

    interface choose_conv_method
        module procedure choose_conv_method_real, choose_conv_method_complex
        module procedure choose_conv_method_real_complex, choose_conv_method_complex_real
    end interface choose_conv_method

    public :: convolve, correlate, correlation_lags, choose_conv_method

contains

    pure subroutine convolution_slice(n, m, mode, first, length, ok)
        integer, intent(in) :: n !! length of first signal, > 0
        integer, intent(in) :: m !! length of second signal, > 0
        character(len=*), intent(in) :: mode !! full, same, or valid
        integer, intent(out) :: first !! one-based offset in full convolution
        integer, intent(out) :: length !! number of returned samples
        logical, intent(out) :: ok !! true for valid inputs and mode

        ok = n > 0 .and. m > 0
        first = 1
        length = 0
        if (.not. ok) return
        select case (mode)
        case ('full')
            length = n + m - 1
        case ('same')
            first = 1 + (m - 1) / 2
            length = n
        case ('valid')
            first = min(n, m)
            length = abs(n - m) + 1
        case default
            ok = .false.
        end select
    end subroutine convolution_slice

    pure subroutine convolution_mode(method, ok)
        character(len=*), intent(in) :: method !! 'auto' or 'direct'
        logical, intent(out) :: ok !! true for supported method
        ok = method == 'direct' .or. method == 'auto' .or. method == 'fft'
    end subroutine convolution_mode

    recursive function convolve_real(in1, in2, mode, method, status) result(y)
        real(dp), intent(in) :: in1(:) !! first real signal
        real(dp), intent(in) :: in2(:) !! second real signal
        character(len=*), intent(in), optional :: mode !! 'full' (default), 'same', 'valid'
        character(len=*), intent(in), optional :: method !! 'auto', 'direct' or 'fft', both use direct
        integer, intent(out), optional :: status !! zero on success, one on invalid input
        real(dp), allocatable :: y(:)

        integer :: i, j, k, first, length, n, m, fft_status
        logical :: ok
        character(len=:), allocatable :: selected, algorithm
        real(dp), allocatable :: full(:)

        if (present(status)) status = convolution_status_invalid_input
        n = size(in1)
        m = size(in2)
        selected = 'full'
        if (present(mode)) selected = mode
        call convolution_slice(n, m, selected, first, length, ok)
        if (present(method)) then
            block
                logical :: valid_method
                call convolution_mode(method, valid_method)
                ok = ok .and. valid_method
            end block
        end if
        allocate(y(0))
        if (.not. ok) return
        algorithm = 'auto'
        if (present(method)) algorithm = trim(method)
        if (algorithm == 'auto') algorithm = choose_conv_method(in1,in2,selected)
        if (algorithm == 'fft') then
            y = fftconvolve(in1,in2,selected,status=fft_status)
            if (present(status)) status = fft_status
            return
        end if
        allocate(full(n + m - 1))
        full = 0.0_dp
        do j = 1, m
            do i = 1, n
                k = i + j - 1
                full(k) = full(k) + in1(i) * in2(j)
            end do
        end do
        y = full(first:first + length - 1)
        if (present(status)) status = convolution_status_success
    end function convolve_real

    recursive function convolve_complex(in1, in2, mode, method, status) result(y)
        complex(dp), intent(in) :: in1(:) !! first complex signal
        complex(dp), intent(in) :: in2(:) !! second complex signal
        character(len=*), intent(in), optional :: mode !! 'full' (default), 'same', 'valid'
        character(len=*), intent(in), optional :: method !! 'auto', 'direct' or 'fft', both use direct
        integer, intent(out), optional :: status !! zero on success, one on invalid input
        complex(dp), allocatable :: y(:)

        integer :: i, j, k, first, length, n, m, fft_status
        logical :: ok
        character(len=:), allocatable :: selected, algorithm
        complex(dp), allocatable :: full(:)

        if (present(status)) status = convolution_status_invalid_input
        n = size(in1)
        m = size(in2)
        selected = 'full'
        if (present(mode)) selected = mode
        call convolution_slice(n, m, selected, first, length, ok)
        if (present(method)) then
            block
                logical :: valid_method
                call convolution_mode(method, valid_method)
                ok = ok .and. valid_method
            end block
        end if
        allocate(y(0))
        if (.not. ok) return
        algorithm = 'auto'
        if (present(method)) algorithm = trim(method)
        if (algorithm == 'auto') algorithm = choose_conv_method(in1,in2,selected)
        if (algorithm == 'fft') then
            y = fftconvolve(in1,in2,selected,status=fft_status)
            if (present(status)) status = fft_status
            return
        end if
        allocate(full(n + m - 1))
        full = (0.0_dp, 0.0_dp)
        do j = 1, m
            do i = 1, n
                k = i + j - 1
                full(k) = full(k) + in1(i) * in2(j)
            end do
        end do
        y = full(first:first + length - 1)
        if (present(status)) status = convolution_status_success
    end function convolve_complex

    ! Mixed-type convolution promotes only the real operand. This retains
    ! the complex FFT/direct kernels and their existing mode/error behavior.
    function convolve_real_complex(in1, in2, mode, method, status) result(y)
        real(dp), intent(in) :: in1(:) !! first real-valued signal
        complex(dp), intent(in) :: in2(:) !! second complex-valued signal
        character(*), intent(in), optional :: mode !! full (default), same, or valid
        character(*), intent(in), optional :: method !! auto (default), direct, or fft
        integer, intent(out), optional :: status !! success or invalid input
        complex(dp), allocatable :: y(:)

        y = convolve_complex(cmplx(in1, 0.0_dp, dp), in2, mode, method, status)
    end function convolve_real_complex

    function convolve_complex_real(in1, in2, mode, method, status) result(y)
        complex(dp), intent(in) :: in1(:) !! first complex-valued signal
        real(dp), intent(in) :: in2(:) !! second real-valued signal
        character(*), intent(in), optional :: mode !! full (default), same, or valid
        character(*), intent(in), optional :: method !! auto (default), direct, or fft
        integer, intent(out), optional :: status !! success or invalid input
        complex(dp), allocatable :: y(:)

        y = convolve_complex(in1, cmplx(in2, 0.0_dp, dp), mode, method, status)
    end function convolve_complex_real

    function correlate_real(in1, in2, mode, method, status) result(y)
        real(dp), intent(in) :: in1(:) !! first real signal
        real(dp), intent(in) :: in2(:) !! second real signal
        character(len=*), intent(in), optional :: mode !! full, same, valid
        character(len=*), intent(in), optional :: method !! 'auto', 'direct' or 'fft'
        integer, intent(out), optional :: status !! zero on success
        real(dp), allocatable :: y(:)
        y = convolve_real(in1, in2(size(in2):1:-1), mode, method, status)
    end function correlate_real

    function correlate_complex(in1, in2, mode, method, status) result(y)
        complex(dp), intent(in) :: in1(:) !! first complex signal
        complex(dp), intent(in) :: in2(:) !! second complex signal
        character(len=*), intent(in), optional :: mode !! full, same, valid
        character(len=*), intent(in), optional :: method !! 'auto', 'direct' or 'fft'
        integer, intent(out), optional :: status !! zero on success
        complex(dp), allocatable :: y(:)
        y = convolve_complex(in1, conjg(in2(size(in2):1:-1)), mode, method, status)
    end function correlate_complex

    function correlate_real_complex(in1, in2, mode, method, status) result(y)
        real(dp), intent(in) :: in1(:) !! first real-valued signal
        complex(dp), intent(in) :: in2(:) !! second complex-valued signal, conjugated
        character(*), intent(in), optional :: mode !! full (default), same, or valid
        character(*), intent(in), optional :: method !! auto (default), direct, or fft
        integer, intent(out), optional :: status !! success or invalid input
        complex(dp), allocatable :: y(:)

        y = correlate_complex(cmplx(in1, 0.0_dp, dp), in2, mode, method, status)
    end function correlate_real_complex

    function correlate_complex_real(in1, in2, mode, method, status) result(y)
        complex(dp), intent(in) :: in1(:) !! first complex-valued signal
        real(dp), intent(in) :: in2(:) !! second real-valued signal, conjugated
        character(*), intent(in), optional :: mode !! full (default), same, or valid
        character(*), intent(in), optional :: method !! auto (default), direct, or fft
        integer, intent(out), optional :: status !! success or invalid input
        complex(dp), allocatable :: y(:)

        y = correlate_complex(in1, cmplx(in2, 0.0_dp, dp), mode, method, status)
    end function correlate_complex_real

    function correlation_lags(in1_len, in2_len, mode, status) result(lags)
        integer, intent(in) :: in1_len !! length of first signal, > 0
        integer, intent(in) :: in2_len !! length of second signal, > 0
        character(len=*), intent(in), optional :: mode !! full, same, or valid
        integer, intent(out), optional :: status !! zero on success
        integer, allocatable :: lags(:)

        integer :: first, length, k
        logical :: ok
        character(len=:), allocatable :: selected

        if (present(status)) status = convolution_status_invalid_input
        selected = 'full'
        if (present(mode)) selected = mode
        call convolution_slice(in1_len, in2_len, selected, first, length, ok)
        allocate(lags(0))
        if (.not. ok) return
        if (selected == 'same') then
            ! SciPy correlation_lags centers the full lag range; its even-
            ! kernel convention may be shifted from correlate(..., 'same').
            first = 1 + (in1_len + in2_len - 1) / 2 - in1_len / 2
        end if
        lags = [(first + k - 1 - in2_len, k = 1, length)]
        if (present(status)) status = convolution_status_success
    end function correlation_lags


    ! Cost estimates are appropriate for SciFort's direct O(n*m) loop and
    ! radix-2 FFT, not SciPy's platform-specific autotuning coefficients.
    pure function estimated_method(n,m) result(method)
        integer, intent(in) :: n !! length of first input
        integer, intent(in) :: m !! length of second input
        character(len=6) :: method
        integer :: nf, lg

        method = 'direct'
        if (min(n,m) < 16 .or. n > huge(n)-m+1) return
        nf = 1
        lg = 0
        do while (nf < n+m-1)
            if (nf > huge(nf)/2) return
            nf = nf*2
            lg = lg+1
        end do
        if (real(n,dp)*real(m,dp) > 8.0_dp*real(nf,dp)*real(max(1,lg),dp)) method = 'fft'
    end function estimated_method

    function choose_conv_method_real(in1,in2,mode,measure,times,status) result(method)
        real(dp), intent(in) :: in1(:) !! first real input
        real(dp), intent(in) :: in2(:) !! second real input
        character(*), intent(in), optional :: mode !! full, same, valid
        logical, intent(in), optional :: measure !! true benchmarks direct versus FFT
        real(dp), intent(out), optional :: times(:) !! measured [direct, fft] CPU seconds
        integer, intent(out), optional :: status !! zero on valid input
        character(len=6) :: method
        character(:), allocatable :: selected
        integer :: first, count
        logical :: ok, do_measure
        real(dp) :: t0,t1,t2
        real(dp), allocatable :: direct_result(:), fft_result(:)

        method = 'direct'
        if (present(status)) status = convolution_status_invalid_input
        if (present(times)) times = -1.0_dp
        selected = 'full'
        if (present(mode)) selected = trim(mode)
        call convolution_slice(size(in1),size(in2),selected,first,count,ok)
        if (.not. ok) return
        if (present(status)) status = convolution_status_success
        if (.not. all(ieee_is_finite(in1)) .or. .not. all(ieee_is_finite(in2))) return
        method = estimated_method(size(in1),size(in2))
        do_measure = .false.
        if (present(measure)) do_measure = measure
        if (.not. do_measure) return
        call cpu_time(t0)
        direct_result = convolve_real(in1,in2,selected,'direct')
        call cpu_time(t1)
        fft_result = fftconvolve(in1,in2,selected)
        call cpu_time(t2)
        method = 'direct'
        if (t2-t1 < t1-t0) method = 'fft'
        if (present(times)) then
            if (size(times)>0) times(1) = t1-t0
            if (size(times)>1) times(2) = t2-t1
        end if
    end function choose_conv_method_real

    function choose_conv_method_complex(in1,in2,mode,measure,times,status) result(method)
        complex(dp), intent(in) :: in1(:) !! first complex input
        complex(dp), intent(in) :: in2(:) !! second complex input
        character(*), intent(in), optional :: mode !! full, same, valid
        logical, intent(in), optional :: measure !! benchmark if true
        real(dp), intent(out), optional :: times(:) !! measured [direct, fft] CPU seconds
        integer, intent(out), optional :: status !! zero on valid input
        character(len=6) :: method
        character(:), allocatable :: selected
        integer :: first,count
        logical :: ok,do_measure
        real(dp) :: t0,t1,t2
        complex(dp), allocatable :: direct_result(:),fft_result(:)

        method = 'direct'
        if (present(status)) status = convolution_status_invalid_input
        if (present(times)) times = -1.0_dp
        selected = 'full'
        if (present(mode)) selected = trim(mode)
        call convolution_slice(size(in1),size(in2),selected,first,count,ok)
        if (.not. ok) return
        if (present(status)) status = convolution_status_success
        if (.not. all(ieee_is_finite(real(in1,dp))) .or. .not. all(ieee_is_finite(aimag(in1))) .or. &
            .not. all(ieee_is_finite(real(in2,dp))) .or. .not. all(ieee_is_finite(aimag(in2)))) return
        method = estimated_method(size(in1),size(in2))
        do_measure = .false.
        if (present(measure)) do_measure = measure
        if (.not. do_measure) return
        call cpu_time(t0)
        direct_result = convolve_complex(in1,in2,selected,'direct')
        call cpu_time(t1)
        fft_result = fftconvolve(in1,in2,selected)
        call cpu_time(t2)
        method = 'direct'
        if (t2-t1 < t1-t0) method = 'fft'
        if (present(times)) then
            if (size(times)>0) times(1) = t1-t0
            if (size(times)>1) times(2) = t2-t1
        end if
    end function choose_conv_method_complex

    ! Match mixed-array promotion for convolution algorithm selection.
    ! Both benchmarking and the returned timing pair use the complex path.
    function choose_conv_method_real_complex(in1,in2,mode,measure,times,status) result(method)
        real(dp), intent(in) :: in1(:) !! first real-valued signal
        complex(dp), intent(in) :: in2(:) !! second complex-valued signal
        character(*), intent(in), optional :: mode !! full (default), same, or valid
        logical, intent(in), optional :: measure !! benchmark direct versus FFT
        real(dp), intent(out), optional :: times(:) !! [direct, FFT] CPU seconds
        integer, intent(out), optional :: status !! success or invalid input
        character(len=6) :: method

        method = choose_conv_method_complex(cmplx(in1, 0.0_dp, dp), in2, &
            mode, measure, times, status)
    end function choose_conv_method_real_complex

    function choose_conv_method_complex_real(in1,in2,mode,measure,times,status) result(method)
        complex(dp), intent(in) :: in1(:) !! first complex-valued signal
        real(dp), intent(in) :: in2(:) !! second real-valued signal
        character(*), intent(in), optional :: mode !! full (default), same, or valid
        logical, intent(in), optional :: measure !! benchmark direct versus FFT
        real(dp), intent(out), optional :: times(:) !! [direct, FFT] CPU seconds
        integer, intent(out), optional :: status !! success or invalid input
        character(len=6) :: method

        method = choose_conv_method_complex(in1, cmplx(in2, 0.0_dp, dp), &
            mode, measure, times, status)
    end function choose_conv_method_complex_real

end module scifort_signal_convolve
