! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Rank-1..15 spectral estimators using the tested 1-D SciFort kernels.
! An assumed-rank interface and a single lane walker allow any signal axis.
! Result data are committed only when all channels have been computed.
module scifort_signal_spectral_nd
    use scifort_kinds, only : dp
    use scifort_signal_spectral, only : spectral_result, cross_spectral_result, &
        spectral_status_success, spectral_status_invalid_input, &
        periodogram_real => periodogram, welch_real => welch, &
        csd_real => csd, coherence_real => coherence, &
        spectral_real_detrend_callback
    use scifort_signal_spectral_complex, only : periodogram_complex, &
        welch_complex, csd_complex, coherence_complex, spectral_complex_detrend_callback
    implicit none
    private

    public :: periodogram_into, welch_into, csd_into, coherence_into
    ! Shared rank-independent staging utilities for time-frequency transforms.
    public :: flatten_real, flatten_complex, scatter_real, scatter_complex

    type :: spectral_job
        integer :: method = 0 ! 1=periodogram, 2=welch, 3=csd, 4=coherence
        logical :: first_real = .true.
        logical :: second_real = .true.
        real(dp), allocatable :: fs
        character(:), allocatable :: window, detrend, scaling, average
        real(dp), allocatable :: window_parameter, window_parameter2
        real(dp), allocatable :: window_coefficients(:), window_values(:)
        integer, allocatable :: nperseg, noverlap, nfft
        logical, allocatable :: return_onesided !! optional, default depends on input type
        procedure(spectral_real_detrend_callback), pointer, nopass :: real_callback => null()
        procedure(spectral_complex_detrend_callback), pointer, nopass :: complex_callback => null()
    end type spectral_job

    interface periodogram_into
        module procedure periodogram_into_real, periodogram_into_complex
    end interface
    interface welch_into
        module procedure welch_into_real, welch_into_complex
    end interface
    interface csd_into
        module procedure csd_into_rr, csd_into_cc, csd_into_rc, csd_into_cr
    end interface
    interface coherence_into
        module procedure coherence_into_rr, coherence_into_cc, &
            coherence_into_rc, coherence_into_cr
    end interface

contains

    subroutine check_shapes(d1, d2, dout, axis, a, status)
        integer, intent(in) :: d1(:) !! first input shape
        integer, intent(in) :: d2(:) !! second input shape
        integer, intent(in) :: dout(:) !! result shape
        integer, intent(in), optional :: axis !! one-based or negative, default last
        integer, intent(out) :: a !! resolved positive axis
        integer, intent(out) :: status !! validation status
        integer :: j, ndim
        ndim = size(d1)
        status = spectral_status_invalid_input
        a = ndim
        if (ndim < 1 .or. ndim > 15) return
        if (size(d2) /= ndim .or. size(dout) /= ndim) return
        if (present(axis)) a = axis
        if (a < 0) a = a + ndim + 1
        if (a < 1 .or. a > ndim) return
        if (any(d1 <= 0) .or. any(d2 <= 0) .or. any(dout <= 0)) return
        do j = 1, ndim
            if (j == a) cycle
            if (d1(j) /= d2(j) .or. d1(j) /= dout(j)) return
        end do
        status = spectral_status_success
    end subroutine check_shapes

    function process_one(u, v, job, frequencies, status) result(power)
        complex(dp), intent(in) :: u(:) !! first channel samples
        complex(dp), intent(in) :: v(:) !! second channel samples, if used
        type(spectral_job), intent(in) :: job !! estimator arguments
        real(dp), allocatable, intent(out) :: frequencies(:) !! spectral grid
        integer, intent(out) :: status !! one-dimensional estimator status
        complex(dp), allocatable :: power(:)
        type(spectral_result) :: r
        type(cross_spectral_result) :: c
        allocate(power(0), frequencies(0))
        status = spectral_status_invalid_input
        select case (job%method)
        case (1)
            if (job%first_real) then
                if (associated(job%real_callback)) then
                    r = periodogram_real(real(u,dp),job%fs,job%window,job%nfft, &
                        job%detrend,job%scaling, &
                        window_parameter=job%window_parameter, &
                        window_parameter2=job%window_parameter2, &
                        window_coefficients=job%window_coefficients, &
                        window_values=job%window_values, &
                        return_onesided=job%return_onesided, &
                        detrend_callback=job%real_callback)
                else
                    r = periodogram_real(real(u,dp),job%fs,job%window,job%nfft, &
                        job%detrend,job%scaling, &
                        window_parameter=job%window_parameter, &
                        window_parameter2=job%window_parameter2, &
                        window_coefficients=job%window_coefficients, &
                        window_values=job%window_values, &
                        return_onesided=job%return_onesided)
                end if
            else
                if (associated(job%complex_callback)) then
                    r = periodogram_complex(u,job%fs,job%window,job%nfft,job%detrend,job%scaling, &
                        window_parameter=job%window_parameter, &
                        window_parameter2=job%window_parameter2, &
                        window_coefficients=job%window_coefficients, &
                        window_values=job%window_values, &
                        return_onesided=job%return_onesided, &
                        detrend_callback=job%complex_callback)
                else
                    r = periodogram_complex(u,job%fs,job%window,job%nfft,job%detrend,job%scaling, &
                        window_parameter=job%window_parameter, &
                        window_parameter2=job%window_parameter2, &
                        window_coefficients=job%window_coefficients, &
                        window_values=job%window_values, &
                        return_onesided=job%return_onesided)
                end if
            end if
        case (2)
            if (job%first_real) then
                if (associated(job%real_callback)) then
                    r = welch_real(real(u,dp),job%fs,job%window,job%nperseg, &
                        job%noverlap,job%nfft,job%detrend,job%scaling,job%average, &
                        window_parameter=job%window_parameter, &
                        window_parameter2=job%window_parameter2, &
                        window_coefficients=job%window_coefficients, &
                        window_values=job%window_values, &
                        return_onesided=job%return_onesided, &
                        detrend_callback=job%real_callback)
                else
                    r = welch_real(real(u,dp),job%fs,job%window,job%nperseg, &
                        job%noverlap,job%nfft,job%detrend,job%scaling,job%average, &
                        window_parameter=job%window_parameter, &
                        window_parameter2=job%window_parameter2, &
                        window_coefficients=job%window_coefficients, &
                        window_values=job%window_values, &
                        return_onesided=job%return_onesided)
                end if
            else
                if (associated(job%complex_callback)) then
                    r = welch_complex(u,job%fs,job%window,job%nperseg, &
                        job%noverlap,job%nfft,job%detrend,job%scaling,job%average, &
                        window_parameter=job%window_parameter, &
                        window_parameter2=job%window_parameter2, &
                        window_coefficients=job%window_coefficients, &
                        window_values=job%window_values, &
                        return_onesided=job%return_onesided, &
                        detrend_callback=job%complex_callback)
                else
                    r = welch_complex(u,job%fs,job%window,job%nperseg, &
                        job%noverlap,job%nfft,job%detrend,job%scaling,job%average, &
                        window_parameter=job%window_parameter, &
                        window_parameter2=job%window_parameter2, &
                        window_coefficients=job%window_coefficients, &
                        window_values=job%window_values, &
                        return_onesided=job%return_onesided)
                end if
            end if
        case (3)
            if (job%first_real .and. job%second_real) then
                if (associated(job%real_callback)) then
                    c = csd_real(real(u,dp),real(v,dp),job%fs,job%window,job%nperseg, &
                        job%noverlap,job%nfft,job%detrend,job%scaling,job%average, &
                        window_parameter=job%window_parameter, &
                        window_parameter2=job%window_parameter2, &
                        window_coefficients=job%window_coefficients, &
                        window_values=job%window_values, &
                        return_onesided=job%return_onesided, &
                        detrend_callback=job%real_callback)
                else
                    c = csd_real(real(u,dp),real(v,dp),job%fs,job%window,job%nperseg, &
                        job%noverlap,job%nfft,job%detrend,job%scaling,job%average, &
                        window_parameter=job%window_parameter, &
                        window_parameter2=job%window_parameter2, &
                        window_coefficients=job%window_coefficients, &
                        window_values=job%window_values, &
                        return_onesided=job%return_onesided)
                end if
            else
                if (associated(job%complex_callback)) then
                    c = csd_complex(u,v,job%fs,job%window,job%nperseg, &
                        job%noverlap,job%nfft,job%detrend,job%scaling,job%average, &
                        window_parameter=job%window_parameter, &
                        window_parameter2=job%window_parameter2, &
                        window_coefficients=job%window_coefficients, &
                        window_values=job%window_values, &
                        return_onesided=job%return_onesided, &
                        detrend_callback=job%complex_callback)
                else
                    c = csd_complex(u,v,job%fs,job%window,job%nperseg, &
                        job%noverlap,job%nfft,job%detrend,job%scaling,job%average, &
                        window_parameter=job%window_parameter, &
                        window_parameter2=job%window_parameter2, &
                        window_coefficients=job%window_coefficients, &
                        window_values=job%window_values, &
                        return_onesided=job%return_onesided)
                end if
            end if
        case (4)
            if (job%first_real .and. job%second_real) then
                if (associated(job%real_callback)) then
                    r = coherence_real(real(u,dp),real(v,dp),job%fs,job%window, &
                        job%nperseg,job%noverlap,job%nfft,job%detrend, &
                        window_parameter=job%window_parameter, &
                        window_parameter2=job%window_parameter2, &
                        window_coefficients=job%window_coefficients, &
                        window_values=job%window_values, &
                        detrend_callback=job%real_callback)
                else
                    r = coherence_real(real(u,dp),real(v,dp),job%fs,job%window, &
                        job%nperseg,job%noverlap,job%nfft,job%detrend, &
                        window_parameter=job%window_parameter, &
                        window_parameter2=job%window_parameter2, &
                        window_coefficients=job%window_coefficients, &
                        window_values=job%window_values)
                end if
            else
                if (associated(job%complex_callback)) then
                    r = coherence_complex(u,v,job%fs,job%window,job%nperseg, &
                        job%noverlap,job%nfft,job%detrend, &
                        window_parameter=job%window_parameter, &
                        window_parameter2=job%window_parameter2, &
                        window_coefficients=job%window_coefficients, &
                        window_values=job%window_values, &
                        detrend_callback=job%complex_callback)
                else
                    r = coherence_complex(u,v,job%fs,job%window,job%nperseg, &
                        job%noverlap,job%nfft,job%detrend, &
                        window_parameter=job%window_parameter, &
                        window_parameter2=job%window_parameter2, &
                        window_coefficients=job%window_coefficients, &
                        window_values=job%window_values)
                end if
            end if
        case default
            return
        end select
        if (job%method == 3) then
            status = c%status
            if (status /= spectral_status_success) return
            power = c%power
            frequencies = c%frequencies
        else
            status = r%status
            if (status /= spectral_status_success) return
            power = cmplx(r%power,0.0_dp,dp)
            frequencies = r%frequencies
        end if
    end function process_one

    subroutine process_lanes(first, second, d1, d2, dout, a, job, result, frequencies, status)
        complex(dp), intent(in) :: first(:) !! flat first input array
        complex(dp), intent(in) :: second(:) !! flat second input array
        integer, intent(in) :: d1(:) !! first input extents
        integer, intent(in) :: d2(:) !! second input extents
        integer, intent(in) :: dout(:) !! output extents
        integer, intent(in) :: a !! selected positive axis
        type(spectral_job), intent(in) :: job !! parameters
        complex(dp), intent(out) :: result(:) !! flat complete output staging
        real(dp), allocatable, intent(out) :: frequencies(:) !! returned grid
        integer, intent(out) :: status !! success only if all channels pass
        complex(dp), allocatable :: x(:), y(:), p(:)
        real(dp), allocatable :: grid(:)
        integer :: inner, outer, n1, n2, nout, i, j, k, st, index
        status = spectral_status_invalid_input
        allocate(frequencies(0))
        inner = product(d1(:a-1))
        outer = product(d1(a+1:))
        n1 = d1(a)
        n2 = d2(a)
        nout = dout(a)
        if (inner <= 0 .or. outer <= 0 .or. n1 <= 0 .or. n2 <= 0 .or. nout <= 0) return
        allocate(x(n1),y(n2))
        do k = 0, outer-1
            do i = 0, inner-1
                do j = 0, n1-1
                    index = k*inner*n1+j*inner+i+1
                    x(j+1) = first(index)
                end do
                if (job%method >= 3) then
                    do j = 0, n2-1
                        index = k*inner*n2+j*inner+i+1
                        y(j+1) = second(index)
                    end do
                end if
                p = process_one(x,y,job,grid,st)
                if (st /= spectral_status_success .or. size(p) /= nout) return
                if (size(grid) /= nout) return
                if (k == 0 .and. i == 0) then
                    frequencies = grid
                else
                    if (any(grid /= frequencies)) return
                end if
                do j = 0, nout-1
                    index = k*inner*nout+j*inner+i+1
                    result(index) = p(j+1)
                end do
            end do
        end do
        status = spectral_status_success
    end subroutine process_lanes

    subroutine flatten_real(x,flat)
        real(dp), intent(in) :: x(..) !! rank-1..15 input
        complex(dp), allocatable, intent(out) :: flat(:) !! column-major flat values
        select rank (x)
        ! BEGIN GENERATED FLATTEN REAL
        rank (1)
            flat = cmplx(reshape(x,[size(x)]),0.0_dp,dp)
        rank (2)
            flat = cmplx(reshape(x,[size(x)]),0.0_dp,dp)
        rank (3)
            flat = cmplx(reshape(x,[size(x)]),0.0_dp,dp)
        rank (4)
            flat = cmplx(reshape(x,[size(x)]),0.0_dp,dp)
        rank (5)
            flat = cmplx(reshape(x,[size(x)]),0.0_dp,dp)
        rank (6)
            flat = cmplx(reshape(x,[size(x)]),0.0_dp,dp)
        rank (7)
            flat = cmplx(reshape(x,[size(x)]),0.0_dp,dp)
        rank (8)
            flat = cmplx(reshape(x,[size(x)]),0.0_dp,dp)
        rank (9)
            flat = cmplx(reshape(x,[size(x)]),0.0_dp,dp)
        rank (10)
            flat = cmplx(reshape(x,[size(x)]),0.0_dp,dp)
        rank (11)
            flat = cmplx(reshape(x,[size(x)]),0.0_dp,dp)
        rank (12)
            flat = cmplx(reshape(x,[size(x)]),0.0_dp,dp)
        rank (13)
            flat = cmplx(reshape(x,[size(x)]),0.0_dp,dp)
        rank (14)
            flat = cmplx(reshape(x,[size(x)]),0.0_dp,dp)
        rank (15)
            flat = cmplx(reshape(x,[size(x)]),0.0_dp,dp)
        ! END GENERATED FLATTEN REAL
        rank default
            allocate(flat(0))
        end select
    end subroutine flatten_real

    subroutine flatten_complex(x,flat)
        complex(dp), intent(in) :: x(..) !! rank-1..15 input
        complex(dp), allocatable, intent(out) :: flat(:) !! column-major flat values
        select rank (x)
        ! BEGIN GENERATED FLATTEN COMPLEX
        rank (1)
            flat = reshape(x,[size(x)])
        rank (2)
            flat = reshape(x,[size(x)])
        rank (3)
            flat = reshape(x,[size(x)])
        rank (4)
            flat = reshape(x,[size(x)])
        rank (5)
            flat = reshape(x,[size(x)])
        rank (6)
            flat = reshape(x,[size(x)])
        rank (7)
            flat = reshape(x,[size(x)])
        rank (8)
            flat = reshape(x,[size(x)])
        rank (9)
            flat = reshape(x,[size(x)])
        rank (10)
            flat = reshape(x,[size(x)])
        rank (11)
            flat = reshape(x,[size(x)])
        rank (12)
            flat = reshape(x,[size(x)])
        rank (13)
            flat = reshape(x,[size(x)])
        rank (14)
            flat = reshape(x,[size(x)])
        rank (15)
            flat = reshape(x,[size(x)])
        ! END GENERATED FLATTEN COMPLEX
        rank default
            allocate(flat(0))
        end select
    end subroutine flatten_complex

    subroutine scatter_real(y,flat)
        real(dp), intent(inout) :: y(..) !! preallocated result, rank 1..15
        complex(dp), intent(in) :: flat(:) !! staged result values
        select rank (y)
        ! BEGIN GENERATED SCATTER REAL
        rank (1)
            y = reshape(real(flat,dp),shape(y))
        rank (2)
            y = reshape(real(flat,dp),shape(y))
        rank (3)
            y = reshape(real(flat,dp),shape(y))
        rank (4)
            y = reshape(real(flat,dp),shape(y))
        rank (5)
            y = reshape(real(flat,dp),shape(y))
        rank (6)
            y = reshape(real(flat,dp),shape(y))
        rank (7)
            y = reshape(real(flat,dp),shape(y))
        rank (8)
            y = reshape(real(flat,dp),shape(y))
        rank (9)
            y = reshape(real(flat,dp),shape(y))
        rank (10)
            y = reshape(real(flat,dp),shape(y))
        rank (11)
            y = reshape(real(flat,dp),shape(y))
        rank (12)
            y = reshape(real(flat,dp),shape(y))
        rank (13)
            y = reshape(real(flat,dp),shape(y))
        rank (14)
            y = reshape(real(flat,dp),shape(y))
        rank (15)
            y = reshape(real(flat,dp),shape(y))
        ! END GENERATED SCATTER REAL
        rank default
        end select
    end subroutine scatter_real

    subroutine scatter_complex(y,flat)
        complex(dp), intent(inout) :: y(..) !! preallocated result, rank 1..15
        complex(dp), intent(in) :: flat(:) !! staged result values
        select rank (y)
        ! BEGIN GENERATED SCATTER COMPLEX
        rank (1)
            y = reshape(flat,shape(y))
        rank (2)
            y = reshape(flat,shape(y))
        rank (3)
            y = reshape(flat,shape(y))
        rank (4)
            y = reshape(flat,shape(y))
        rank (5)
            y = reshape(flat,shape(y))
        rank (6)
            y = reshape(flat,shape(y))
        rank (7)
            y = reshape(flat,shape(y))
        rank (8)
            y = reshape(flat,shape(y))
        rank (9)
            y = reshape(flat,shape(y))
        rank (10)
            y = reshape(flat,shape(y))
        rank (11)
            y = reshape(flat,shape(y))
        rank (12)
            y = reshape(flat,shape(y))
        rank (13)
            y = reshape(flat,shape(y))
        rank (14)
            y = reshape(flat,shape(y))
        rank (15)
            y = reshape(flat,shape(y))
        ! END GENERATED SCATTER COMPLEX
        rank default
        end select
    end subroutine scatter_complex

    subroutine apply_to_real(d1,d2,u,v,result,frequencies,job,axis,status)
        integer, intent(in) :: d1(:) !! first input shape
        integer, intent(in) :: d2(:) !! second input shape
        complex(dp), intent(in) :: u(:) !! flat first input
        complex(dp), intent(in) :: v(:) !! flat second input or empty
        real(dp), intent(inout) :: result(..) !! caller-provided rank-1..15 result
        real(dp), allocatable, intent(out) :: frequencies(:) !! frequencies of output
        type(spectral_job), intent(in) :: job !! fixed estimator configuration
        integer, intent(in), optional :: axis !! one-based/negative; default last
        integer, intent(out) :: status !! zero on success
        integer :: a
        complex(dp), allocatable :: staged(:)
        call check_shapes(d1,d2,shape(result),axis,a,status)
        if (status /= spectral_status_success) then
            allocate(frequencies(0))
            return
        end if
        allocate(staged(size(result)))
        call process_lanes(u,v,d1,d2,shape(result),a,job,staged,frequencies,status)
        if (status == spectral_status_success) then
            call scatter_real(result,staged)
        else
            if (allocated(frequencies)) deallocate(frequencies)
            allocate(frequencies(0))
        end if
    end subroutine apply_to_real

    subroutine apply_to_complex(d1,d2,u,v,result,frequencies,job,axis,status)
        integer, intent(in) :: d1(:) !! first input shape
        integer, intent(in) :: d2(:) !! second input shape
        complex(dp), intent(in) :: u(:) !! flat first input
        complex(dp), intent(in) :: v(:) !! flat second input or empty
        complex(dp), intent(inout) :: result(..) !! caller-provided rank-1..15 result
        real(dp), allocatable, intent(out) :: frequencies(:) !! frequencies of output
        type(spectral_job), intent(in) :: job !! fixed estimator configuration
        integer, intent(in), optional :: axis !! one-based/negative; default last
        integer, intent(out) :: status !! zero on success
        integer :: a
        complex(dp), allocatable :: staged(:)
        call check_shapes(d1,d2,shape(result),axis,a,status)
        if (status /= spectral_status_success) then
            allocate(frequencies(0))
            return
        end if
        allocate(staged(size(result)))
        call process_lanes(u,v,d1,d2,shape(result),a,job,staged,frequencies,status)
        if (status == spectral_status_success) then
            call scatter_complex(result,staged)
        else
            if (allocated(frequencies)) deallocate(frequencies)
            allocate(frequencies(0))
        end if
    end subroutine apply_to_complex

    ! BEGIN GENERATED PUBLIC WRAPPERS
    subroutine periodogram_into_real(x,result,frequencies,fs,window,nfft,detrend,scaling,axis,status, &
        window_parameter,window_parameter2,window_coefficients,window_values,return_onesided, detrend_callback)
        real(dp), intent(in) :: x(..) !! first input, rank 1-15
        real(dp), intent(inout) :: result(..) !! caller-allocated spectrum
        real(dp), allocatable, intent(out) :: frequencies(:) !! sample frequencies
        real(dp), intent(in), optional :: fs !! SciPy spectral option
        character(*), intent(in), optional :: window !! SciPy spectral option
        integer, intent(in), optional :: nfft !! SciPy spectral option
        character(*), intent(in), optional :: detrend !! SciPy spectral option
        character(*), intent(in), optional :: scaling !! SciPy spectral option
        integer, intent(in), optional :: axis !! one-based or negative; default last
        integer, intent(out), optional :: status !! zero on success
        real(dp), intent(in), optional :: window_parameter !! named window shape parameter
        real(dp), intent(in), optional :: window_parameter2 !! second named window parameter
        real(dp), intent(in), optional :: window_coefficients(:) !! general cosine coefficients
        real(dp), intent(in), optional :: window_values(:) !! explicit real window samples
        logical, intent(in), optional :: return_onesided !! false retains FFT negative frequencies
        procedure(spectral_real_detrend_callback), optional :: detrend_callback !! segment callback; excludes detrend=
        type(spectral_job) :: job
        complex(dp), allocatable :: u(:), v(:)
        integer :: st
        job%method = 1
        job%first_real = .true.
        job%second_real = .true.
        if (present(fs)) job%fs = fs
        if (present(window)) job%window = window
        if (present(nfft)) job%nfft = nfft
        if (present(detrend)) job%detrend = detrend
        if (present(detrend_callback)) job%real_callback => detrend_callback
        if (present(scaling)) job%scaling = scaling
        if (present(window_parameter)) job%window_parameter = window_parameter
        if (present(window_parameter2)) job%window_parameter2 = window_parameter2
        if (present(window_coefficients)) job%window_coefficients = window_coefficients
        if (present(window_values)) job%window_values = window_values
        if (present(return_onesided)) job%return_onesided = return_onesided
        call flatten_real(x,u)
        allocate(v(0))
        call apply_to_real(shape(x),shape(x),u,v,result,frequencies,job,axis,st)
        if (present(status)) status = st
    end subroutine periodogram_into_real

    subroutine periodogram_into_complex(x,result,frequencies,fs,window,nfft,detrend,scaling,axis,status, &
        window_parameter,window_parameter2,window_coefficients,window_values,return_onesided, detrend_callback)
        complex(dp), intent(in) :: x(..) !! first input, rank 1-15
        real(dp), intent(inout) :: result(..) !! caller-allocated spectrum
        real(dp), allocatable, intent(out) :: frequencies(:) !! sample frequencies
        real(dp), intent(in), optional :: fs !! SciPy spectral option
        character(*), intent(in), optional :: window !! SciPy spectral option
        integer, intent(in), optional :: nfft !! SciPy spectral option
        character(*), intent(in), optional :: detrend !! SciPy spectral option
        character(*), intent(in), optional :: scaling !! SciPy spectral option
        integer, intent(in), optional :: axis !! one-based or negative; default last
        integer, intent(out), optional :: status !! zero on success
        real(dp), intent(in), optional :: window_parameter !! named window shape parameter
        real(dp), intent(in), optional :: window_parameter2 !! second named window parameter
        real(dp), intent(in), optional :: window_coefficients(:) !! general cosine coefficients
        real(dp), intent(in), optional :: window_values(:) !! explicit real window samples
        logical, intent(in), optional :: return_onesided !! false retains FFT negative frequencies
        procedure(spectral_complex_detrend_callback), optional :: detrend_callback !! segment callback; excludes detrend=
        type(spectral_job) :: job
        complex(dp), allocatable :: u(:), v(:)
        integer :: st
        job%method = 1
        job%first_real = .false.
        job%second_real = .true.
        if (present(fs)) job%fs = fs
        if (present(window)) job%window = window
        if (present(nfft)) job%nfft = nfft
        if (present(detrend)) job%detrend = detrend
        if (present(detrend_callback)) job%complex_callback => detrend_callback
        if (present(scaling)) job%scaling = scaling
        if (present(window_parameter)) job%window_parameter = window_parameter
        if (present(window_parameter2)) job%window_parameter2 = window_parameter2
        if (present(window_coefficients)) job%window_coefficients = window_coefficients
        if (present(window_values)) job%window_values = window_values
        if (present(return_onesided)) job%return_onesided = return_onesided
        call flatten_complex(x,u)
        allocate(v(0))
        call apply_to_real(shape(x),shape(x),u,v,result,frequencies,job,axis,st)
        if (present(status)) status = st
    end subroutine periodogram_into_complex

    subroutine welch_into_real(x,result,frequencies,fs,window,nperseg,noverlap,nfft,detrend,scaling,average,axis,status, &
        window_parameter,window_parameter2,window_coefficients,window_values,return_onesided, detrend_callback)
        real(dp), intent(in) :: x(..) !! first input, rank 1-15
        real(dp), intent(inout) :: result(..) !! caller-allocated spectrum
        real(dp), allocatable, intent(out) :: frequencies(:) !! sample frequencies
        real(dp), intent(in), optional :: fs !! SciPy spectral option
        character(*), intent(in), optional :: window !! SciPy spectral option
        integer, intent(in), optional :: nperseg !! SciPy spectral option
        integer, intent(in), optional :: noverlap !! SciPy spectral option
        integer, intent(in), optional :: nfft !! SciPy spectral option
        character(*), intent(in), optional :: detrend !! SciPy spectral option
        character(*), intent(in), optional :: scaling !! SciPy spectral option
        character(*), intent(in), optional :: average !! SciPy spectral option
        integer, intent(in), optional :: axis !! one-based or negative; default last
        integer, intent(out), optional :: status !! zero on success
        real(dp), intent(in), optional :: window_parameter !! named window shape parameter
        real(dp), intent(in), optional :: window_parameter2 !! second named window parameter
        real(dp), intent(in), optional :: window_coefficients(:) !! general cosine coefficients
        real(dp), intent(in), optional :: window_values(:) !! explicit real window samples
        logical, intent(in), optional :: return_onesided !! false retains FFT negative frequencies
        procedure(spectral_real_detrend_callback), optional :: detrend_callback !! segment callback; excludes detrend=
        type(spectral_job) :: job
        complex(dp), allocatable :: u(:), v(:)
        integer :: st
        job%method = 2
        job%first_real = .true.
        job%second_real = .true.
        if (present(fs)) job%fs = fs
        if (present(window)) job%window = window
        if (present(nperseg)) job%nperseg = nperseg
        if (present(noverlap)) job%noverlap = noverlap
        if (present(nfft)) job%nfft = nfft
        if (present(detrend)) job%detrend = detrend
        if (present(detrend_callback)) job%real_callback => detrend_callback
        if (present(scaling)) job%scaling = scaling
        if (present(average)) job%average = average
        if (present(window_parameter)) job%window_parameter = window_parameter
        if (present(window_parameter2)) job%window_parameter2 = window_parameter2
        if (present(window_coefficients)) job%window_coefficients = window_coefficients
        if (present(window_values)) job%window_values = window_values
        if (present(return_onesided)) job%return_onesided = return_onesided
        call flatten_real(x,u)
        allocate(v(0))
        call apply_to_real(shape(x),shape(x),u,v,result,frequencies,job,axis,st)
        if (present(status)) status = st
    end subroutine welch_into_real

    subroutine welch_into_complex(x,result,frequencies,fs,window,nperseg,noverlap,nfft,detrend,scaling,average,axis,status, &
        window_parameter,window_parameter2,window_coefficients,window_values,return_onesided, detrend_callback)
        complex(dp), intent(in) :: x(..) !! first input, rank 1-15
        real(dp), intent(inout) :: result(..) !! caller-allocated spectrum
        real(dp), allocatable, intent(out) :: frequencies(:) !! sample frequencies
        real(dp), intent(in), optional :: fs !! SciPy spectral option
        character(*), intent(in), optional :: window !! SciPy spectral option
        integer, intent(in), optional :: nperseg !! SciPy spectral option
        integer, intent(in), optional :: noverlap !! SciPy spectral option
        integer, intent(in), optional :: nfft !! SciPy spectral option
        character(*), intent(in), optional :: detrend !! SciPy spectral option
        character(*), intent(in), optional :: scaling !! SciPy spectral option
        character(*), intent(in), optional :: average !! SciPy spectral option
        integer, intent(in), optional :: axis !! one-based or negative; default last
        integer, intent(out), optional :: status !! zero on success
        real(dp), intent(in), optional :: window_parameter !! named window shape parameter
        real(dp), intent(in), optional :: window_parameter2 !! second named window parameter
        real(dp), intent(in), optional :: window_coefficients(:) !! general cosine coefficients
        real(dp), intent(in), optional :: window_values(:) !! explicit real window samples
        logical, intent(in), optional :: return_onesided !! false retains FFT negative frequencies
        procedure(spectral_complex_detrend_callback), optional :: detrend_callback !! segment callback; excludes detrend=
        type(spectral_job) :: job
        complex(dp), allocatable :: u(:), v(:)
        integer :: st
        job%method = 2
        job%first_real = .false.
        job%second_real = .true.
        if (present(fs)) job%fs = fs
        if (present(window)) job%window = window
        if (present(nperseg)) job%nperseg = nperseg
        if (present(noverlap)) job%noverlap = noverlap
        if (present(nfft)) job%nfft = nfft
        if (present(detrend)) job%detrend = detrend
        if (present(detrend_callback)) job%complex_callback => detrend_callback
        if (present(scaling)) job%scaling = scaling
        if (present(average)) job%average = average
        if (present(window_parameter)) job%window_parameter = window_parameter
        if (present(window_parameter2)) job%window_parameter2 = window_parameter2
        if (present(window_coefficients)) job%window_coefficients = window_coefficients
        if (present(window_values)) job%window_values = window_values
        if (present(return_onesided)) job%return_onesided = return_onesided
        call flatten_complex(x,u)
        allocate(v(0))
        call apply_to_real(shape(x),shape(x),u,v,result,frequencies,job,axis,st)
        if (present(status)) status = st
    end subroutine welch_into_complex

    subroutine csd_into_rr(x,y,result,frequencies,fs,window,nperseg,noverlap,nfft,detrend,scaling,average,axis,status, &
        window_parameter,window_parameter2,window_coefficients,window_values,return_onesided, detrend_callback)
        real(dp), intent(in) :: x(..) !! first input, rank 1-15
        real(dp), intent(in) :: y(..) !! second input, rank 1-15
        complex(dp), intent(inout) :: result(..) !! caller-allocated spectrum
        real(dp), allocatable, intent(out) :: frequencies(:) !! sample frequencies
        real(dp), intent(in), optional :: fs !! SciPy spectral option
        character(*), intent(in), optional :: window !! SciPy spectral option
        integer, intent(in), optional :: nperseg !! SciPy spectral option
        integer, intent(in), optional :: noverlap !! SciPy spectral option
        integer, intent(in), optional :: nfft !! SciPy spectral option
        character(*), intent(in), optional :: detrend !! SciPy spectral option
        character(*), intent(in), optional :: scaling !! SciPy spectral option
        character(*), intent(in), optional :: average !! SciPy spectral option
        integer, intent(in), optional :: axis !! one-based or negative; default last
        integer, intent(out), optional :: status !! zero on success
        real(dp), intent(in), optional :: window_parameter !! named window shape parameter
        real(dp), intent(in), optional :: window_parameter2 !! second named window parameter
        real(dp), intent(in), optional :: window_coefficients(:) !! general cosine coefficients
        real(dp), intent(in), optional :: window_values(:) !! explicit real window samples
        logical, intent(in), optional :: return_onesided !! false retains FFT negative frequencies
        procedure(spectral_real_detrend_callback), optional :: detrend_callback !! segment callback; excludes detrend=
        type(spectral_job) :: job
        complex(dp), allocatable :: u(:), v(:)
        integer :: st
        job%method = 3
        job%first_real = .true.
        job%second_real = .true.
        if (present(fs)) job%fs = fs
        if (present(window)) job%window = window
        if (present(nperseg)) job%nperseg = nperseg
        if (present(noverlap)) job%noverlap = noverlap
        if (present(nfft)) job%nfft = nfft
        if (present(detrend)) job%detrend = detrend
        if (present(detrend_callback)) job%real_callback => detrend_callback
        if (present(scaling)) job%scaling = scaling
        if (present(average)) job%average = average
        if (present(window_parameter)) job%window_parameter = window_parameter
        if (present(window_parameter2)) job%window_parameter2 = window_parameter2
        if (present(window_coefficients)) job%window_coefficients = window_coefficients
        if (present(window_values)) job%window_values = window_values
        if (present(return_onesided)) job%return_onesided = return_onesided
        call flatten_real(x,u)
        call flatten_real(y,v)
        call apply_to_complex(shape(x),shape(y),u,v,result,frequencies,job,axis,st)
        if (present(status)) status = st
    end subroutine csd_into_rr

    subroutine csd_into_cc(x,y,result,frequencies,fs,window,nperseg,noverlap,nfft,detrend,scaling,average,axis,status, &
        window_parameter,window_parameter2,window_coefficients,window_values,return_onesided, detrend_callback)
        complex(dp), intent(in) :: x(..) !! first input, rank 1-15
        complex(dp), intent(in) :: y(..) !! second input, rank 1-15
        complex(dp), intent(inout) :: result(..) !! caller-allocated spectrum
        real(dp), allocatable, intent(out) :: frequencies(:) !! sample frequencies
        real(dp), intent(in), optional :: fs !! SciPy spectral option
        character(*), intent(in), optional :: window !! SciPy spectral option
        integer, intent(in), optional :: nperseg !! SciPy spectral option
        integer, intent(in), optional :: noverlap !! SciPy spectral option
        integer, intent(in), optional :: nfft !! SciPy spectral option
        character(*), intent(in), optional :: detrend !! SciPy spectral option
        character(*), intent(in), optional :: scaling !! SciPy spectral option
        character(*), intent(in), optional :: average !! SciPy spectral option
        integer, intent(in), optional :: axis !! one-based or negative; default last
        integer, intent(out), optional :: status !! zero on success
        real(dp), intent(in), optional :: window_parameter !! named window shape parameter
        real(dp), intent(in), optional :: window_parameter2 !! second named window parameter
        real(dp), intent(in), optional :: window_coefficients(:) !! general cosine coefficients
        real(dp), intent(in), optional :: window_values(:) !! explicit real window samples
        logical, intent(in), optional :: return_onesided !! false retains FFT negative frequencies
        procedure(spectral_complex_detrend_callback), optional :: detrend_callback !! segment callback; excludes detrend=
        type(spectral_job) :: job
        complex(dp), allocatable :: u(:), v(:)
        integer :: st
        job%method = 3
        job%first_real = .false.
        job%second_real = .false.
        if (present(fs)) job%fs = fs
        if (present(window)) job%window = window
        if (present(nperseg)) job%nperseg = nperseg
        if (present(noverlap)) job%noverlap = noverlap
        if (present(nfft)) job%nfft = nfft
        if (present(detrend)) job%detrend = detrend
        if (present(detrend_callback)) job%complex_callback => detrend_callback
        if (present(scaling)) job%scaling = scaling
        if (present(average)) job%average = average
        if (present(window_parameter)) job%window_parameter = window_parameter
        if (present(window_parameter2)) job%window_parameter2 = window_parameter2
        if (present(window_coefficients)) job%window_coefficients = window_coefficients
        if (present(window_values)) job%window_values = window_values
        if (present(return_onesided)) job%return_onesided = return_onesided
        call flatten_complex(x,u)
        call flatten_complex(y,v)
        call apply_to_complex(shape(x),shape(y),u,v,result,frequencies,job,axis,st)
        if (present(status)) status = st
    end subroutine csd_into_cc

    subroutine csd_into_rc(x,y,result,frequencies,fs,window,nperseg,noverlap,nfft,detrend,scaling,average,axis,status, &
        window_parameter,window_parameter2,window_coefficients,window_values,return_onesided, detrend_callback)
        real(dp), intent(in) :: x(..) !! first input, rank 1-15
        complex(dp), intent(in) :: y(..) !! second input, rank 1-15
        complex(dp), intent(inout) :: result(..) !! caller-allocated spectrum
        real(dp), allocatable, intent(out) :: frequencies(:) !! sample frequencies
        real(dp), intent(in), optional :: fs !! SciPy spectral option
        character(*), intent(in), optional :: window !! SciPy spectral option
        integer, intent(in), optional :: nperseg !! SciPy spectral option
        integer, intent(in), optional :: noverlap !! SciPy spectral option
        integer, intent(in), optional :: nfft !! SciPy spectral option
        character(*), intent(in), optional :: detrend !! SciPy spectral option
        character(*), intent(in), optional :: scaling !! SciPy spectral option
        character(*), intent(in), optional :: average !! SciPy spectral option
        integer, intent(in), optional :: axis !! one-based or negative; default last
        integer, intent(out), optional :: status !! zero on success
        real(dp), intent(in), optional :: window_parameter !! named window shape parameter
        real(dp), intent(in), optional :: window_parameter2 !! second named window parameter
        real(dp), intent(in), optional :: window_coefficients(:) !! general cosine coefficients
        real(dp), intent(in), optional :: window_values(:) !! explicit real window samples
        logical, intent(in), optional :: return_onesided !! false retains FFT negative frequencies
        procedure(spectral_complex_detrend_callback), optional :: detrend_callback !! segment callback; excludes detrend=
        type(spectral_job) :: job
        complex(dp), allocatable :: u(:), v(:)
        integer :: st
        job%method = 3
        job%first_real = .true.
        job%second_real = .false.
        if (present(fs)) job%fs = fs
        if (present(window)) job%window = window
        if (present(nperseg)) job%nperseg = nperseg
        if (present(noverlap)) job%noverlap = noverlap
        if (present(nfft)) job%nfft = nfft
        if (present(detrend)) job%detrend = detrend
        if (present(detrend_callback)) job%complex_callback => detrend_callback
        if (present(scaling)) job%scaling = scaling
        if (present(average)) job%average = average
        if (present(window_parameter)) job%window_parameter = window_parameter
        if (present(window_parameter2)) job%window_parameter2 = window_parameter2
        if (present(window_coefficients)) job%window_coefficients = window_coefficients
        if (present(window_values)) job%window_values = window_values
        if (present(return_onesided)) job%return_onesided = return_onesided
        call flatten_real(x,u)
        call flatten_complex(y,v)
        call apply_to_complex(shape(x),shape(y),u,v,result,frequencies,job,axis,st)
        if (present(status)) status = st
    end subroutine csd_into_rc

    subroutine csd_into_cr(x,y,result,frequencies,fs,window,nperseg,noverlap,nfft,detrend,scaling,average,axis,status, &
        window_parameter,window_parameter2,window_coefficients,window_values,return_onesided, detrend_callback)
        complex(dp), intent(in) :: x(..) !! first input, rank 1-15
        real(dp), intent(in) :: y(..) !! second input, rank 1-15
        complex(dp), intent(inout) :: result(..) !! caller-allocated spectrum
        real(dp), allocatable, intent(out) :: frequencies(:) !! sample frequencies
        real(dp), intent(in), optional :: fs !! SciPy spectral option
        character(*), intent(in), optional :: window !! SciPy spectral option
        integer, intent(in), optional :: nperseg !! SciPy spectral option
        integer, intent(in), optional :: noverlap !! SciPy spectral option
        integer, intent(in), optional :: nfft !! SciPy spectral option
        character(*), intent(in), optional :: detrend !! SciPy spectral option
        character(*), intent(in), optional :: scaling !! SciPy spectral option
        character(*), intent(in), optional :: average !! SciPy spectral option
        integer, intent(in), optional :: axis !! one-based or negative; default last
        integer, intent(out), optional :: status !! zero on success
        real(dp), intent(in), optional :: window_parameter !! named window shape parameter
        real(dp), intent(in), optional :: window_parameter2 !! second named window parameter
        real(dp), intent(in), optional :: window_coefficients(:) !! general cosine coefficients
        real(dp), intent(in), optional :: window_values(:) !! explicit real window samples
        logical, intent(in), optional :: return_onesided !! false retains FFT negative frequencies
        procedure(spectral_complex_detrend_callback), optional :: detrend_callback !! segment callback; excludes detrend=
        type(spectral_job) :: job
        complex(dp), allocatable :: u(:), v(:)
        integer :: st
        job%method = 3
        job%first_real = .false.
        job%second_real = .true.
        if (present(fs)) job%fs = fs
        if (present(window)) job%window = window
        if (present(nperseg)) job%nperseg = nperseg
        if (present(noverlap)) job%noverlap = noverlap
        if (present(nfft)) job%nfft = nfft
        if (present(detrend)) job%detrend = detrend
        if (present(detrend_callback)) job%complex_callback => detrend_callback
        if (present(scaling)) job%scaling = scaling
        if (present(average)) job%average = average
        if (present(window_parameter)) job%window_parameter = window_parameter
        if (present(window_parameter2)) job%window_parameter2 = window_parameter2
        if (present(window_coefficients)) job%window_coefficients = window_coefficients
        if (present(window_values)) job%window_values = window_values
        if (present(return_onesided)) job%return_onesided = return_onesided
        call flatten_complex(x,u)
        call flatten_real(y,v)
        call apply_to_complex(shape(x),shape(y),u,v,result,frequencies,job,axis,st)
        if (present(status)) status = st
    end subroutine csd_into_cr

    subroutine coherence_into_rr(x,y,result,frequencies,fs,window,nperseg,noverlap,nfft,detrend,axis,status, &
        window_parameter,window_parameter2,window_coefficients,window_values, detrend_callback)
        real(dp), intent(in) :: x(..) !! first input, rank 1-15
        real(dp), intent(in) :: y(..) !! second input, rank 1-15
        real(dp), intent(inout) :: result(..) !! caller-allocated spectrum
        real(dp), allocatable, intent(out) :: frequencies(:) !! sample frequencies
        real(dp), intent(in), optional :: fs !! SciPy spectral option
        character(*), intent(in), optional :: window !! SciPy spectral option
        integer, intent(in), optional :: nperseg !! SciPy spectral option
        integer, intent(in), optional :: noverlap !! SciPy spectral option
        integer, intent(in), optional :: nfft !! SciPy spectral option
        character(*), intent(in), optional :: detrend !! SciPy spectral option
        integer, intent(in), optional :: axis !! one-based or negative; default last
        integer, intent(out), optional :: status !! zero on success
        real(dp), intent(in), optional :: window_parameter !! named window shape parameter
        real(dp), intent(in), optional :: window_parameter2 !! second named window parameter
        real(dp), intent(in), optional :: window_coefficients(:) !! general cosine coefficients
        real(dp), intent(in), optional :: window_values(:) !! explicit real window samples
        procedure(spectral_real_detrend_callback), optional :: detrend_callback !! segment callback; excludes detrend=
        type(spectral_job) :: job
        complex(dp), allocatable :: u(:), v(:)
        integer :: st
        job%method = 4
        job%first_real = .true.
        job%second_real = .true.
        if (present(fs)) job%fs = fs
        if (present(window)) job%window = window
        if (present(nperseg)) job%nperseg = nperseg
        if (present(noverlap)) job%noverlap = noverlap
        if (present(nfft)) job%nfft = nfft
        if (present(detrend)) job%detrend = detrend
        if (present(detrend_callback)) job%real_callback => detrend_callback
        if (present(window_parameter)) job%window_parameter = window_parameter
        if (present(window_parameter2)) job%window_parameter2 = window_parameter2
        if (present(window_coefficients)) job%window_coefficients = window_coefficients
        if (present(window_values)) job%window_values = window_values
        call flatten_real(x,u)
        call flatten_real(y,v)
        call apply_to_real(shape(x),shape(y),u,v,result,frequencies,job,axis,st)
        if (present(status)) status = st
    end subroutine coherence_into_rr

    subroutine coherence_into_cc(x,y,result,frequencies,fs,window,nperseg,noverlap,nfft,detrend,axis,status, &
        window_parameter,window_parameter2,window_coefficients,window_values, detrend_callback)
        complex(dp), intent(in) :: x(..) !! first input, rank 1-15
        complex(dp), intent(in) :: y(..) !! second input, rank 1-15
        real(dp), intent(inout) :: result(..) !! caller-allocated spectrum
        real(dp), allocatable, intent(out) :: frequencies(:) !! sample frequencies
        real(dp), intent(in), optional :: fs !! SciPy spectral option
        character(*), intent(in), optional :: window !! SciPy spectral option
        integer, intent(in), optional :: nperseg !! SciPy spectral option
        integer, intent(in), optional :: noverlap !! SciPy spectral option
        integer, intent(in), optional :: nfft !! SciPy spectral option
        character(*), intent(in), optional :: detrend !! SciPy spectral option
        integer, intent(in), optional :: axis !! one-based or negative; default last
        integer, intent(out), optional :: status !! zero on success
        real(dp), intent(in), optional :: window_parameter !! named window shape parameter
        real(dp), intent(in), optional :: window_parameter2 !! second named window parameter
        real(dp), intent(in), optional :: window_coefficients(:) !! general cosine coefficients
        real(dp), intent(in), optional :: window_values(:) !! explicit real window samples
        procedure(spectral_complex_detrend_callback), optional :: detrend_callback !! segment callback; excludes detrend=
        type(spectral_job) :: job
        complex(dp), allocatable :: u(:), v(:)
        integer :: st
        job%method = 4
        job%first_real = .false.
        job%second_real = .false.
        if (present(fs)) job%fs = fs
        if (present(window)) job%window = window
        if (present(nperseg)) job%nperseg = nperseg
        if (present(noverlap)) job%noverlap = noverlap
        if (present(nfft)) job%nfft = nfft
        if (present(detrend)) job%detrend = detrend
        if (present(detrend_callback)) job%complex_callback => detrend_callback
        if (present(window_parameter)) job%window_parameter = window_parameter
        if (present(window_parameter2)) job%window_parameter2 = window_parameter2
        if (present(window_coefficients)) job%window_coefficients = window_coefficients
        if (present(window_values)) job%window_values = window_values
        call flatten_complex(x,u)
        call flatten_complex(y,v)
        call apply_to_real(shape(x),shape(y),u,v,result,frequencies,job,axis,st)
        if (present(status)) status = st
    end subroutine coherence_into_cc

    subroutine coherence_into_rc(x,y,result,frequencies,fs,window,nperseg,noverlap,nfft,detrend,axis,status, &
        window_parameter,window_parameter2,window_coefficients,window_values, detrend_callback)
        real(dp), intent(in) :: x(..) !! first input, rank 1-15
        complex(dp), intent(in) :: y(..) !! second input, rank 1-15
        real(dp), intent(inout) :: result(..) !! caller-allocated spectrum
        real(dp), allocatable, intent(out) :: frequencies(:) !! sample frequencies
        real(dp), intent(in), optional :: fs !! SciPy spectral option
        character(*), intent(in), optional :: window !! SciPy spectral option
        integer, intent(in), optional :: nperseg !! SciPy spectral option
        integer, intent(in), optional :: noverlap !! SciPy spectral option
        integer, intent(in), optional :: nfft !! SciPy spectral option
        character(*), intent(in), optional :: detrend !! SciPy spectral option
        integer, intent(in), optional :: axis !! one-based or negative; default last
        integer, intent(out), optional :: status !! zero on success
        real(dp), intent(in), optional :: window_parameter !! named window shape parameter
        real(dp), intent(in), optional :: window_parameter2 !! second named window parameter
        real(dp), intent(in), optional :: window_coefficients(:) !! general cosine coefficients
        real(dp), intent(in), optional :: window_values(:) !! explicit real window samples
        procedure(spectral_complex_detrend_callback), optional :: detrend_callback !! segment callback; excludes detrend=
        type(spectral_job) :: job
        complex(dp), allocatable :: u(:), v(:)
        integer :: st
        job%method = 4
        job%first_real = .true.
        job%second_real = .false.
        if (present(fs)) job%fs = fs
        if (present(window)) job%window = window
        if (present(nperseg)) job%nperseg = nperseg
        if (present(noverlap)) job%noverlap = noverlap
        if (present(nfft)) job%nfft = nfft
        if (present(detrend)) job%detrend = detrend
        if (present(detrend_callback)) job%complex_callback => detrend_callback
        if (present(window_parameter)) job%window_parameter = window_parameter
        if (present(window_parameter2)) job%window_parameter2 = window_parameter2
        if (present(window_coefficients)) job%window_coefficients = window_coefficients
        if (present(window_values)) job%window_values = window_values
        call flatten_real(x,u)
        call flatten_complex(y,v)
        call apply_to_real(shape(x),shape(y),u,v,result,frequencies,job,axis,st)
        if (present(status)) status = st
    end subroutine coherence_into_rc

    subroutine coherence_into_cr(x,y,result,frequencies,fs,window,nperseg,noverlap,nfft,detrend,axis,status, &
        window_parameter,window_parameter2,window_coefficients,window_values, detrend_callback)
        complex(dp), intent(in) :: x(..) !! first input, rank 1-15
        real(dp), intent(in) :: y(..) !! second input, rank 1-15
        real(dp), intent(inout) :: result(..) !! caller-allocated spectrum
        real(dp), allocatable, intent(out) :: frequencies(:) !! sample frequencies
        real(dp), intent(in), optional :: fs !! SciPy spectral option
        character(*), intent(in), optional :: window !! SciPy spectral option
        integer, intent(in), optional :: nperseg !! SciPy spectral option
        integer, intent(in), optional :: noverlap !! SciPy spectral option
        integer, intent(in), optional :: nfft !! SciPy spectral option
        character(*), intent(in), optional :: detrend !! SciPy spectral option
        integer, intent(in), optional :: axis !! one-based or negative; default last
        integer, intent(out), optional :: status !! zero on success
        real(dp), intent(in), optional :: window_parameter !! named window shape parameter
        real(dp), intent(in), optional :: window_parameter2 !! second named window parameter
        real(dp), intent(in), optional :: window_coefficients(:) !! general cosine coefficients
        real(dp), intent(in), optional :: window_values(:) !! explicit real window samples
        procedure(spectral_complex_detrend_callback), optional :: detrend_callback !! segment callback; excludes detrend=
        type(spectral_job) :: job
        complex(dp), allocatable :: u(:), v(:)
        integer :: st
        job%method = 4
        job%first_real = .false.
        job%second_real = .true.
        if (present(fs)) job%fs = fs
        if (present(window)) job%window = window
        if (present(nperseg)) job%nperseg = nperseg
        if (present(noverlap)) job%noverlap = noverlap
        if (present(nfft)) job%nfft = nfft
        if (present(detrend)) job%detrend = detrend
        if (present(detrend_callback)) job%complex_callback => detrend_callback
        if (present(window_parameter)) job%window_parameter = window_parameter
        if (present(window_parameter2)) job%window_parameter2 = window_parameter2
        if (present(window_coefficients)) job%window_coefficients = window_coefficients
        if (present(window_values)) job%window_values = window_values
        call flatten_complex(x,u)
        call flatten_real(y,v)
        call apply_to_real(shape(x),shape(y),u,v,result,frequencies,job,axis,st)
        if (present(status)) status = st
    end subroutine coherence_into_cr

    ! END GENERATED PUBLIC WRAPPERS

end module scifort_signal_spectral_nd
