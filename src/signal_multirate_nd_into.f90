! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Rank-independent, caller-allocated multirate transforms.  A single lane
! traversal applies the established one-dimensional kernels along any axis.
! Entire output is committed only after all lanes succeed.
module scifort_signal_multirate_nd_into
    use scifort_kinds, only : dp
    use scifort_signal_resample, only : upfirdn, resample_poly, decimate, &
        resample_status_success, resample_status_invalid_input
    use scifort_signal_fourier_resample, only : resample, &
        fourier_resample_status_success, fourier_resample_status_invalid_input, &
        real_frequency_window_callback, complex_frequency_window_callback, &
        evaluate_real_frequency_window, evaluate_complex_frequency_window
    implicit none
    private
    public :: resample_into, resample_poly_into, upfirdn_into, decimate_into

    ! Immutable job arguments replace nested callbacks. This permits native
    ! module procedures without compiler-generated executable-stack trampolines.
    type :: multirate_job
        integer :: method = 0
        integer :: num = 0
        integer :: up = 1
        integer :: down = 1
        integer, allocatable :: n
        real(dp), allocatable :: beta, window_parameter, window_parameter2
        real(dp), allocatable :: window_coefficients(:)
        complex(dp) :: cval = cmplx(0.0_dp,0.0_dp,dp)
        logical :: has_cval = .false.
        real(dp), allocatable :: taps_real(:)
        real(dp), allocatable :: sos_real(:,:)
        complex(dp), allocatable :: sos_complex(:,:)
        complex(dp), allocatable :: taps_complex(:)
        complex(dp), allocatable :: numerator(:), denominator(:)
        real(dp), allocatable :: frequency_window(:)
        complex(dp), allocatable :: complex_window(:)
        logical :: real_input = .false.
        real(dp), allocatable :: t(:)
        character(:), allocatable :: window, padtype, mode, ftype, domain
        logical, allocatable :: zero_phase
    end type multirate_job

    interface resample_into
        module procedure resample_into_real, resample_into_complex
    end interface
    interface resample_poly_into
        module procedure poly_into_rr, poly_into_cc, poly_into_rc, poly_into_cr
    end interface
    interface upfirdn_into
        module procedure upfirdn_into_rr, upfirdn_into_rc, upfirdn_into_cr, upfirdn_into_cc
    end interface
    interface decimate_into
        module procedure decimate_into_real, decimate_into_complex
        module procedure decimate_into_real_complex_tf
        module procedure decimate_into_complex_sos, decimate_into_real_complex_sos
    end interface

contains

    subroutine check_shapes(x, y, axis, default_axis, a, status)
        integer, intent(in) :: default_axis !! 1 for resample, -1 for filter/decimate
        ! Validate extent vectors before any lane is processed.
        integer, intent(in) :: x(:) !! input dimensions
        integer, intent(in) :: y(:) !! output dimensions
        integer, intent(in), optional :: axis !! one-based or negative axis
        integer, intent(out) :: a !! normalized processing axis
        integer, intent(out) :: status !! zero if shapes are compatible
        integer :: j, ndim
        ndim = size(x)
        a = default_axis
        status = resample_status_invalid_input
        if (ndim < 1 .or. ndim > 15 .or. size(y) /= ndim) return
        if (any(x <= 0) .or. any(y <= 0)) return
        if (present(axis)) a = axis
        if (a < 0) a = a + ndim + 1
        if (a < 1 .or. a > ndim) return
        do j = 1, ndim
            if (j /= a .and. x(j) /= y(j)) return
        end do
        status = resample_status_success
    end subroutine check_shapes

    function run_lane(v, job, status) result(w)
        complex(dp), intent(in) :: v(:) !! input complex signal lane
        type(multirate_job), intent(in) :: job !! fixed processing parameters
        integer, intent(out) :: status !! zero on success
        complex(dp), allocatable :: w(:)
        select case (job%method)
        case (1) ! Periodic Fourier resampling.
            w = resample(v,job%num,job%window,job%t,domain=job%domain,status=status, &
                frequency_window=job%frequency_window,complex_window=job%complex_window, &
                real_input=job%real_input, &
                window_parameter=job%window_parameter, window_parameter2=job%window_parameter2, &
                window_coefficients=job%window_coefficients)
        case (2) ! Polyphase resampling.
            if (job%has_cval) then
                if (allocated(job%taps_complex)) then
                    w = resample_poly(v,job%up,job%down,job%taps_complex, &
                        job%window,job%beta,job%padtype,job%cval,status=status, &
                        window_parameter=job%window_parameter, &
                        window_parameter2=job%window_parameter2, &
                        window_coefficients=job%window_coefficients)
                else
                    w = resample_poly(v,job%up,job%down,job%taps_real, &
                        job%window,job%beta,job%padtype,job%cval,status=status, &
                        window_parameter=job%window_parameter, &
                        window_parameter2=job%window_parameter2, &
                        window_coefficients=job%window_coefficients)
                end if
            else
                if (allocated(job%taps_complex)) then
                    w = resample_poly(v,job%up,job%down,job%taps_complex, &
                        job%window,job%beta,job%padtype,status=status, &
                        window_parameter=job%window_parameter, &
                        window_parameter2=job%window_parameter2, &
                        window_coefficients=job%window_coefficients)
                else
                    w = resample_poly(v,job%up,job%down,job%taps_real, &
                        job%window,job%beta,job%padtype,status=status, &
                        window_parameter=job%window_parameter, &
                        window_parameter2=job%window_parameter2, &
                        window_coefficients=job%window_coefficients)
                end if
            end if
        case (3) ! FIR interpolation and decimation.
            if (allocated(job%taps_complex)) then
                w = upfirdn(job%taps_complex,v,job%up,job%down,job%mode,job%cval,status=status)
            else
                w = upfirdn(job%taps_real,v,job%up,job%down,job%mode,job%cval,status=status)
            end if
        case (4) ! FIR or IIR decimation.
            if (allocated(job%sos_complex)) then
                w = decimate(v,job%down,sos=job%sos_complex, &
                    zero_phase=job%zero_phase,status=status)
            else
                w = decimate(v,job%down,job%n,job%ftype,job%zero_phase,status=status, &
                    taps=job%taps_real,sos=job%sos_real, &
                    b=job%numerator,a=job%denominator)
            end if
        case default
            status = resample_status_invalid_input
            allocate(w(0))
        end select
    end function run_lane

    subroutine map_lanes(src, dst, dims, a, job, status)
        complex(dp), intent(in) :: src(:) !! flat Fortran-order input
        complex(dp), intent(out) :: dst(:) !! flat output staging buffer
        integer, intent(in) :: dims(:) !! original input dimensions
        integer, intent(in) :: a !! positive processing axis
        type(multirate_job), intent(in) :: job !! transform arguments
        integer, intent(out) :: status !! zero only if every lane succeeds
        integer :: inner, outer, nx, ny, i, j, k, st
        integer :: il, io
        complex(dp), allocatable :: v(:), w(:)
        status = resample_status_invalid_input
        inner = product(dims(:a-1))
        outer = product(dims(a+1:))
        nx = dims(a)
        if (inner <= 0 .or. outer <= 0 .or. nx <= 0) return
        if (mod(size(src), inner*outer) /= 0) return
        if (mod(size(dst), inner*outer) /= 0) return
        ny = size(dst)/(inner*outer)
        if (ny <= 0) return
        allocate(v(nx))
        do k = 0, outer-1
            do i = 0, inner-1
                do j = 0, nx-1
                    il = k*inner*nx + j*inner + i + 1
                    v(j+1) = src(il)
                end do
                w = run_lane(v, job, st)
                if (st /= resample_status_success) return
                if (size(w) /= ny) return
                do j = 0, ny-1
                    io = k*inner*ny + j*inner + i + 1
                    dst(io) = w(j+1)
                end do
            end do
        end do
        status = resample_status_success
    end subroutine map_lanes

    subroutine flatten_real(x, v)
        real(dp), intent(in) :: x(..) !! input of any rank 1-15
        complex(dp), allocatable, intent(out) :: v(:) !! flat complex copy
        select rank(x)
        ! BEGIN GENERATED FLATTEN REAL
        rank (1)
            v = cmplx(reshape(x,[size(x)]),0.0_dp,dp)
        rank (2)
            v = cmplx(reshape(x,[size(x)]),0.0_dp,dp)
        rank (3)
            v = cmplx(reshape(x,[size(x)]),0.0_dp,dp)
        rank (4)
            v = cmplx(reshape(x,[size(x)]),0.0_dp,dp)
        rank (5)
            v = cmplx(reshape(x,[size(x)]),0.0_dp,dp)
        rank (6)
            v = cmplx(reshape(x,[size(x)]),0.0_dp,dp)
        rank (7)
            v = cmplx(reshape(x,[size(x)]),0.0_dp,dp)
        rank (8)
            v = cmplx(reshape(x,[size(x)]),0.0_dp,dp)
        rank (9)
            v = cmplx(reshape(x,[size(x)]),0.0_dp,dp)
        rank (10)
            v = cmplx(reshape(x,[size(x)]),0.0_dp,dp)
        rank (11)
            v = cmplx(reshape(x,[size(x)]),0.0_dp,dp)
        rank (12)
            v = cmplx(reshape(x,[size(x)]),0.0_dp,dp)
        rank (13)
            v = cmplx(reshape(x,[size(x)]),0.0_dp,dp)
        rank (14)
            v = cmplx(reshape(x,[size(x)]),0.0_dp,dp)
        rank (15)
            v = cmplx(reshape(x,[size(x)]),0.0_dp,dp)
                ! END GENERATED FLATTEN REAL
        rank default
            allocate(v(0))
        end select
    end subroutine flatten_real

    subroutine flatten_complex(x, v)
        complex(dp), intent(in) :: x(..) !! input of any rank 1-15
        complex(dp), allocatable, intent(out) :: v(:) !! flat complex copy
        select rank(x)
        ! BEGIN GENERATED FLATTEN COMPLEX
        rank (1)
            v = reshape(x,[size(x)])
        rank (2)
            v = reshape(x,[size(x)])
        rank (3)
            v = reshape(x,[size(x)])
        rank (4)
            v = reshape(x,[size(x)])
        rank (5)
            v = reshape(x,[size(x)])
        rank (6)
            v = reshape(x,[size(x)])
        rank (7)
            v = reshape(x,[size(x)])
        rank (8)
            v = reshape(x,[size(x)])
        rank (9)
            v = reshape(x,[size(x)])
        rank (10)
            v = reshape(x,[size(x)])
        rank (11)
            v = reshape(x,[size(x)])
        rank (12)
            v = reshape(x,[size(x)])
        rank (13)
            v = reshape(x,[size(x)])
        rank (14)
            v = reshape(x,[size(x)])
        rank (15)
            v = reshape(x,[size(x)])
                ! END GENERATED FLATTEN COMPLEX
        rank default
            allocate(v(0))
        end select
    end subroutine flatten_complex

    subroutine scatter_real(y, v)
        real(dp), intent(inout) :: y(..) !! complete caller-provided result array
        complex(dp), intent(in) :: v(:) !! staged output data
        select rank(y)
        ! BEGIN GENERATED SCATTER REAL
        rank (1)
            y = reshape(real(v,dp),shape(y))
        rank (2)
            y = reshape(real(v,dp),shape(y))
        rank (3)
            y = reshape(real(v,dp),shape(y))
        rank (4)
            y = reshape(real(v,dp),shape(y))
        rank (5)
            y = reshape(real(v,dp),shape(y))
        rank (6)
            y = reshape(real(v,dp),shape(y))
        rank (7)
            y = reshape(real(v,dp),shape(y))
        rank (8)
            y = reshape(real(v,dp),shape(y))
        rank (9)
            y = reshape(real(v,dp),shape(y))
        rank (10)
            y = reshape(real(v,dp),shape(y))
        rank (11)
            y = reshape(real(v,dp),shape(y))
        rank (12)
            y = reshape(real(v,dp),shape(y))
        rank (13)
            y = reshape(real(v,dp),shape(y))
        rank (14)
            y = reshape(real(v,dp),shape(y))
        rank (15)
            y = reshape(real(v,dp),shape(y))
                ! END GENERATED SCATTER REAL
        rank default
            return
        end select
    end subroutine scatter_real

    subroutine scatter_complex(y, v)
        complex(dp), intent(inout) :: y(..) !! complete caller-provided result array
        complex(dp), intent(in) :: v(:) !! staged output data
        select rank(y)
        ! BEGIN GENERATED SCATTER COMPLEX
        rank (1)
            y = reshape(v,shape(y))
        rank (2)
            y = reshape(v,shape(y))
        rank (3)
            y = reshape(v,shape(y))
        rank (4)
            y = reshape(v,shape(y))
        rank (5)
            y = reshape(v,shape(y))
        rank (6)
            y = reshape(v,shape(y))
        rank (7)
            y = reshape(v,shape(y))
        rank (8)
            y = reshape(v,shape(y))
        rank (9)
            y = reshape(v,shape(y))
        rank (10)
            y = reshape(v,shape(y))
        rank (11)
            y = reshape(v,shape(y))
        rank (12)
            y = reshape(v,shape(y))
        rank (13)
            y = reshape(v,shape(y))
        rank (14)
            y = reshape(v,shape(y))
        rank (15)
            y = reshape(v,shape(y))
                ! END GENERATED SCATTER COMPLEX
        rank default
            return
        end select
    end subroutine scatter_complex

    subroutine apply_rr(x, y, axis, default_axis, job, status)
        real(dp), intent(in) :: x(..) !! real input array, rank 1-15
        real(dp), intent(inout) :: y(..) !! preallocated real output array
        integer, intent(in), optional :: axis !! one-based or negative transform axis
        integer, intent(in) :: default_axis !! default processing axis
        type(multirate_job), intent(in) :: job !! transform arguments
        integer, intent(out) :: status !! zero on success
        complex(dp), allocatable :: src(:), dst(:)
        integer :: a
        call check_shapes(shape(x),shape(y),axis,default_axis,a,status)
        if (status /= resample_status_success) return
        call flatten_real(x,src)
        allocate(dst(size(y)))
        call map_lanes(src,dst,shape(x),a,job,status)
        if (status == resample_status_success) call scatter_real(y,dst)
    end subroutine apply_rr

    subroutine apply_rc(x, y, axis, default_axis, job, status)
        real(dp), intent(in) :: x(..) !! real input array, rank 1-15
        complex(dp), intent(inout) :: y(..) !! preallocated complex output
        integer, intent(in), optional :: axis !! processing axis
        integer, intent(in) :: default_axis !! default processing axis
        type(multirate_job), intent(in) :: job !! transform arguments
        integer, intent(out) :: status !! zero on success
        complex(dp), allocatable :: src(:), dst(:)
        integer :: a
        call check_shapes(shape(x),shape(y),axis,default_axis,a,status)
        if (status /= resample_status_success) return
        call flatten_real(x,src)
        allocate(dst(size(y)))
        call map_lanes(src,dst,shape(x),a,job,status)
        if (status == resample_status_success) call scatter_complex(y,dst)
    end subroutine apply_rc

    subroutine apply_cc(x, y, axis, default_axis, job, status)
        complex(dp), intent(in) :: x(..) !! complex input array, rank 1-15
        complex(dp), intent(inout) :: y(..) !! preallocated complex output
        integer, intent(in), optional :: axis !! processing axis
        integer, intent(in) :: default_axis !! default processing axis
        type(multirate_job), intent(in) :: job !! transform arguments
        integer, intent(out) :: status !! zero on success
        complex(dp), allocatable :: src(:), dst(:)
        integer :: a
        call check_shapes(shape(x),shape(y),axis,default_axis,a,status)
        if (status /= resample_status_success) return
        call flatten_complex(x,src)
        allocate(dst(size(y)))
        call map_lanes(src,dst,shape(x),a,job,status)
        if (status == resample_status_success) call scatter_complex(y,dst)
    end subroutine apply_cc

    subroutine resample_into_real(x, y, num, window, t, output_times, axis, status, frequency_window, window_function, &
            window_parameter, window_parameter2, window_coefficients)
        real(dp), intent(in) :: x(..) !! input array, rank 1-15
        real(dp), intent(inout) :: y(..) !! caller-allocated output array
        integer, intent(in) :: num !! desired output length along axis
        character(*), intent(in), optional :: window !! spectral taper
        real(dp), intent(in), optional :: t(:) !! sample-time grid
        real(dp), allocatable, intent(out), optional :: output_times(:) !! resampled times
        integer, intent(in), optional :: axis !! one-based or negative axis; default 1
        integer, intent(out), optional :: status !! zero on success
        real(dp), intent(in), optional :: frequency_window(:) !! FFT-order real bin weights
        procedure(real_frequency_window_callback), optional :: window_function
        real(dp), intent(in), optional :: window_parameter !! first named window shape parameter
        real(dp), intent(in), optional :: window_parameter2 !! second named window shape parameter
        real(dp), intent(in), optional :: window_coefficients(:) !! general-cosine series coefficients
        type(multirate_job) :: job
        integer :: st, a
        job%method = 1
        job%real_input = .true.
        job%num = num
        if (present(window)) job%window = window
        if (present(window_parameter)) job%window_parameter = window_parameter
        if (present(window_parameter2)) job%window_parameter2 = window_parameter2
        if (present(window_coefficients)) job%window_coefficients = window_coefficients
        if (present(t)) job%t = t
        if (present(frequency_window)) job%frequency_window = frequency_window
        call check_shapes(shape(x),shape(y),axis,1,a,st)
        if (st == resample_status_success) then
            if (size(y,a) /= num) st = resample_status_invalid_input
        end if
        if (st == resample_status_success .and. present(window_function)) then
            if (present(window) .or. present(frequency_window)) then
                st = resample_status_invalid_input
            else
                call evaluate_real_frequency_window(size(x,a),window_function,job%frequency_window,st)
            end if
        end if
        if (st == resample_status_success) call apply_rr(x,y,axis,1,job,st)
        if (present(output_times)) then
            allocate(output_times(0))
            if (st == fourier_resample_status_success .and. present(t)) then
                deallocate(output_times)
                allocate(output_times(num))
                block
                    integer :: j, a
                    a = 1
                    if (present(axis)) a = axis
                    if (a < 0) a = a + rank(x) + 1
                    do j = 1, num
                        output_times(j) = t(1) + (t(2)-t(1))* &
                            real(size(x,a),dp)/real(num,dp)*real(j-1,dp)
                    end do
                end block
            end if
        end if
        if (present(status)) status = st
    end subroutine resample_into_real

    subroutine resample_into_complex(x, y, num, window, t, output_times, domain, axis, status, &
            frequency_window, complex_window, window_function, complex_window_function, &
            window_parameter, window_parameter2, window_coefficients)
        complex(dp), intent(in) :: x(..) !! input array, rank 1-15
        complex(dp), intent(inout) :: y(..) !! caller-allocated output array
        integer, intent(in) :: num !! desired output length along axis
        character(*), intent(in), optional :: window !! spectral taper
        real(dp), intent(in), optional :: t(:) !! sample-time grid
        real(dp), allocatable, intent(out), optional :: output_times(:) !! resampled times
        character(*), intent(in), optional :: domain !! time or freq
        integer, intent(in), optional :: axis !! one-based or negative axis; default 1
        integer, intent(out), optional :: status !! zero on success
        real(dp), intent(in), optional :: frequency_window(:) !! FFT-order real bin weights
        complex(dp), intent(in), optional :: complex_window(:) !! FFT-order complex bin weights
        procedure(real_frequency_window_callback), optional :: window_function
        procedure(complex_frequency_window_callback), optional :: complex_window_function
        real(dp), intent(in), optional :: window_parameter !! first named window shape parameter
        real(dp), intent(in), optional :: window_parameter2 !! second named window shape parameter
        real(dp), intent(in), optional :: window_coefficients(:) !! general-cosine series coefficients
        type(multirate_job) :: job
        integer :: st, a
        job%method = 1
        job%num = num
        if (present(window)) job%window = window
        if (present(window_parameter)) job%window_parameter = window_parameter
        if (present(window_parameter2)) job%window_parameter2 = window_parameter2
        if (present(window_coefficients)) job%window_coefficients = window_coefficients
        if (present(t)) job%t = t
        if (present(frequency_window)) job%frequency_window = frequency_window
        if (present(domain)) job%domain = domain
        if (present(complex_window)) job%complex_window = complex_window
        call check_shapes(shape(x),shape(y),axis,1,a,st)
        if (st == resample_status_success) then
            if (size(y,a) /= num) st = resample_status_invalid_input
        end if
        if (st == resample_status_success .and. present(window_function)) then
            if (present(window) .or. present(frequency_window) .or. &
                present(complex_window) .or. present(complex_window_function)) then
                st = resample_status_invalid_input
            else
                call evaluate_real_frequency_window(size(x,a),window_function,job%frequency_window,st)
            end if
        else if (st == resample_status_success .and. present(complex_window_function)) then
            if (present(window) .or. present(frequency_window) .or. present(complex_window)) then
                st = resample_status_invalid_input
            else
                call evaluate_complex_frequency_window(size(x,a),complex_window_function,job%complex_window,st)
            end if
        end if
        if (st == resample_status_success) call apply_cc(x,y,axis,1,job,st)
        if (present(output_times)) then
            allocate(output_times(0))
            if (st == fourier_resample_status_success .and. present(t)) then
                deallocate(output_times)
                allocate(output_times(num))
                block
                    integer :: j, a
                    a = 1
                    if (present(axis)) a = axis
                    if (a < 0) a = a + rank(x) + 1
                    do j = 1, num
                        output_times(j) = t(1) + (t(2)-t(1))* &
                            real(size(x,a),dp)/real(num,dp)*real(j-1,dp)
                    end do
                end block
            end if
        end if
        if (present(status)) status = st
    end subroutine resample_into_complex

    subroutine poly_into_rr(x, y, up, down, taps, window, beta, padtype, cval, axis, status, &
            window_parameter, window_parameter2, window_coefficients)
        real(dp), intent(in) :: x(..) !! input array, rank 1-15
        real(dp), intent(inout) :: y(..) !! preallocated result
        integer, intent(in) :: up !! upsampling factor
        integer, intent(in) :: down !! downsampling factor
        real(dp), intent(in), optional :: taps(:) !! FIR coefficients
        character(*), intent(in), optional :: window !! FIR window
        real(dp), intent(in), optional :: beta !! Kaiser beta
        character(*), intent(in), optional :: padtype !! padding mode
        real(dp), intent(in), optional :: cval !! constant padding
        integer, intent(in), optional :: axis !! default first axis
        integer, intent(out), optional :: status !! zero on success
        real(dp), intent(in), optional :: window_parameter !! first window shape parameter
        real(dp), intent(in), optional :: window_parameter2 !! second window shape parameter
        real(dp), intent(in), optional :: window_coefficients(:) !! general cosine window weights
        type(multirate_job) :: job
        integer :: st
        job%method = 2
        job%up = up
        job%down = down
        if (present(taps)) job%taps_real = taps
        if (present(window)) job%window = window
        if (present(beta)) job%beta = beta
        if (present(window_parameter)) job%window_parameter = window_parameter
        if (present(window_parameter2)) job%window_parameter2 = window_parameter2
        if (present(window_coefficients)) job%window_coefficients = window_coefficients
        if (present(padtype)) job%padtype = padtype
        if (present(cval)) then
            job%has_cval = .true.
            job%cval = cmplx(cval,0.0_dp,dp)
        end if
        call apply_rr(x,y,axis,1,job,st)
        if (present(status)) status = st
    end subroutine poly_into_rr

    subroutine poly_into_cc(x, y, up, down, taps, window, beta, padtype, cval, axis, status, &
            window_parameter, window_parameter2, window_coefficients)
        complex(dp), intent(in) :: x(..) !! input array, rank 1-15
        complex(dp), intent(inout) :: y(..) !! preallocated result
        integer, intent(in) :: up !! upsampling factor
        integer, intent(in) :: down !! downsampling factor
        real(dp), intent(in), optional :: taps(:) !! FIR coefficients
        character(*), intent(in), optional :: window !! FIR window
        real(dp), intent(in), optional :: beta !! Kaiser beta
        character(*), intent(in), optional :: padtype !! padding mode
        complex(dp), intent(in), optional :: cval !! constant padding
        integer, intent(in), optional :: axis !! default first axis
        integer, intent(out), optional :: status !! zero on success
        real(dp), intent(in), optional :: window_parameter !! first window shape parameter
        real(dp), intent(in), optional :: window_parameter2 !! second window shape parameter
        real(dp), intent(in), optional :: window_coefficients(:) !! general cosine window weights
        type(multirate_job) :: job
        integer :: st
        job%method = 2
        job%up = up
        job%down = down
        if (present(taps)) job%taps_real = taps
        if (present(window)) job%window = window
        if (present(beta)) job%beta = beta
        if (present(window_parameter)) job%window_parameter = window_parameter
        if (present(window_parameter2)) job%window_parameter2 = window_parameter2
        if (present(window_coefficients)) job%window_coefficients = window_coefficients
        if (present(padtype)) job%padtype = padtype
        if (present(cval)) then
            job%has_cval = .true.
            job%cval = cval
        end if
        call apply_cc(x,y,axis,1,job,st)
        if (present(status)) status = st
    end subroutine poly_into_cc

    subroutine poly_into_rc(x, y, up, down, taps, window, beta, padtype, cval, axis, status, &
            window_parameter, window_parameter2, window_coefficients)
        real(dp), intent(in) :: x(..) !! input array, rank 1-15
        complex(dp), intent(inout) :: y(..) !! preallocated result
        integer, intent(in) :: up !! upsampling factor
        integer, intent(in) :: down !! downsampling factor
        complex(dp), intent(in) :: taps(:) !! FIR coefficients
        character(*), intent(in), optional :: window !! FIR window
        real(dp), intent(in), optional :: beta !! Kaiser beta
        character(*), intent(in), optional :: padtype !! padding mode
        complex(dp), intent(in), optional :: cval !! constant padding
        integer, intent(in), optional :: axis !! default first axis
        integer, intent(out), optional :: status !! zero on success
        real(dp), intent(in), optional :: window_parameter !! first window shape parameter
        real(dp), intent(in), optional :: window_parameter2 !! second window shape parameter
        real(dp), intent(in), optional :: window_coefficients(:) !! general cosine window weights
        type(multirate_job) :: job
        integer :: st
        job%method = 2
        job%up = up
        job%down = down
        job%taps_complex = taps
        if (present(window)) job%window = window
        if (present(beta)) job%beta = beta
        if (present(window_parameter)) job%window_parameter = window_parameter
        if (present(window_parameter2)) job%window_parameter2 = window_parameter2
        if (present(window_coefficients)) job%window_coefficients = window_coefficients
        if (present(padtype)) job%padtype = padtype
        if (present(cval)) then
            job%has_cval = .true.
            job%cval = cval
        end if
        call apply_rc(x,y,axis,1,job,st)
        if (present(status)) status = st
    end subroutine poly_into_rc

    subroutine poly_into_cr(x, y, up, down, taps, window, beta, padtype, cval, axis, status, &
            window_parameter, window_parameter2, window_coefficients)
        complex(dp), intent(in) :: x(..) !! input array, rank 1-15
        complex(dp), intent(inout) :: y(..) !! preallocated result
        integer, intent(in) :: up !! upsampling factor
        integer, intent(in) :: down !! downsampling factor
        complex(dp), intent(in) :: taps(:) !! FIR coefficients
        character(*), intent(in), optional :: window !! FIR window
        real(dp), intent(in), optional :: beta !! Kaiser beta
        character(*), intent(in), optional :: padtype !! padding mode
        complex(dp), intent(in), optional :: cval !! constant padding
        integer, intent(in), optional :: axis !! default first axis
        integer, intent(out), optional :: status !! zero on success
        real(dp), intent(in), optional :: window_parameter !! first window shape parameter
        real(dp), intent(in), optional :: window_parameter2 !! second window shape parameter
        real(dp), intent(in), optional :: window_coefficients(:) !! general cosine window weights
        type(multirate_job) :: job
        integer :: st
        job%method = 2
        job%up = up
        job%down = down
        job%taps_complex = taps
        if (present(window)) job%window = window
        if (present(beta)) job%beta = beta
        if (present(window_parameter)) job%window_parameter = window_parameter
        if (present(window_parameter2)) job%window_parameter2 = window_parameter2
        if (present(window_coefficients)) job%window_coefficients = window_coefficients
        if (present(padtype)) job%padtype = padtype
        if (present(cval)) then
            job%has_cval = .true.
            job%cval = cval
        end if
        call apply_cc(x,y,axis,1,job,st)
        if (present(status)) status = st
    end subroutine poly_into_cr

    subroutine upfirdn_into_rr(h, x, y, up, down, mode, cval, axis, status)
        real(dp), intent(in) :: h(:) !! FIR taps
        real(dp), intent(in) :: x(..) !! input array, rank 1-15
        real(dp), intent(inout) :: y(..) !! preallocated result
        integer, intent(in), optional :: up !! upsampling factor
        integer, intent(in), optional :: down !! downsampling factor
        character(*), intent(in), optional :: mode !! padding mode
        real(dp), intent(in), optional :: cval !! constant padding
        integer, intent(in), optional :: axis !! default last axis
        integer, intent(out), optional :: status !! zero on success
        type(multirate_job) :: job
        integer :: st
        job%method = 3
        job%taps_real = h
        if (present(up)) job%up = up
        if (present(down)) job%down = down
        if (present(mode)) job%mode = mode
        if (present(cval)) job%cval = cmplx(cval,0.0_dp,dp)
        call apply_rr(x,y,axis,-1,job,st)
        if (present(status)) status = st
    end subroutine upfirdn_into_rr

    subroutine upfirdn_into_rc(h, x, y, up, down, mode, cval, axis, status)
        real(dp), intent(in) :: h(:) !! FIR taps
        complex(dp), intent(in) :: x(..) !! input array, rank 1-15
        complex(dp), intent(inout) :: y(..) !! preallocated result
        integer, intent(in), optional :: up !! upsampling factor
        integer, intent(in), optional :: down !! downsampling factor
        character(*), intent(in), optional :: mode !! padding mode
        complex(dp), intent(in), optional :: cval !! constant padding
        integer, intent(in), optional :: axis !! default last axis
        integer, intent(out), optional :: status !! zero on success
        type(multirate_job) :: job
        integer :: st
        job%method = 3
        job%taps_real = h
        if (present(up)) job%up = up
        if (present(down)) job%down = down
        if (present(mode)) job%mode = mode
        if (present(cval)) job%cval = cval
        call apply_cc(x,y,axis,-1,job,st)
        if (present(status)) status = st
    end subroutine upfirdn_into_rc

    subroutine upfirdn_into_cr(h, x, y, up, down, mode, cval, axis, status)
        complex(dp), intent(in) :: h(:) !! FIR taps
        real(dp), intent(in) :: x(..) !! input array, rank 1-15
        complex(dp), intent(inout) :: y(..) !! preallocated result
        integer, intent(in), optional :: up !! upsampling factor
        integer, intent(in), optional :: down !! downsampling factor
        character(*), intent(in), optional :: mode !! padding mode
        complex(dp), intent(in), optional :: cval !! constant padding
        integer, intent(in), optional :: axis !! default last axis
        integer, intent(out), optional :: status !! zero on success
        type(multirate_job) :: job
        integer :: st
        job%method = 3
        job%taps_complex = h
        if (present(up)) job%up = up
        if (present(down)) job%down = down
        if (present(mode)) job%mode = mode
        if (present(cval)) job%cval = cval
        call apply_rc(x,y,axis,-1,job,st)
        if (present(status)) status = st
    end subroutine upfirdn_into_cr

    subroutine upfirdn_into_cc(h, x, y, up, down, mode, cval, axis, status)
        complex(dp), intent(in) :: h(:) !! FIR taps
        complex(dp), intent(in) :: x(..) !! input array, rank 1-15
        complex(dp), intent(inout) :: y(..) !! preallocated result
        integer, intent(in), optional :: up !! upsampling factor
        integer, intent(in), optional :: down !! downsampling factor
        character(*), intent(in), optional :: mode !! padding mode
        complex(dp), intent(in), optional :: cval !! constant padding
        integer, intent(in), optional :: axis !! default last axis
        integer, intent(out), optional :: status !! zero on success
        type(multirate_job) :: job
        integer :: st
        job%method = 3
        job%taps_complex = h
        if (present(up)) job%up = up
        if (present(down)) job%down = down
        if (present(mode)) job%mode = mode
        if (present(cval)) job%cval = cval
        call apply_cc(x,y,axis,-1,job,st)
        if (present(status)) status = st
    end subroutine upfirdn_into_cc

    subroutine decimate_into_real(x, y, q, n, ftype, zero_phase, axis, status, taps, sos, b, a)
        real(dp), intent(in) :: x(..) !! input array, rank 1-15
        real(dp), intent(inout) :: y(..) !! preallocated result
        integer, intent(in) :: q !! downsampling factor
        integer, intent(in), optional :: n !! filter order
        character(*), intent(in), optional :: ftype !! fir or iir
        logical, intent(in), optional :: zero_phase !! filtering phase
        integer, intent(in), optional :: axis !! default last axis
        integer, intent(out), optional :: status !! zero on success
        real(dp), intent(in), optional :: taps(:) !! optional real FIR taps
        real(dp), intent(in), optional :: sos(:,:) !! optional real IIR SOS array
        real(dp), intent(in), optional :: b(:) !! transfer-function numerator coefficients
        real(dp), intent(in), optional :: a(:) !! transfer-function denominator coefficients
        type(multirate_job) :: job
        integer :: st
        job%method = 4
        job%down = q
        if (present(n)) job%n = n
        if (present(ftype)) job%ftype = ftype
        if (present(zero_phase)) job%zero_phase = zero_phase
        if (present(taps)) job%taps_real = taps
        if (present(sos)) job%sos_real = sos
        if (present(b)) job%numerator = cmplx(b,0.0_dp,dp)
        if (present(a)) job%denominator = cmplx(a,0.0_dp,dp)
        call apply_rr(x,y,axis,-1,job,st)
        if (present(status)) status = st
    end subroutine decimate_into_real

    subroutine decimate_into_complex(x, y, q, n, ftype, zero_phase, axis, status, taps, sos, b, a)
        complex(dp), intent(in) :: x(..) !! input array, rank 1-15
        complex(dp), intent(inout) :: y(..) !! preallocated result
        integer, intent(in) :: q !! downsampling factor
        integer, intent(in), optional :: n !! filter order
        character(*), intent(in), optional :: ftype !! fir or iir
        logical, intent(in), optional :: zero_phase !! filtering phase
        integer, intent(in), optional :: axis !! default last axis
        integer, intent(out), optional :: status !! zero on success
        real(dp), intent(in), optional :: taps(:) !! optional real FIR taps
        real(dp), intent(in), optional :: sos(:,:) !! optional real IIR SOS array
        complex(dp), intent(in), optional :: b(:) !! transfer-function numerator coefficients
        complex(dp), intent(in), optional :: a(:) !! transfer-function denominator coefficients
        type(multirate_job) :: job
        integer :: st
        job%method = 4
        job%down = q
        if (present(n)) job%n = n
        if (present(ftype)) job%ftype = ftype
        if (present(zero_phase)) job%zero_phase = zero_phase
        if (present(taps)) job%taps_real = taps
        if (present(sos)) job%sos_real = sos
        if (present(b)) job%numerator = b
        if (present(a)) job%denominator = a
        call apply_cc(x,y,axis,-1,job,st)
        if (present(status)) status = st
    end subroutine decimate_into_complex
    subroutine decimate_into_real_complex_tf(x, y, q, b, a, zero_phase, axis, status)
        real(dp), intent(in) :: x(..) !! real input array, rank 1-15
        complex(dp), intent(inout) :: y(..) !! caller-allocated complex output
        integer, intent(in) :: q !! positive decimation factor
        complex(dp), intent(in) :: b(:) !! numerator coefficients, ascending lag
        complex(dp), intent(in) :: a(:) !! denominator coefficients, nonzero leading element
        logical, intent(in), optional :: zero_phase !! default true
        integer, intent(in), optional :: axis !! default last axis
        integer, intent(out), optional :: status !! zero on success
        type(multirate_job) :: job
        integer :: st
        job%method = 4
        job%down = q
        job%numerator = b
        job%denominator = a
        if (present(zero_phase)) job%zero_phase = zero_phase
        call apply_rc(x,y,axis,-1,job,st)
        if (present(status)) status = st
    end subroutine decimate_into_real_complex_tf

    subroutine decimate_into_complex_sos(x, y, q, sos, zero_phase, axis, status)
        complex(dp), intent(in) :: x(..) !! input array, rank 1 through 15
        complex(dp), intent(inout) :: y(..) !! preallocated complex result, unchanged on failure
        integer, intent(in) :: q !! positive decimation factor
        complex(dp), intent(in) :: sos(:,:) !! complex SOS (nsections,6), a0=1
        logical, intent(in), optional :: zero_phase !! forward-backward filtering (default true)
        integer, intent(in), optional :: axis !! selected axis (default last)
        integer, intent(out), optional :: status !! zero on success
        type(multirate_job) :: job
        integer :: st

        job%method = 4
        job%down = q
        job%sos_complex = sos
        if (present(zero_phase)) job%zero_phase = zero_phase
        call apply_cc(x,y,axis,-1,job,st)
        if (present(status)) status = st
    end subroutine decimate_into_complex_sos

    subroutine decimate_into_real_complex_sos(x, y, q, sos, zero_phase, axis, status)
        real(dp), intent(in) :: x(..) !! real input array, rank 1 through 15
        complex(dp), intent(inout) :: y(..) !! preallocated complex result, unchanged on failure
        integer, intent(in) :: q !! positive decimation factor
        complex(dp), intent(in) :: sos(:,:) !! complex SOS (nsections,6), a0=1
        logical, intent(in), optional :: zero_phase !! forward-backward filtering (default true)
        integer, intent(in), optional :: axis !! selected axis (default last)
        integer, intent(out), optional :: status !! zero on success
        type(multirate_job) :: job
        integer :: st

        job%method = 4
        job%down = q
        job%sos_complex = sos
        if (present(zero_phase)) job%zero_phase = zero_phase
        call apply_rc(x,y,axis,-1,job,st)
        if (present(status)) status = st
    end subroutine decimate_into_real_complex_sos

end module scifort_signal_multirate_nd_into
