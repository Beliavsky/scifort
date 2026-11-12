! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Arbitrary-rank FFT analytic signals on selected axes.  Scalar FFT lanes
! use SciFort's reusable transform plans.  All output is staged before commit.
module scifort_signal_hilbert_nd
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_kinds, only : dp
    use scifort_signal_spectral, only : fft_plan, spectral_status_success
    use scifort_signal_hilbert, only : hilbert_status_success, hilbert_status_invalid_input
    implicit none
    private
    public :: hilbert_into, hilbert2_into

contains

    subroutine flatten_input(x, values)
        real(dp), intent(in) :: x(..) !! input of rank 1 to 15
        real(dp), allocatable, intent(out) :: values(:) !! flat column-major samples
        select rank(x)
        rank (1)
            values = reshape(x, [size(x)])
        rank (2)
            values = reshape(x, [size(x)])
        rank (3)
            values = reshape(x, [size(x)])
        rank (4)
            values = reshape(x, [size(x)])
        rank (5)
            values = reshape(x, [size(x)])
        rank (6)
            values = reshape(x, [size(x)])
        rank (7)
            values = reshape(x, [size(x)])
        rank (8)
            values = reshape(x, [size(x)])
        rank (9)
            values = reshape(x, [size(x)])
        rank (10)
            values = reshape(x, [size(x)])
        rank (11)
            values = reshape(x, [size(x)])
        rank (12)
            values = reshape(x, [size(x)])
        rank (13)
            values = reshape(x, [size(x)])
        rank (14)
            values = reshape(x, [size(x)])
        rank (15)
            values = reshape(x, [size(x)])
        rank default
            allocate(values(0))
        end select
    end subroutine flatten_input

    subroutine commit_output(y, values)
        complex(dp), intent(inout) :: y(..) !! validated output array of rank 1 to 15
        complex(dp), intent(in) :: values(:) !! complete column-major output samples
        select rank(y)
        rank (1)
            y = reshape(values, shape(y))
        rank (2)
            y = reshape(values, shape(y))
        rank (3)
            y = reshape(values, shape(y))
        rank (4)
            y = reshape(values, shape(y))
        rank (5)
            y = reshape(values, shape(y))
        rank (6)
            y = reshape(values, shape(y))
        rank (7)
            y = reshape(values, shape(y))
        rank (8)
            y = reshape(values, shape(y))
        rank (9)
            y = reshape(values, shape(y))
        rank (10)
            y = reshape(values, shape(y))
        rank (11)
            y = reshape(values, shape(y))
        rank (12)
            y = reshape(values, shape(y))
        rank (13)
            y = reshape(values, shape(y))
        rank (14)
            y = reshape(values, shape(y))
        rank (15)
            y = reshape(values, shape(y))
        rank default
            continue
        end select
    end subroutine commit_output

    subroutine validate_shapes(indims, outdims, axes, lengths, code)
        integer, intent(in) :: indims(:) !! original shape
        integer, intent(in) :: outdims(:) !! allocated destination shape
        integer, intent(in) :: axes(:) !! requested one-based or negative axes
        integer, intent(in) :: lengths(:) !! requested transform lengths
        integer, intent(out) :: code !! status, zero if destination matches
        integer :: d, j, ax

        code = hilbert_status_invalid_input
        if (size(indims) < size(axes) .or. size(indims) > 15) return
        if (size(indims) /= size(outdims) .or. size(axes) /= size(lengths)) return
        if (any(lengths < 1)) return
        do d = 1, size(indims)
            j = indims(d)
            do ax = 1, size(axes)
                if (d == axes(ax)) j = lengths(ax)
            end do
            if (outdims(d) /= j) return
        end do
        code = hilbert_status_success
    end subroutine validate_shapes

    subroutine make_strides(dims, strides)
        integer, intent(in) :: dims(:) !! tensor shape
        integer, intent(out) :: strides(:) !! Fortran column-major strides
        integer :: d
        strides(1) = 1
        do d = 2, size(dims)
            strides(d) = strides(d-1)*dims(d-1)
        end do
    end subroutine make_strides

    subroutine pad_or_truncate(x, indims, outdims, values)
        real(dp), intent(in) :: x(:) !! flat original signal
        integer, intent(in) :: indims(:) !! input shape
        integer, intent(in) :: outdims(:) !! destination shape
        complex(dp), intent(out) :: values(:) !! padded/truncated tensor as flat values
        integer :: ostride(size(outdims))
        integer :: coords(size(indims)), i, d, q, dest
        logical :: inside

        values = (0.0_dp, 0.0_dp)
        call make_strides(outdims, ostride)
        do i = 1, size(x)
            q = i-1
            inside = .true.
            dest = 1
            do d = 1, size(indims)
                coords(d) = mod(q,indims(d)) + 1
                q = q / indims(d)
                if (coords(d) > outdims(d)) then
                    inside = .false.
                    exit
                end if
                dest = dest + (coords(d)-1)*ostride(d)
            end do
            if (inside) values(dest) = cmplx(x(i),0.0_dp,dp)
        end do
    end subroutine pad_or_truncate

    subroutine transform_axis(values, dims, axis, plan, forward, code)
        complex(dp), intent(inout) :: values(:) !! flat spectrum/samples
        integer, intent(in) :: dims(:) !! tensor shape
        integer, intent(in) :: axis !! positive one-based transform dimension
        type(fft_plan), intent(inout) :: plan !! reusable FFT plan for this axis
        logical, intent(in) :: forward !! forward FFT if true, inverse otherwise
        integer, intent(out) :: code !! zero after all lanes succeed
        complex(dp), allocatable :: line(:)
        integer :: strides(size(dims)), n, lanes, lane, j, d, q, base

        code = hilbert_status_success
        if (size(values) == 0) return
        call make_strides(dims,strides)
        n = dims(axis)
        lanes = size(values)/n
        allocate(line(n))
        do lane = 0, lanes-1
            q = lane
            base = 1
            do d = 1, size(dims)
                if (d == axis) cycle
                base = base + mod(q,dims(d))*strides(d)
                q = q/dims(d)
            end do
            do j = 1,n
                line(j) = values(base+(j-1)*strides(axis))
            end do
            if (forward) then
                call plan%execute(line,code)
            else
                call plan%inverse(line,code)
            end if
            if (code /= spectral_status_success) return
            do j = 1,n
                values(base+(j-1)*strides(axis)) = line(j)
            end do
        end do
    end subroutine transform_axis

    subroutine apply_positive_mask(values, dims, axis, two_dimensional)
        complex(dp), intent(inout) :: values(:) !! flat spectrum to multiply in place
        integer, intent(in) :: dims(:) !! tensor dimensions
        integer, intent(in) :: axis !! selected transform dimension
        logical, intent(in) :: two_dimensional !! use SciPy hilbert2 Nyquist convention
        integer :: n, i, k, stride
        real(dp), allocatable :: weights(:)

        n = dims(axis)
        allocate(weights(n))
        weights = 0.0_dp
        weights(1) = 1.0_dp
        if (mod(n,2) == 0) then
            weights(2:n/2) = 2.0_dp
            if (.not. two_dimensional) weights(n/2+1) = 1.0_dp
        else
            weights(2:(n+1)/2) = 2.0_dp
        end if
        stride = 1
        do k = 1,axis-1
            stride = stride*dims(k)
        end do
        do i = 1,size(values)
            k = mod((i-1)/stride,n)+1
            values(i) = values(i)*weights(k)
        end do
    end subroutine apply_positive_mask

    subroutine analytic_transform(x, y, axes, code)
        real(dp), intent(in) :: x(..) !! input real tensor, at most rank 15
        complex(dp), intent(inout) :: y(..) !! destination, not touched until success
        integer, intent(in) :: axes(:) !! distinct, normalized axes
        integer, intent(out) :: code !! zero on successful transforms
        integer, allocatable :: indims(:), outdims(:)
        real(dp), allocatable :: flat_in(:)
        complex(dp), allocatable :: work(:)
        type(fft_plan) :: plans(size(axes))
        integer :: d, ierr

        code = hilbert_status_invalid_input
        indims = shape(x)
        outdims = shape(y)
        if (size(indims) /= size(outdims)) return
        call validate_shapes(indims,outdims,axes,outdims(axes),code)
        if (code /= hilbert_status_success) return
        call flatten_input(x,flat_in)
        if (size(flat_in) /= size(x)) then
            code = hilbert_status_invalid_input
            return
        end if
        if (.not. all(ieee_is_finite(flat_in))) then
            code = hilbert_status_invalid_input
            return
        end if
        allocate(work(size(y)))
        call pad_or_truncate(flat_in,indims,outdims,work)
        if (size(work) > 0) then
            do d = 1,size(axes)
                call plans(d)%initialize(outdims(axes(d)),ierr)
                if (ierr /= spectral_status_success) then
                    code = hilbert_status_invalid_input
                    return
                end if
            end do
            do d = 1,size(axes)
                call transform_axis(work,outdims,axes(d),plans(d),.true.,ierr)
                if (ierr /= spectral_status_success) then
                    code = hilbert_status_invalid_input
                    return
                end if
            end do
            do d = 1,size(axes)
                call apply_positive_mask(work,outdims,axes(d),size(axes) == 2)
            end do
            do d = size(axes),1,-1
                call transform_axis(work,outdims,axes(d),plans(d),.false.,ierr)
                if (ierr /= spectral_status_success) then
                    code = hilbert_status_invalid_input
                    return
                end if
            end do
            if (.not. all(ieee_is_finite(real(work,dp))) .or. &
                .not. all(ieee_is_finite(aimag(work)))) then
                code = hilbert_status_invalid_input
                return
            end if
        end if
        call commit_output(y,work)
        code = hilbert_status_success
    end subroutine analytic_transform

    subroutine hilbert_into(x, y, axis, n, status)
        real(dp), intent(in) :: x(..) !! real samples, rank 1 to 15
        complex(dp), intent(inout) :: y(..) !! analytic signal; same batch dimensions
        integer, intent(in), optional :: axis !! one-based or negative axis, default last
        integer, intent(in), optional :: n !! transform length, default input axis length
        integer, intent(out), optional :: status !! zero if destination updated
        integer, allocatable :: dims(:), odims(:)
        integer :: ax, code, length

        code = hilbert_status_invalid_input
        dims = shape(x)
        odims = shape(y)
        if (size(dims) < 1 .or. size(dims) > 15 .or. &
            size(odims) /= size(dims)) goto 900
        ax = size(dims)
        if (present(axis)) ax = axis
        if (ax < 0) ax = ax+size(dims)+1
        if (ax < 1 .or. ax > size(dims)) goto 900
        length = dims(ax)
        if (present(n)) length = n
        call validate_shapes(dims,odims,[ax],[length],code)
        if (code /= hilbert_status_success) goto 900
        call analytic_transform(x,y,[ax],code)
900     if (present(status)) status = code
    end subroutine hilbert_into

    subroutine hilbert2_into(x, y, axes, n, status)
        real(dp), intent(in) :: x(..) !! real samples, rank 2 to 15
        complex(dp), intent(inout) :: y(..) !! two-axis analytic signal
        integer, intent(in), optional :: axes(:) !! two one-based or negative axes, default last two
        integer, intent(in), optional :: n(:) !! two positive transform lengths
        integer, intent(out), optional :: status !! zero if destination updated
        integer, allocatable :: dims(:), odims(:)
        integer :: ax(2), lengths(2), code, d

        code = hilbert_status_invalid_input
        dims = shape(x)
        odims = shape(y)
        if (size(dims) < 2 .or. size(dims) > 15 .or. &
            size(odims) /= size(dims)) goto 900
        ax = [size(dims)-1,size(dims)]
        if (present(axes)) then
            if (size(axes) /= 2) goto 900
            ax = axes
        end if
        do d = 1,2
            if (ax(d) < 0) ax(d) = ax(d)+size(dims)+1
        end do
        if (any(ax < 1) .or. any(ax > size(dims))) goto 900
        if (ax(1) == ax(2)) goto 900
        lengths = dims(ax)
        if (present(n)) then
            if (size(n) /= 2) goto 900
            lengths = n
        end if
        call validate_shapes(dims,odims,ax,lengths,code)
        if (code /= hilbert_status_success) goto 900
        call analytic_transform(x,y,ax,code)
900     if (present(status)) status = code
    end subroutine hilbert2_into

end module scifort_signal_hilbert_nd
