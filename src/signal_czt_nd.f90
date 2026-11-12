! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Arbitrary-rank chirp Z-transform and zoom FFT adapters.  Calls share one
! immutable CZT/zoom plan across all lanes along the chosen signal axis.
! Output is staged and committed only after every lane has succeeded.
module scifort_signal_czt_nd
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_kinds, only : dp
    use scifort_signal_czt, only : czt_plan, zoom_fft_plan, &
        czt_status_success, czt_status_invalid_input
    implicit none
    private

    interface czt_into
        module procedure czt_into_real, czt_into_complex
    end interface
    interface czt_plan_into
        module procedure czt_plan_into_real, czt_plan_into_complex
    end interface
    interface zoom_fft_into
        module procedure zoom_into_real_pair, zoom_into_complex_pair
        module procedure zoom_into_real_scalar, zoom_into_complex_scalar
    end interface
    interface zoom_fft_plan_into
        module procedure zoom_plan_into_real, zoom_plan_into_complex
    end interface
    public :: czt_into, czt_plan_into, zoom_fft_into, zoom_fft_plan_into

contains

    subroutine check_dimensions(in_dims, out_dims, axis, m, ax, status)
        integer, intent(in) :: in_dims(:) !! input shape, Fortran ordering
        integer, intent(in) :: out_dims(:) !! output shape, Fortran ordering
        integer, intent(in), optional :: axis !! 1-based or negative signal axis, default last
        integer, intent(in), optional :: m !! requested number of transformed samples
        integer, intent(out) :: ax !! normalized positive signal axis
        integer, intent(out) :: status !! zero when all shape constraints hold
        integer :: j, n

        status = czt_status_invalid_input
        n = size(in_dims)
        ax = 0
        if (n < 1 .or. n > 15 .or. size(out_dims) /= n) return
        ax = n
        if (present(axis)) ax = axis
        if (ax < 0) ax = ax+n+1
        if (ax < 1 .or. ax > n) return
        if (in_dims(ax) < 1 .or. out_dims(ax) < 1) return
        if (present(m)) then
            if (m /= out_dims(ax)) return
        else
            if (in_dims(ax) /= out_dims(ax)) return
        end if
        do j = 1,n
            if (j == ax) cycle
            if (in_dims(j) /= out_dims(j)) return
        end do
        status = czt_status_success
    end subroutine check_dimensions

    subroutine flatten_real(x, v)
        real(dp), intent(in) :: x(..) !! input array, rank 1 through 15
        complex(dp), allocatable, intent(out) :: v(:) !! contiguous Fortran-order copy
        select rank(x)
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
        rank default
            allocate(v(0))
        end select
    end subroutine flatten_real

    subroutine flatten_complex(x, v)
        complex(dp), intent(in) :: x(..) !! input array, rank 1 through 15
        complex(dp), allocatable, intent(out) :: v(:) !! contiguous Fortran-order copy
        select rank(x)
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
        rank default
            allocate(v(0))
        end select
    end subroutine flatten_complex

    subroutine scatter_complex(y, v)
        complex(dp), intent(inout) :: y(..) !! output written only after success
        complex(dp), intent(in) :: v(:) !! fully computed flat output
        select rank(y)
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
        rank default
            error stop 'unreachable rank in scatter_complex'
        end select
    end subroutine scatter_complex

    subroutine transform_lanes(src, dst, in_dims, ax, status, cplan, zplan)
        complex(dp), intent(in) :: src(:) !! flat Fortran-order input
        complex(dp), intent(out) :: dst(:) !! flat staged output
        integer, intent(in) :: in_dims(:) !! input shape, rank 1 through 15
        integer, intent(in) :: ax !! positive transform axis
        integer, intent(out) :: status !! status of first unsuccessful lane
        type(czt_plan), intent(in), optional :: cplan !! optional prepared CZT plan
        type(zoom_fft_plan), intent(in), optional :: zplan !! optional prepared zoom plan
        integer(int64) :: inner, outer, nx, ny, channel, batch, j, ii, io
        integer :: st
        complex(dp), allocatable :: line(:), result(:)

        status = czt_status_invalid_input
        if (present(cplan) .eqv. present(zplan)) return
        nx = int(in_dims(ax),int64)
        inner = product(int(in_dims(:ax-1),int64))
        outer = product(int(in_dims(ax+1:),int64))
        if (nx < 1 .or. inner < 0 .or. outer < 0) return
        if (size(src,kind=int64) == 0_int64 .and. size(dst,kind=int64) == 0_int64) then
            status = czt_status_success
            return
        end if
        if (inner < 1 .or. outer < 1) return
        if (size(src,kind=int64) /= inner*nx*outer) return
        if (mod(size(dst,kind=int64),inner*outer) /= 0_int64) return
        ny = size(dst,kind=int64)/(inner*outer)
        allocate(line(nx))
        do batch = 0_int64,outer-1_int64
            do channel = 0_int64,inner-1_int64
                do j = 0_int64,nx-1_int64
                    ii = batch*inner*nx + j*inner + channel + 1_int64
                    line(j+1) = src(ii)
                end do
                if (present(cplan)) then
                    result = cplan%transform(line,st)
                else
                    result = zplan%transform(line,st)
                end if
                if (st /= czt_status_success) then
                    status = st
                    return
                end if
                if (size(result,kind=int64) /= ny) return
                do j = 0_int64,ny-1_int64
                    io = batch*inner*ny + j*inner + channel + 1_int64
                    dst(io) = result(j+1)
                end do
            end do
        end do
        status = czt_status_success
    end subroutine transform_lanes

    subroutine czt_plan_into_complex(plan, x, y, axis, status)
        type(czt_plan), intent(in) :: plan !! reusable, initialized chirp-Z plan
        complex(dp), intent(in) :: x(..) !! input signals, rank 1 through 15
        complex(dp), intent(inout) :: y(..) !! caller-provided complex output
        integer, intent(in), optional :: axis !! 1-based or negative signal axis, default last
        integer, intent(out), optional :: status !! success or first lane failure
        complex(dp), allocatable :: src(:), dst(:)
        integer :: ax, st

        if (present(status)) status = czt_status_invalid_input
        call check_dimensions(shape(x),shape(y),axis,plan%output_length(),ax,st)
        if (st /= czt_status_success .or. .not. plan%is_valid()) return
        if (size(x,dim=ax) /= plan%input_length()) return
        call flatten_complex(x,src)
        allocate(dst(size(y)))
        call transform_lanes(src,dst,shape(x),ax,st,cplan=plan)
        if (st == czt_status_success) call scatter_complex(y,dst)
        if (present(status)) status = st
    end subroutine czt_plan_into_complex

    subroutine czt_plan_into_real(plan, x, y, axis, status)
        type(czt_plan), intent(in) :: plan !! reusable, initialized chirp-Z plan
        real(dp), intent(in) :: x(..) !! real input signals, rank 1 through 15
        complex(dp), intent(inout) :: y(..) !! caller-provided complex output
        integer, intent(in), optional :: axis !! 1-based or negative signal axis, default last
        integer, intent(out), optional :: status !! success or first lane failure
        complex(dp), allocatable :: src(:), dst(:)
        integer :: ax, st

        if (present(status)) status = czt_status_invalid_input
        call check_dimensions(shape(x),shape(y),axis,plan%output_length(),ax,st)
        if (st /= czt_status_success .or. .not. plan%is_valid()) return
        if (size(x,dim=ax) /= plan%input_length()) return
        call flatten_real(x,src)
        allocate(dst(size(y)))
        call transform_lanes(src,dst,shape(x),ax,st,cplan=plan)
        if (st == czt_status_success) call scatter_complex(y,dst)
        if (present(status)) status = st
    end subroutine czt_plan_into_real

    subroutine zoom_plan_into_complex(plan, x, y, axis, status)
        type(zoom_fft_plan), intent(in) :: plan !! reusable, initialized zoom FFT plan
        complex(dp), intent(in) :: x(..) !! complex signals, rank 1 through 15
        complex(dp), intent(inout) :: y(..) !! caller-provided complex output
        integer, intent(in), optional :: axis !! 1-based or negative signal axis, default last
        integer, intent(out), optional :: status !! success or first lane failure
        complex(dp), allocatable :: src(:), dst(:)
        integer :: ax, st

        if (present(status)) status = czt_status_invalid_input
        call check_dimensions(shape(x),shape(y),axis,plan%output_length(),ax,st)
        if (st /= czt_status_success .or. .not. plan%is_valid()) return
        if (size(x,dim=ax) /= plan%input_length()) return
        call flatten_complex(x,src)
        allocate(dst(size(y)))
        call transform_lanes(src,dst,shape(x),ax,st,zplan=plan)
        if (st == czt_status_success) call scatter_complex(y,dst)
        if (present(status)) status = st
    end subroutine zoom_plan_into_complex

    subroutine zoom_plan_into_real(plan, x, y, axis, status)
        type(zoom_fft_plan), intent(in) :: plan !! reusable, initialized zoom FFT plan
        real(dp), intent(in) :: x(..) !! real signals, rank 1 through 15
        complex(dp), intent(inout) :: y(..) !! caller-provided complex output
        integer, intent(in), optional :: axis !! 1-based or negative signal axis, default last
        integer, intent(out), optional :: status !! success or first lane failure
        complex(dp), allocatable :: src(:), dst(:)
        integer :: ax, st

        if (present(status)) status = czt_status_invalid_input
        call check_dimensions(shape(x),shape(y),axis,plan%output_length(),ax,st)
        if (st /= czt_status_success .or. .not. plan%is_valid()) return
        if (size(x,dim=ax) /= plan%input_length()) return
        call flatten_real(x,src)
        allocate(dst(size(y)))
        call transform_lanes(src,dst,shape(x),ax,st,zplan=plan)
        if (st == czt_status_success) call scatter_complex(y,dst)
        if (present(status)) status = st
    end subroutine zoom_plan_into_real

    subroutine czt_into_real(x, y, m, w, a, axis, status)
        real(dp), intent(in) :: x(..) !! input signals, rank 1 through 15
        complex(dp), intent(inout) :: y(..) !! complex output, unchanged on error
        integer, intent(in), optional :: m !! output length; default input-axis length
        complex(dp), intent(in), optional :: w !! CZT spiral ratio, default unit circle
        complex(dp), intent(in), optional :: a !! first spiral point, default one
        integer, intent(in), optional :: axis !! one-based or negative axis, default last
        integer, intent(out), optional :: status !! success or error code
        type(czt_plan) :: plan
        integer :: ax, st

        if (present(status)) status = czt_status_invalid_input
        call check_dimensions(shape(x),shape(y),axis,m,ax,st)
        if (st /= czt_status_success) return
        call plan%initialize(size(x,dim=ax),m=size(y,dim=ax),w=w,a=a,status=st)
        if (st /= czt_status_success) then
            if (present(status)) status = st
            return
        end if
        call czt_plan_into(plan,x,y,axis=axis,status=st)
        if (present(status)) status = st
    end subroutine czt_into_real

    subroutine czt_into_complex(x, y, m, w, a, axis, status)
        complex(dp), intent(in) :: x(..) !! input signals, rank 1 through 15
        complex(dp), intent(inout) :: y(..) !! complex output, unchanged on error
        integer, intent(in), optional :: m !! output length; default input-axis length
        complex(dp), intent(in), optional :: w !! CZT spiral ratio, default unit circle
        complex(dp), intent(in), optional :: a !! first spiral point, default one
        integer, intent(in), optional :: axis !! one-based or negative axis, default last
        integer, intent(out), optional :: status !! success or error code
        type(czt_plan) :: plan
        integer :: ax, st

        if (present(status)) status = czt_status_invalid_input
        call check_dimensions(shape(x),shape(y),axis,m,ax,st)
        if (st /= czt_status_success) return
        call plan%initialize(size(x,dim=ax),m=size(y,dim=ax),w=w,a=a,status=st)
        if (st /= czt_status_success) then
            if (present(status)) status = st
            return
        end if
        call czt_plan_into(plan,x,y,axis=axis,status=st)
        if (present(status)) status = st
    end subroutine czt_into_complex

    subroutine zoom_into_real_pair(x, y, fn, m, fs, endpoint, axis, status)
        real(dp), intent(in) :: x(..) !! input signals, rank 1 through 15
        complex(dp), intent(inout) :: y(..) !! complex output, unchanged on error
        real(dp), intent(in) :: fn(2) !! first and upper frequencies
        integer, intent(in), optional :: m !! output length; default input-axis length
        real(dp), intent(in), optional :: fs !! sampling rate, default 2
        logical, intent(in), optional :: endpoint !! include upper endpoint, default false
        integer, intent(in), optional :: axis !! one-based or negative axis, default last
        integer, intent(out), optional :: status !! success or error code
        type(zoom_fft_plan) :: plan
        integer :: ax, st

        if (present(status)) status = czt_status_invalid_input
        call check_dimensions(shape(x),shape(y),axis,m,ax,st)
        if (st /= czt_status_success) return
        call plan%initialize(size(x,dim=ax),fn,m=size(y,dim=ax),fs=fs,endpoint=endpoint,status=st)
        if (st /= czt_status_success) then
            if (present(status)) status = st
            return
        end if
        call zoom_fft_plan_into(plan,x,y,axis=axis,status=st)
        if (present(status)) status = st
    end subroutine zoom_into_real_pair

    subroutine zoom_into_real_scalar(x, y, fn, m, fs, endpoint, axis, status)
        real(dp), intent(in) :: x(..) !! input signals, rank 1 through 15
        complex(dp), intent(inout) :: y(..) !! complex output, unchanged on error
        real(dp), intent(in) :: fn !! upper frequency, starting at zero
        integer, intent(in), optional :: m !! output length; default input-axis length
        real(dp), intent(in), optional :: fs !! sampling rate, default 2
        logical, intent(in), optional :: endpoint !! include upper endpoint, default false
        integer, intent(in), optional :: axis !! one-based or negative axis, default last
        integer, intent(out), optional :: status !! success or error code
        type(zoom_fft_plan) :: plan
        integer :: ax, st

        if (present(status)) status = czt_status_invalid_input
        call check_dimensions(shape(x),shape(y),axis,m,ax,st)
        if (st /= czt_status_success) return
        call plan%initialize(size(x,dim=ax),fn,m=size(y,dim=ax),fs=fs,endpoint=endpoint,status=st)
        if (st /= czt_status_success) then
            if (present(status)) status = st
            return
        end if
        call zoom_fft_plan_into(plan,x,y,axis=axis,status=st)
        if (present(status)) status = st
    end subroutine zoom_into_real_scalar

    subroutine zoom_into_complex_pair(x, y, fn, m, fs, endpoint, axis, status)
        complex(dp), intent(in) :: x(..) !! input signals, rank 1 through 15
        complex(dp), intent(inout) :: y(..) !! complex output, unchanged on error
        real(dp), intent(in) :: fn(2) !! first and upper frequencies
        integer, intent(in), optional :: m !! output length; default input-axis length
        real(dp), intent(in), optional :: fs !! sampling rate, default 2
        logical, intent(in), optional :: endpoint !! include upper endpoint, default false
        integer, intent(in), optional :: axis !! one-based or negative axis, default last
        integer, intent(out), optional :: status !! success or error code
        type(zoom_fft_plan) :: plan
        integer :: ax, st

        if (present(status)) status = czt_status_invalid_input
        call check_dimensions(shape(x),shape(y),axis,m,ax,st)
        if (st /= czt_status_success) return
        call plan%initialize(size(x,dim=ax),fn,m=size(y,dim=ax),fs=fs,endpoint=endpoint,status=st)
        if (st /= czt_status_success) then
            if (present(status)) status = st
            return
        end if
        call zoom_fft_plan_into(plan,x,y,axis=axis,status=st)
        if (present(status)) status = st
    end subroutine zoom_into_complex_pair

    subroutine zoom_into_complex_scalar(x, y, fn, m, fs, endpoint, axis, status)
        complex(dp), intent(in) :: x(..) !! input signals, rank 1 through 15
        complex(dp), intent(inout) :: y(..) !! complex output, unchanged on error
        real(dp), intent(in) :: fn !! upper frequency, starting at zero
        integer, intent(in), optional :: m !! output length; default input-axis length
        real(dp), intent(in), optional :: fs !! sampling rate, default 2
        logical, intent(in), optional :: endpoint !! include upper endpoint, default false
        integer, intent(in), optional :: axis !! one-based or negative axis, default last
        integer, intent(out), optional :: status !! success or error code
        type(zoom_fft_plan) :: plan
        integer :: ax, st

        if (present(status)) status = czt_status_invalid_input
        call check_dimensions(shape(x),shape(y),axis,m,ax,st)
        if (st /= czt_status_success) return
        call plan%initialize(size(x,dim=ax),fn,m=size(y,dim=ax),fs=fs,endpoint=endpoint,status=st)
        if (st /= czt_status_success) then
            if (present(status)) status = st
            return
        end if
        call zoom_fft_plan_into(plan,x,y,axis=axis,status=st)
        if (present(status)) status = st
    end subroutine zoom_into_complex_scalar

end module scifort_signal_czt_nd
