! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Cascaded second-order-section filters; ranks 1 through 15.
! The SOS matrix is (number_of_sections, 6), with a0 = 1 in column 4.
module scifort_signal_sos
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
    use, intrinsic :: iso_c_binding, only : c_loc, c_f_pointer
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_kinds, only : dp
    use scifort_signal_lfilter, only : lfilter_zi, lfilter_status_success, &
        lfilter_status_no_steady_state
    implicit none
    private

    integer, parameter, public :: sos_status_success = 0
    integer, parameter, public :: sos_status_invalid_input = 1
    integer, parameter, public :: sos_status_no_steady_state = 2

    interface sosfilt
        module procedure sosfilt_real
        module procedure sosfilt_complex
        module procedure sosfilt_matrix_real, sosfilt_matrix_complex
        module procedure sosfilt_rank3_real, sosfilt_rank3_complex
    end interface
    interface sosfilt_zi
        module procedure sosfilt_zi_real
        module procedure sosfilt_zi_complex
    end interface
    interface sosfilt_nd
        module procedure sosfilt_nd_real, sosfilt_nd_complex
        module procedure sosfilt_nd_real_sos_complex_input, sosfilt_nd_complex_sos_real_input
    end interface sosfilt_nd
    public :: sosfilt, sosfilt_nd, sosfilt_zi

contains

    pure logical function sos_real_valid(sos) result(ok)
        real(dp), intent(in) :: sos(:,:) !! sections with six coefficients per row
        ok = .false.
        if (size(sos, 1) < 1 .or. size(sos, 2) /= 6) return
        if (.not. all(ieee_is_finite(sos))) return
        if (any(sos(:,4) /= 1.0_dp)) return
        ok = .true.
    end function

    pure logical function sos_complex_valid(sos) result(ok)
        complex(dp), intent(in) :: sos(:,:) !! sections with six coefficients per row
        ok = .false.
        if (size(sos, 1) < 1 .or. size(sos, 2) /= 6) return
        if (.not. all(ieee_is_finite(real(sos, dp)))) return
        if (.not. all(ieee_is_finite(aimag(sos)))) return
        if (any(sos(:,4) /= (1.0_dp, 0.0_dp))) return
        ok = .true.
    end function

    function sosfilt_real(sos, x, zi, zf, status) result(y)
        real(dp), intent(in) :: sos(:,:) !! SOS coefficients (nsections,6); a0 must equal one
        real(dp), intent(in) :: x(:) !! input samples
        real(dp), intent(in), optional :: zi(:,:) !! initial delays (nsections,2)
        real(dp), allocatable, intent(out), optional :: zf(:,:) !! final delays (nsections,2)
        integer, intent(out), optional :: status !! success or invalid input
        real(dp), allocatable :: y(:)
        real(dp), allocatable :: state(:,:)
        real(dp) :: sample, tmp
        integer :: i, j, n

        allocate(y(size(x)))
        y = ieee_value(0.0_dp, ieee_quiet_nan)
        if (present(zf)) allocate(zf(0,0))
        if (present(status)) status = sos_status_invalid_input
        if (.not. sos_real_valid(sos)) return
        n = size(sos,1)
        if (present(zi)) then
            if (size(zi,1) /= n .or. size(zi,2) /= 2) return
            if (.not. all(ieee_is_finite(zi))) return
        end if
        allocate(state(n,2))
        state = 0.0_dp
        if (present(zi)) state = zi
        do i = 1, size(x)
            sample = x(i)
            do j = 1, n
                tmp = sos(j,1) * sample + state(j,1)
                state(j,1) = sos(j,2) * sample - sos(j,5) * tmp + state(j,2)
                state(j,2) = sos(j,3) * sample - sos(j,6) * tmp
                sample = tmp
            end do
            y(i) = sample
        end do
        if (present(zf)) zf = state
        if (present(status)) status = sos_status_success
    end function

    function sosfilt_complex(sos, x, zi, zf, status) result(y)
        complex(dp), intent(in) :: sos(:,:) !! SOS coefficients (nsections,6); a0 must equal one
        complex(dp), intent(in) :: x(:) !! input samples
        complex(dp), intent(in), optional :: zi(:,:) !! initial delays (nsections,2)
        complex(dp), allocatable, intent(out), optional :: zf(:,:) !! final delays (nsections,2)
        integer, intent(out), optional :: status !! success or invalid input
        complex(dp), allocatable :: y(:)
        complex(dp), allocatable :: state(:,:)
        complex(dp) :: sample, tmp
        real(dp) :: nanval
        integer :: i, j, n

        nanval = ieee_value(0.0_dp, ieee_quiet_nan)
        allocate(y(size(x)))
        y = cmplx(nanval, nanval, dp)
        if (present(zf)) allocate(zf(0,0))
        if (present(status)) status = sos_status_invalid_input
        if (.not. sos_complex_valid(sos)) return
        n = size(sos,1)
        if (present(zi)) then
            if (size(zi,1) /= n .or. size(zi,2) /= 2) return
            if (.not. all(ieee_is_finite(real(zi, dp)))) return
            if (.not. all(ieee_is_finite(aimag(zi)))) return
        end if
        allocate(state(n,2))
        state = (0.0_dp, 0.0_dp)
        if (present(zi)) state = zi
        do i = 1, size(x)
            sample = x(i)
            do j = 1, n
                tmp = sos(j,1) * sample + state(j,1)
                state(j,1) = sos(j,2) * sample - sos(j,5) * tmp + state(j,2)
                state(j,2) = sos(j,3) * sample - sos(j,6) * tmp
                sample = tmp
            end do
            y(i) = sample
        end do
        if (present(zf)) zf = state
        if (present(status)) status = sos_status_success
    end function

    function sosfilt_zi_real(sos, status) result(zi)
        real(dp), intent(in) :: sos(:,:) !! section coefficients (nsections,6)
        integer, intent(out), optional :: status !! success, invalid, or no steady state
        real(dp), allocatable :: zi(:,:)
        real(dp), allocatable :: state(:)
        real(dp) :: scale, denominator
        integer :: i, ls

        allocate(zi(size(sos,1), 2))
        zi = ieee_value(0.0_dp, ieee_quiet_nan)
        if (present(status)) status = sos_status_invalid_input
        if (.not. sos_real_valid(sos)) return
        scale = 1.0_dp
        do i = 1, size(sos,1)
            state = lfilter_zi(sos(i,1:3), sos(i,4:6), status=ls)
            if (ls /= lfilter_status_success) then
                if (present(status)) then
                    if (ls == lfilter_status_no_steady_state) status = sos_status_no_steady_state
                end if
                return
            end if
            denominator = sum(sos(i,4:6))
            zi(i,:) = scale * state
            scale = scale * (sum(sos(i,1:3)) / denominator)
            if (.not. ieee_is_finite(scale)) return
        end do
        if (present(status)) status = sos_status_success
    end function

    function sosfilt_zi_complex(sos, status) result(zi)
        complex(dp), intent(in) :: sos(:,:) !! section coefficients (nsections,6)
        integer, intent(out), optional :: status !! success, invalid, or no steady state
        complex(dp), allocatable :: zi(:,:)
        complex(dp), allocatable :: state(:)
        complex(dp) :: scale, denominator
        real(dp) :: nanval
        integer :: i, ls

        nanval = ieee_value(0.0_dp, ieee_quiet_nan)
        allocate(zi(size(sos,1), 2))
        zi = cmplx(nanval, nanval, dp)
        if (present(status)) status = sos_status_invalid_input
        if (.not. sos_complex_valid(sos)) return
        scale = (1.0_dp, 0.0_dp)
        do i = 1, size(sos,1)
            state = lfilter_zi(sos(i,1:3), sos(i,4:6), status=ls)
            if (ls /= lfilter_status_success) then
                if (present(status)) then
                    if (ls == lfilter_status_no_steady_state) status = sos_status_no_steady_state
                end if
                return
            end if
            denominator = sum(sos(i,4:6))
            zi(i,:) = scale * state
            scale = scale * (sum(sos(i,1:3)) / denominator)
            if (.not. (ieee_is_finite(real(scale, dp)) .and. ieee_is_finite(aimag(scale)))) return
        end do
        if (present(status)) status = sos_status_success
    end function

    function sosfilt_matrix_real(sos, x, zi, zf, status, axis) result(y)
        real(dp), intent(in) :: sos(:,:) !! coefficients, one section per row
        real(dp), intent(in) :: x(:,:) !! input data
        real(dp), intent(in), optional :: zi(:,:,:) !! initial SOS delays
        real(dp), allocatable, intent(out), optional :: zf(:,:,:) !! final SOS delays
        integer, intent(out), optional :: status !! success or invalid input
        integer, intent(in), optional :: axis !! one-based axis, default 2
        real(dp), allocatable :: y(:,:)
        real(dp), allocatable :: line(:)
        real(dp) :: nanval
        integer :: along, i, j, code
        real(dp), allocatable :: last_state(:,:)
        integer :: state_shape(3)

        nanval = ieee_value(0.0_dp, ieee_quiet_nan)
        allocate(y(size(x,1),size(x,2)))
        y = nanval
        if (present(zf)) allocate(zf(0,0,0))
        if (present(status)) status = sos_status_invalid_input
        along = 2
        if (present(axis)) along = axis
        if (along < 1 .or. along > 2) return
        if (.not. sos_real_valid(sos)) return
        state_shape = [size(sos,1), size(x,1), size(x,2)]
        state_shape(along+1) = 2
        if (present(zi)) then
            if (any(shape(zi) /= state_shape)) return
            if (.not. all(ieee_is_finite(zi))) return
        end if
        if (present(zf)) then
            deallocate(zf)
            allocate(zf(state_shape(1),state_shape(2),state_shape(3)))
        end if
        select case (along)
        case (1)
            do i = 1, size(x,2)
                if (present(zi)) then
                    line = sosfilt_real(sos, x(:,i), zi=zi(:,:,i), &
                        zf=last_state, status=code)
                else
                    line = sosfilt_real(sos, x(:,i), zf=last_state, status=code)
                end if
                if (code /= sos_status_success) return
                y(:,i) = line
                if (present(zf)) zf(:,:,i) = last_state
            end do
        case (2)
            do i = 1, size(x,1)
                if (present(zi)) then
                    line = sosfilt_real(sos, x(i,:), zi=zi(:,i,:), &
                        zf=last_state, status=code)
                else
                    line = sosfilt_real(sos, x(i,:), zf=last_state, status=code)
                end if
                if (code /= sos_status_success) return
                y(i,:) = line
                if (present(zf)) zf(:,i,:) = last_state
            end do
        end select
        if (present(status)) status = sos_status_success
    end function sosfilt_matrix_real

    function sosfilt_matrix_complex(sos, x, zi, zf, status, axis) result(y)
        complex(dp), intent(in) :: sos(:,:) !! coefficients, one section per row
        complex(dp), intent(in) :: x(:,:) !! input data
        complex(dp), intent(in), optional :: zi(:,:,:) !! initial SOS delays
        complex(dp), allocatable, intent(out), optional :: zf(:,:,:) !! final SOS delays
        integer, intent(out), optional :: status !! success or invalid input
        integer, intent(in), optional :: axis !! one-based axis, default 2
        complex(dp), allocatable :: y(:,:)
        complex(dp), allocatable :: line(:)
        real(dp) :: nanval
        integer :: along, i, j, code
        complex(dp), allocatable :: last_state(:,:)
        integer :: state_shape(3)

        nanval = ieee_value(0.0_dp, ieee_quiet_nan)
        allocate(y(size(x,1),size(x,2)))
        y = cmplx(nanval,nanval,dp)
        if (present(zf)) allocate(zf(0,0,0))
        if (present(status)) status = sos_status_invalid_input
        along = 2
        if (present(axis)) along = axis
        if (along < 1 .or. along > 2) return
        if (.not. sos_complex_valid(sos)) return
        state_shape = [size(sos,1), size(x,1), size(x,2)]
        state_shape(along+1) = 2
        if (present(zi)) then
            if (any(shape(zi) /= state_shape)) return
            if (.not. all(ieee_is_finite(real(zi,dp)))) return
            if (.not. all(ieee_is_finite(aimag(zi)))) return
        end if
        if (present(zf)) then
            deallocate(zf)
            allocate(zf(state_shape(1),state_shape(2),state_shape(3)))
        end if
        select case (along)
        case (1)
            do i = 1, size(x,2)
                if (present(zi)) then
                    line = sosfilt_complex(sos, x(:,i), zi=zi(:,:,i), &
                        zf=last_state, status=code)
                else
                    line = sosfilt_complex(sos, x(:,i), zf=last_state, status=code)
                end if
                if (code /= sos_status_success) return
                y(:,i) = line
                if (present(zf)) zf(:,:,i) = last_state
            end do
        case (2)
            do i = 1, size(x,1)
                if (present(zi)) then
                    line = sosfilt_complex(sos, x(i,:), zi=zi(:,i,:), &
                        zf=last_state, status=code)
                else
                    line = sosfilt_complex(sos, x(i,:), zf=last_state, status=code)
                end if
                if (code /= sos_status_success) return
                y(i,:) = line
                if (present(zf)) zf(:,i,:) = last_state
            end do
        end select
        if (present(status)) status = sos_status_success
    end function sosfilt_matrix_complex

    function sosfilt_rank3_real(sos, x, zi, zf, status, axis) result(y)
        real(dp), intent(in) :: sos(:,:) !! coefficients, one section per row
        real(dp), intent(in) :: x(:,:,:) !! input data
        real(dp), intent(in), optional :: zi(:,:,:,:) !! initial SOS delays
        real(dp), allocatable, intent(out), optional :: zf(:,:,:,:) !! final SOS delays
        integer, intent(out), optional :: status !! success or invalid input
        integer, intent(in), optional :: axis !! one-based axis, default 3
        real(dp), allocatable :: y(:,:,:)
        real(dp), allocatable :: line(:)
        real(dp) :: nanval
        integer :: along, i, j, code
        real(dp), allocatable :: last_state(:,:)
        integer :: state_shape(4)

        nanval = ieee_value(0.0_dp, ieee_quiet_nan)
        allocate(y(size(x,1),size(x,2),size(x,3)))
        y = nanval
        if (present(zf)) allocate(zf(0,0,0,0))
        if (present(status)) status = sos_status_invalid_input
        along = 3
        if (present(axis)) along = axis
        if (along < 1 .or. along > 3) return
        if (.not. sos_real_valid(sos)) return
        state_shape = [size(sos,1), size(x,1), size(x,2), size(x,3)]
        state_shape(along+1) = 2
        if (present(zi)) then
            if (any(shape(zi) /= state_shape)) return
            if (.not. all(ieee_is_finite(zi))) return
        end if
        if (present(zf)) then
            deallocate(zf)
            allocate(zf(state_shape(1),state_shape(2),state_shape(3),state_shape(4)))
        end if
        select case (along)
        case (1)
            do j = 1, size(x,3)
                do i = 1, size(x,2)
                    if (present(zi)) then
                        line = sosfilt_real(sos, x(:,i,j), zi=zi(:,:,i,j), &
                            zf=last_state, status=code)
                    else
                        line = sosfilt_real(sos, x(:,i,j), zf=last_state, status=code)
                    end if
                    if (code /= sos_status_success) return
                    y(:,i,j) = line
                    if (present(zf)) zf(:,:,i,j) = last_state
                end do
            end do
        case (2)
            do j = 1, size(x,3)
                do i = 1, size(x,1)
                    if (present(zi)) then
                        line = sosfilt_real(sos, x(i,:,j), zi=zi(:,i,:,j), &
                            zf=last_state, status=code)
                    else
                        line = sosfilt_real(sos, x(i,:,j), zf=last_state, status=code)
                    end if
                    if (code /= sos_status_success) return
                    y(i,:,j) = line
                    if (present(zf)) zf(:,i,:,j) = last_state
                end do
            end do
        case (3)
            do j = 1, size(x,2)
                do i = 1, size(x,1)
                    if (present(zi)) then
                        line = sosfilt_real(sos, x(i,j,:), zi=zi(:,i,j,:), &
                            zf=last_state, status=code)
                    else
                        line = sosfilt_real(sos, x(i,j,:), zf=last_state, status=code)
                    end if
                    if (code /= sos_status_success) return
                    y(i,j,:) = line
                    if (present(zf)) zf(:,i,j,:) = last_state
                end do
            end do
        end select
        if (present(status)) status = sos_status_success
    end function sosfilt_rank3_real

    function sosfilt_rank3_complex(sos, x, zi, zf, status, axis) result(y)
        complex(dp), intent(in) :: sos(:,:) !! coefficients, one section per row
        complex(dp), intent(in) :: x(:,:,:) !! input data
        complex(dp), intent(in), optional :: zi(:,:,:,:) !! initial SOS delays
        complex(dp), allocatable, intent(out), optional :: zf(:,:,:,:) !! final SOS delays
        integer, intent(out), optional :: status !! success or invalid input
        integer, intent(in), optional :: axis !! one-based axis, default 3
        complex(dp), allocatable :: y(:,:,:)
        complex(dp), allocatable :: line(:)
        real(dp) :: nanval
        integer :: along, i, j, code
        complex(dp), allocatable :: last_state(:,:)
        integer :: state_shape(4)

        nanval = ieee_value(0.0_dp, ieee_quiet_nan)
        allocate(y(size(x,1),size(x,2),size(x,3)))
        y = cmplx(nanval,nanval,dp)
        if (present(zf)) allocate(zf(0,0,0,0))
        if (present(status)) status = sos_status_invalid_input
        along = 3
        if (present(axis)) along = axis
        if (along < 1 .or. along > 3) return
        if (.not. sos_complex_valid(sos)) return
        state_shape = [size(sos,1), size(x,1), size(x,2), size(x,3)]
        state_shape(along+1) = 2
        if (present(zi)) then
            if (any(shape(zi) /= state_shape)) return
            if (.not. all(ieee_is_finite(real(zi,dp)))) return
            if (.not. all(ieee_is_finite(aimag(zi)))) return
        end if
        if (present(zf)) then
            deallocate(zf)
            allocate(zf(state_shape(1),state_shape(2),state_shape(3),state_shape(4)))
        end if
        select case (along)
        case (1)
            do j = 1, size(x,3)
                do i = 1, size(x,2)
                    if (present(zi)) then
                        line = sosfilt_complex(sos, x(:,i,j), zi=zi(:,:,i,j), &
                            zf=last_state, status=code)
                    else
                        line = sosfilt_complex(sos, x(:,i,j), zf=last_state, status=code)
                    end if
                    if (code /= sos_status_success) return
                    y(:,i,j) = line
                    if (present(zf)) zf(:,:,i,j) = last_state
                end do
            end do
        case (2)
            do j = 1, size(x,3)
                do i = 1, size(x,1)
                    if (present(zi)) then
                        line = sosfilt_complex(sos, x(i,:,j), zi=zi(:,i,:,j), &
                            zf=last_state, status=code)
                    else
                        line = sosfilt_complex(sos, x(i,:,j), zf=last_state, status=code)
                    end if
                    if (code /= sos_status_success) return
                    y(i,:,j) = line
                    if (present(zf)) zf(:,i,:,j) = last_state
                end do
            end do
        case (3)
            do j = 1, size(x,2)
                do i = 1, size(x,1)
                    if (present(zi)) then
                        line = sosfilt_complex(sos, x(i,j,:), zi=zi(:,i,j,:), &
                            zf=last_state, status=code)
                    else
                        line = sosfilt_complex(sos, x(i,j,:), zf=last_state, status=code)
                    end if
                    if (code /= sos_status_success) return
                    y(i,j,:) = line
                    if (present(zf)) zf(:,i,j,:) = last_state
                end do
            end do
        end select
        if (present(status)) status = sos_status_success
    end function sosfilt_rank3_complex

    ! Rank-independent SOS state layout:
    !   (nsections, d1, ..., d_axis=2, ..., d_rank).
    ! For rank-15 inputs a rank-16 state would exceed Fortran 2018's limit,
    ! so rank-one packed states are accepted for all ranks.  Their ordering
    ! is the column-major flattening of the above conceptual state layout.
    pure logical function sos_nd_state_valid(state_dims, state_rank, dims, axis, &
        nsections) result(ok)
        integer, intent(in) :: state_dims(:) !! extents of caller-supplied state storage
        integer, intent(in) :: state_rank !! actual state-array rank
        integer, intent(in) :: dims(:) !! extents of the input data
        integer, intent(in) :: axis !! one-based filtering axis
        integer, intent(in) :: nsections !! number of biquad sections
        integer :: expected(size(dims)+1)

        expected = [nsections, dims]
        expected(axis+1) = 2
        ok = .false.
        if (state_rank == 1) then
            ok = int(state_dims(1),int64) == product(int(expected,int64))
        else if (state_rank == size(expected)) then
            ok = all(state_dims == expected)
        end if
    end function sos_nd_state_valid

    subroutine sosfilt_nd_real(sos, x, y, zi, zf, status, axis)
        real(dp), intent(in) :: sos(:,:) !! SOS matrix (nsections,6), with a0=1
        real(dp), contiguous, target, intent(in) :: x(..) !! data, rank 1 through 15
        real(dp), contiguous, target, intent(out) :: y(..) !! output with shape(x)
        real(dp), contiguous, target, intent(in), optional :: zi(..) !! per-channel or packed states
        real(dp), contiguous, target, intent(out), optional :: zf(..) !! final SOS delays
        integer, intent(out), optional :: status !! zero on success, one for invalid input
        integer, intent(in), optional :: axis !! one-based filtered axis, default last
        real(dp), pointer :: xin(:) => null(), yout(:) => null()
        real(dp), pointer :: zin(:) => null(), zout(:) => null()
        real(dp), allocatable :: state(:,:)
        real(dp) :: sample, output, nanval
        integer :: ax, nsec, dims(rank(x)), k, t
        integer(int64) :: stride, block, state_stride, state_block, outer, i, j, base, pos

        nanval = ieee_value(0.0_dp,ieee_quiet_nan)
        if (present(status)) status = sos_status_invalid_input
        if (size(y,kind=int64) > 0_int64) then
            call c_f_pointer(c_loc(y),yout,[size(y,kind=int64)])
            yout = nanval
        end if
        if (present(zf)) then
            if (size(zf,kind=int64) > 0_int64) then
                call c_f_pointer(c_loc(zf),zout,[size(zf,kind=int64)])
                zout = nanval
            end if
        end if
        if (rank(x) < 1 .or. rank(x) /= rank(y)) return
        if (any(shape(x) /= shape(y))) return
        ax = rank(x)
        if (present(axis)) ax = axis
        if (ax < 1 .or. ax > rank(x)) return
        if (.not. sos_real_valid(sos)) return
        nsec = size(sos,1)
        dims = shape(x)
        if (present(zi)) then
            if (.not. sos_nd_state_valid(shape(zi),rank(zi),dims,ax,nsec)) return
            if (size(zi,kind=int64) > 0_int64) then
                call c_f_pointer(c_loc(zi),zin,[size(zi,kind=int64)])
                if (.not. all(ieee_is_finite(zin))) return
            end if
        end if
        if (present(zf)) then
            if (.not. sos_nd_state_valid(shape(zf),rank(zf),dims,ax,nsec)) return
        end if
        if (size(x,kind=int64) > 0_int64) call c_f_pointer(c_loc(x),xin,[size(x,kind=int64)])
        allocate(state(nsec,2))
        stride = product(int(dims(:ax-1),int64))
        block = stride*int(dims(ax),int64)
        state_stride = int(nsec,int64)*stride
        state_block = 2_int64*state_stride
        outer = product(int(dims(ax+1:),int64))
        do j = 0_int64, outer-1_int64
            do i = 0_int64, stride-1_int64
                state = 0.0_dp
                base = j*state_block+i*int(nsec,int64)+1_int64
                if (present(zi)) then
                    do k = 1, nsec
                        state(k,1) = zin(base+k-1)
                        state(k,2) = zin(base+state_stride+k-1)
                    end do
                end if
                do t = 0, dims(ax)-1
                    pos = j*block+i+1_int64+int(t,int64)*stride
                    sample = xin(pos)
                    do k = 1, nsec
                        output = sos(k,1)*sample+state(k,1)
                        state(k,1) = sos(k,2)*sample-sos(k,5)*output+state(k,2)
                        state(k,2) = sos(k,3)*sample-sos(k,6)*output
                        sample = output
                    end do
                    yout(pos) = sample
                end do
                if (present(zf)) then
                    do k = 1, nsec
                        zout(base+k-1) = state(k,1)
                        zout(base+state_stride+k-1) = state(k,2)
                    end do
                end if
            end do
        end do
        if (present(status)) status = sos_status_success
    end subroutine sosfilt_nd_real

    subroutine sosfilt_nd_complex(sos, x, y, zi, zf, status, axis)
        complex(dp), intent(in) :: sos(:,:) !! SOS matrix (nsections,6), with a0=1
        complex(dp), contiguous, target, intent(in) :: x(..) !! data, rank 1 through 15
        complex(dp), contiguous, target, intent(out) :: y(..) !! output with shape(x)
        complex(dp), contiguous, target, intent(in), optional :: zi(..) !! initial SOS delays
        complex(dp), contiguous, target, intent(out), optional :: zf(..) !! final SOS delays
        integer, intent(out), optional :: status !! zero on success, one for invalid input
        integer, intent(in), optional :: axis !! one-based filtered axis, default last
        complex(dp), pointer :: xin(:) => null(), yout(:) => null()
        complex(dp), pointer :: zin(:) => null(), zout(:) => null()
        complex(dp), allocatable :: state(:,:)
        complex(dp) :: sample, output
        real(dp) :: nanval
        integer :: ax, nsec, dims(rank(x)), k, t
        integer(int64) :: stride, block, state_stride, state_block, outer, i, j, base, pos

        nanval = ieee_value(0.0_dp,ieee_quiet_nan)
        if (present(status)) status = sos_status_invalid_input
        if (size(y,kind=int64) > 0_int64) then
            call c_f_pointer(c_loc(y),yout,[size(y,kind=int64)])
            yout = cmplx(nanval,nanval,dp)
        end if
        if (present(zf)) then
            if (size(zf,kind=int64) > 0_int64) then
                call c_f_pointer(c_loc(zf),zout,[size(zf,kind=int64)])
                zout = cmplx(nanval,nanval,dp)
            end if
        end if
        if (rank(x) < 1 .or. rank(x) /= rank(y)) return
        if (any(shape(x) /= shape(y))) return
        ax = rank(x)
        if (present(axis)) ax = axis
        if (ax < 1 .or. ax > rank(x)) return
        if (.not. sos_complex_valid(sos)) return
        nsec = size(sos,1)
        dims = shape(x)
        if (present(zi)) then
            if (.not. sos_nd_state_valid(shape(zi),rank(zi),dims,ax,nsec)) return
            if (size(zi,kind=int64) > 0_int64) then
                call c_f_pointer(c_loc(zi),zin,[size(zi,kind=int64)])
                if (.not. (all(ieee_is_finite(real(zin,dp))) .and. &
                    all(ieee_is_finite(aimag(zin))))) return
            end if
        end if
        if (present(zf)) then
            if (.not. sos_nd_state_valid(shape(zf),rank(zf),dims,ax,nsec)) return
        end if
        if (size(x,kind=int64) > 0_int64) call c_f_pointer(c_loc(x),xin,[size(x,kind=int64)])
        allocate(state(nsec,2))
        stride = product(int(dims(:ax-1),int64))
        block = stride*int(dims(ax),int64)
        state_stride = int(nsec,int64)*stride
        state_block = 2_int64*state_stride
        outer = product(int(dims(ax+1:),int64))
        do j = 0_int64, outer-1_int64
            do i = 0_int64, stride-1_int64
                state = (0.0_dp,0.0_dp)
                base = j*state_block+i*int(nsec,int64)+1_int64
                if (present(zi)) then
                    do k = 1, nsec
                        state(k,1) = zin(base+k-1)
                        state(k,2) = zin(base+state_stride+k-1)
                    end do
                end if
                do t = 0, dims(ax)-1
                    pos = j*block+i+1_int64+int(t,int64)*stride
                    sample = xin(pos)
                    do k = 1, nsec
                        output = sos(k,1)*sample+state(k,1)
                        state(k,1) = sos(k,2)*sample-sos(k,5)*output+state(k,2)
                        state(k,2) = sos(k,3)*sample-sos(k,6)*output
                        sample = output
                    end do
                    yout(pos) = sample
                end do
                if (present(zf)) then
                    do k = 1, nsec
                        zout(base+k-1) = state(k,1)
                        zout(base+state_stride+k-1) = state(k,2)
                    end do
                end if
            end do
        end do
        if (present(status)) status = sos_status_success
    end subroutine sosfilt_nd_complex

    ! Convert only the small coefficient matrix, never the signal array.
    subroutine sosfilt_nd_real_sos_complex_input(sos, x, y, zi, zf, status, axis)
        real(dp), intent(in) :: sos(:,:) !! real SOS matrix, with a0=1
        complex(dp), contiguous, target, intent(in) :: x(..) !! complex input, rank 1 through 15
        complex(dp), contiguous, target, intent(out) :: y(..) !! complex filtered output
        complex(dp), contiguous, target, intent(in), optional :: zi(..) !! complex initial states
        complex(dp), contiguous, target, intent(out), optional :: zf(..) !! complex final states
        integer, intent(out), optional :: status !! zero on success, one for invalid input
        integer, intent(in), optional :: axis !! one-based filtered axis, default last

        call sosfilt_nd_complex(cmplx(sos,0.0_dp,dp),x,y,zi,zf,status,axis)
    end subroutine sosfilt_nd_real_sos_complex_input

    ! Keep real samples real in memory; promote one sample per filter step.
    subroutine sosfilt_nd_complex_sos_real_input(sos, x, y, zi, zf, status, axis)
        complex(dp), intent(in) :: sos(:,:) !! complex SOS matrix, with a0=1
        real(dp), contiguous, target, intent(in) :: x(..) !! real input, rank 1 through 15
        complex(dp), contiguous, target, intent(out) :: y(..) !! complex filtered output
        complex(dp), contiguous, target, intent(in), optional :: zi(..) !! complex initial states
        complex(dp), contiguous, target, intent(out), optional :: zf(..) !! complex final states
        integer, intent(out), optional :: status !! zero on success, one for invalid input
        integer, intent(in), optional :: axis !! one-based filtered axis, default last
        real(dp), pointer :: xin(:) => null()
        complex(dp), pointer :: yout(:) => null(), zin(:) => null(), zout(:) => null()
        complex(dp), allocatable :: state(:,:)
        complex(dp) :: sample, output
        real(dp) :: nanval
        integer :: ax, nsec, dims(rank(x)), k, t
        integer(int64) :: stride, block, state_stride, state_block, outer, i, j, base, pos

        nanval = ieee_value(0.0_dp,ieee_quiet_nan)
        if (present(status)) status = sos_status_invalid_input
        if (size(y,kind=int64) > 0_int64) then
            call c_f_pointer(c_loc(y),yout,[size(y,kind=int64)])
            yout = cmplx(nanval,nanval,dp)
        end if
        if (present(zf)) then
            if (size(zf,kind=int64) > 0_int64) then
                call c_f_pointer(c_loc(zf),zout,[size(zf,kind=int64)])
                zout = cmplx(nanval,nanval,dp)
            end if
        end if
        if (rank(x) < 1 .or. rank(x) /= rank(y)) return
        if (any(shape(x) /= shape(y))) return
        ax = rank(x)
        if (present(axis)) ax = axis
        if (ax < 1 .or. ax > rank(x)) return
        if (.not. sos_complex_valid(sos)) return
        nsec = size(sos,1)
        dims = shape(x)
        if (present(zi)) then
            if (.not. sos_nd_state_valid(shape(zi),rank(zi),dims,ax,nsec)) return
            if (size(zi,kind=int64) > 0_int64) then
                call c_f_pointer(c_loc(zi),zin,[size(zi,kind=int64)])
                if (.not. (all(ieee_is_finite(real(zin,dp))) .and. &
                    all(ieee_is_finite(aimag(zin))))) return
            end if
        end if
        if (present(zf)) then
            if (.not. sos_nd_state_valid(shape(zf),rank(zf),dims,ax,nsec)) return
        end if
        if (size(x,kind=int64) > 0_int64) call c_f_pointer(c_loc(x),xin,[size(x,kind=int64)])
        allocate(state(nsec,2))
        stride = product(int(dims(:ax-1),int64))
        block = stride*int(dims(ax),int64)
        state_stride = int(nsec,int64)*stride
        state_block = 2_int64*state_stride
        outer = product(int(dims(ax+1:),int64))
        do j = 0_int64, outer-1_int64
            do i = 0_int64, stride-1_int64
                state = (0.0_dp,0.0_dp)
                base = j*state_block+i*int(nsec,int64)+1_int64
                if (present(zi)) then
                    do k = 1, nsec
                        state(k,1) = zin(base+k-1)
                        state(k,2) = zin(base+state_stride+k-1)
                    end do
                end if
                do t = 0, dims(ax)-1
                    pos = j*block+i+1_int64+int(t,int64)*stride
                    sample = cmplx(xin(pos),0.0_dp,dp)
                    do k = 1, nsec
                        output = sos(k,1)*sample+state(k,1)
                        state(k,1) = sos(k,2)*sample-sos(k,5)*output+state(k,2)
                        state(k,2) = sos(k,3)*sample-sos(k,6)*output
                        sample = output
                    end do
                    yout(pos) = sample
                end do
                if (present(zf)) then
                    do k = 1, nsec
                        zout(base+k-1) = state(k,1)
                        zout(base+state_stride+k-1) = state(k,2)
                    end do
                end if
            end do
        end do
        if (present(status)) status = sos_status_success
    end subroutine sosfilt_nd_complex_sos_real_input
end module scifort_signal_sos
