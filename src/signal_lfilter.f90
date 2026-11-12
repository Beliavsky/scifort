! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Real and complex 1D transposed direct-form-II digital filtering.
module scifort_signal_lfilter
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
    use, intrinsic :: iso_c_binding, only : c_loc, c_f_pointer
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_kinds, only : dp
    implicit none
    private

    integer, parameter, public :: lfilter_status_success = 0
    integer, parameter, public :: lfilter_status_invalid_input = 1
    integer, parameter, public :: lfilter_status_no_steady_state = 2

    interface lfilter
        module procedure lfilter_real
        module procedure lfilter_complex
        module procedure lfilter_matrix_real, lfilter_matrix_complex
        module procedure lfilter_rank3_real, lfilter_rank3_complex
        module procedure lfilter_mixed_rrc_1d
        module procedure lfilter_mixed_rcr_1d
        module procedure lfilter_mixed_rcc_1d
        module procedure lfilter_mixed_crr_1d
        module procedure lfilter_mixed_crc_1d
        module procedure lfilter_mixed_ccr_1d
        module procedure lfilter_mixed_rrc_2d
        module procedure lfilter_mixed_rcr_2d
        module procedure lfilter_mixed_rcc_2d
        module procedure lfilter_mixed_crr_2d
        module procedure lfilter_mixed_crc_2d
        module procedure lfilter_mixed_ccr_2d
        module procedure lfilter_mixed_rrc_3d
        module procedure lfilter_mixed_rcr_3d
        module procedure lfilter_mixed_rcc_3d
        module procedure lfilter_mixed_crr_3d
        module procedure lfilter_mixed_crc_3d
        module procedure lfilter_mixed_ccr_3d
    end interface lfilter

    interface lfilter_nd
        module procedure lfilter_nd_real, lfilter_nd_complex
        module procedure lfilter_nd_mixed_rrc, lfilter_nd_mixed_rcr
        module procedure lfilter_nd_mixed_rcc, lfilter_nd_mixed_crr
        module procedure lfilter_nd_mixed_crc, lfilter_nd_mixed_ccr
    end interface lfilter_nd

    interface lfilter_zi
        module procedure lfilter_zi_real
        module procedure lfilter_zi_complex
    end interface lfilter_zi

    interface lfiltic
        module procedure lfiltic_real
        module procedure lfiltic_complex
    end interface lfiltic

    public :: lfilter, lfilter_nd, lfilter_zi, lfiltic

contains

    function lfilter_real(b, a, x, zi, zf, status) result(y)
        real(dp), intent(in) :: b(:) !! numerator transfer-function coefficients, lowest lag first
        real(dp), intent(in) :: a(:) !! denominator coefficients, nonzero leading coefficient
        real(dp), intent(in) :: x(:) !! input samples in chronological order
        real(dp), intent(in), optional :: zi(:) !! starting direct-form-II delay values
        real(dp), allocatable, intent(out), optional :: zf(:) !! ending delay values
        integer, intent(out), optional :: status !! zero on success; nonzero for invalid coefficients
        real(dp), allocatable :: y(:)

        real(dp), allocatable :: aa(:), bb(:), state(:)
        real(dp) :: nanval
        integer :: nstate, i, k

        nanval = ieee_value(0.0_dp, ieee_quiet_nan)
        allocate(y(size(x)))
        y = nanval
        if (present(zf)) allocate(zf(0))
        if (present(status)) status = lfilter_status_invalid_input
        if (size(a) < 1 .or. size(b) < 1) return
        if (abs(a(1)) < tiny(nanval)) return
        if (.not. (all(ieee_is_finite(b)) .and. all(ieee_is_finite(a)))) return
        nstate = max(size(a), size(b)) - 1
        if (present(zi)) then
            if (size(zi) /= nstate) return
            if (.not. (all(ieee_is_finite(zi)))) return
        end if
        allocate(aa(nstate + 1), bb(nstate + 1), state(nstate))
        aa = 0.0_dp
        bb = 0.0_dp
        aa(1:size(a)) = a / a(1)
        bb(1:size(b)) = b / a(1)
        state = 0.0_dp
        if (present(zi)) state = zi

        do k = 1, size(x)
            y(k) = bb(1) * x(k)
            if (nstate > 0) y(k) = y(k) + state(1)
            do i = 1, nstate - 1
                state(i) = state(i + 1) + bb(i + 1) * x(k) - aa(i + 1) * y(k)
            end do
            if (nstate > 0) then
                state(nstate) = bb(nstate + 1) * x(k) - aa(nstate + 1) * y(k)
            end if
        end do
        if (present(zf)) zf = state
        if (present(status)) status = lfilter_status_success
    end function lfilter_real

    function lfilter_zi_real(b, a, status) result(zi)
        real(dp), intent(in) :: b(:) !! numerator transfer-function coefficients
        real(dp), intent(in) :: a(:) !! denominator coefficients, nonzero leading coefficient
        integer, intent(out), optional :: status !! zero on success; nonzero if invalid or pole at one
        real(dp), allocatable :: zi(:)

        real(dp), allocatable :: aa(:), bb(:)
        real(dp) :: steady, denominator
        real(dp) :: nanval
        integer :: nstate, i

        nanval = ieee_value(0.0_dp, ieee_quiet_nan)
        if (present(status)) status = lfilter_status_invalid_input
        nstate = max(0, max(size(a), size(b)) - 1)
        allocate(zi(nstate))
        zi = nanval
        if (size(a) < 1 .or. size(b) < 1) return
        if (abs(a(1)) < tiny(nanval)) return
        if (.not. (all(ieee_is_finite(b)) .and. all(ieee_is_finite(a)))) return
        allocate(aa(nstate + 1), bb(nstate + 1))
        aa = 0.0_dp
        bb = 0.0_dp
        aa(1:size(a)) = a / a(1)
        bb(1:size(b)) = b / a(1)
        denominator = sum(aa)
        if (abs(denominator) <= 32.0_dp * epsilon(nanval) * sum(abs(aa))) then
            if (present(status)) status = lfilter_status_no_steady_state
            return
        end if
        steady = sum(bb) / denominator
        if (.not. ieee_is_finite(steady)) return
        if (nstate > 0) then
            zi(nstate) = bb(nstate + 1) - aa(nstate + 1) * steady
            do i = nstate - 1, 1, -1
                zi(i) = zi(i + 1) + bb(i + 1) - aa(i + 1) * steady
            end do
        end if
        if (present(status)) status = lfilter_status_success
    end function lfilter_zi_real

    function lfilter_complex(b, a, x, zi, zf, status) result(y)
        complex(dp), intent(in) :: b(:) !! numerator transfer-function coefficients, lowest lag first
        complex(dp), intent(in) :: a(:) !! denominator coefficients, nonzero leading coefficient
        complex(dp), intent(in) :: x(:) !! input samples in chronological order
        complex(dp), intent(in), optional :: zi(:) !! starting direct-form-II delay values
        complex(dp), allocatable, intent(out), optional :: zf(:) !! ending delay values
        integer, intent(out), optional :: status !! zero on success; nonzero for invalid coefficients
        complex(dp), allocatable :: y(:)

        complex(dp), allocatable :: aa(:), bb(:), state(:)
        real(dp) :: nanval
        integer :: nstate, i, k

        nanval = ieee_value(0.0_dp, ieee_quiet_nan)
        allocate(y(size(x)))
        y = cmplx(nanval, nanval, dp)
        if (present(zf)) allocate(zf(0))
        if (present(status)) status = lfilter_status_invalid_input
        if (size(a) < 1 .or. size(b) < 1) return
        if (abs(a(1)) < tiny(nanval)) return
        if (.not. (all(ieee_is_finite(real(b, dp))) .and. all(ieee_is_finite(aimag(b))) .and. &
            all(ieee_is_finite(real(a, dp))) .and. all(ieee_is_finite(aimag(a))))) return
        nstate = max(size(a), size(b)) - 1
        if (present(zi)) then
            if (size(zi) /= nstate) return
            if (.not. (all(ieee_is_finite(real(zi, dp))) .and. all(ieee_is_finite(aimag(zi))))) return
        end if
        allocate(aa(nstate + 1), bb(nstate + 1), state(nstate))
        aa = (0.0_dp, 0.0_dp)
        bb = (0.0_dp, 0.0_dp)
        aa(1:size(a)) = a / a(1)
        bb(1:size(b)) = b / a(1)
        state = (0.0_dp, 0.0_dp)
        if (present(zi)) state = zi

        do k = 1, size(x)
            y(k) = bb(1) * x(k)
            if (nstate > 0) y(k) = y(k) + state(1)
            do i = 1, nstate - 1
                state(i) = state(i + 1) + bb(i + 1) * x(k) - aa(i + 1) * y(k)
            end do
            if (nstate > 0) then
                state(nstate) = bb(nstate + 1) * x(k) - aa(nstate + 1) * y(k)
            end if
        end do
        if (present(zf)) zf = state
        if (present(status)) status = lfilter_status_success
    end function lfilter_complex

    function lfilter_zi_complex(b, a, status) result(zi)
        complex(dp), intent(in) :: b(:) !! numerator transfer-function coefficients
        complex(dp), intent(in) :: a(:) !! denominator coefficients, nonzero leading coefficient
        integer, intent(out), optional :: status !! zero on success; nonzero if invalid or pole at one
        complex(dp), allocatable :: zi(:)

        complex(dp), allocatable :: aa(:), bb(:)
        complex(dp) :: steady, denominator
        real(dp) :: nanval
        integer :: nstate, i

        nanval = ieee_value(0.0_dp, ieee_quiet_nan)
        if (present(status)) status = lfilter_status_invalid_input
        nstate = max(0, max(size(a), size(b)) - 1)
        allocate(zi(nstate))
        zi = cmplx(nanval, nanval, dp)
        if (size(a) < 1 .or. size(b) < 1) return
        if (abs(a(1)) < tiny(nanval)) return
        if (.not. (all(ieee_is_finite(real(b, dp))) .and. all(ieee_is_finite(aimag(b))) .and. &
            all(ieee_is_finite(real(a, dp))) .and. all(ieee_is_finite(aimag(a))))) return
        allocate(aa(nstate + 1), bb(nstate + 1))
        aa = (0.0_dp, 0.0_dp)
        bb = (0.0_dp, 0.0_dp)
        aa(1:size(a)) = a / a(1)
        bb(1:size(b)) = b / a(1)
        denominator = sum(aa)
        if (abs(denominator) <= 32.0_dp * epsilon(nanval) * sum(abs(aa))) then
            if (present(status)) status = lfilter_status_no_steady_state
            return
        end if
        steady = sum(bb) / denominator
        if (.not. (ieee_is_finite(real(steady, dp)) .and. ieee_is_finite(aimag(steady)))) return
        if (nstate > 0) then
            zi(nstate) = bb(nstate + 1) - aa(nstate + 1) * steady
            do i = nstate - 1, 1, -1
                zi(i) = zi(i + 1) + bb(i + 1) - aa(i + 1) * steady
            end do
        end if
        if (present(status)) status = lfilter_status_success
    end function lfilter_zi_complex


    function lfiltic_real(b, a, y, x, status) result(zi)
        ! Initialize the transposed direct-form-II state from past observations.
        ! Histories are newest-first: y(1)=y[-1], x(1)=x[-1].
        real(dp), intent(in) :: b(:) !! numerator coefficients in increasing lag order
        real(dp), intent(in) :: a(:) !! denominator coefficients with nonzero a(1)
        real(dp), intent(in) :: y(:) !! previous outputs, most recent first
        real(dp), intent(in), optional :: x(:) !! previous inputs, most recent first; default zeros
        integer, intent(out), optional :: status !! zero on success, nonzero on invalid inputs
        real(dp), allocatable :: zi(:)

        real(dp), allocatable :: aa(:), bb(:)
        real(dp) :: nanval, history_x, history_y
        integer :: nstate, i, j, lag

        nanval = ieee_value(0.0_dp, ieee_quiet_nan)
        if (present(status)) status = lfilter_status_invalid_input
        nstate = max(0, max(size(a), size(b)) - 1)
        allocate(zi(nstate))
        zi = nanval
        if (size(a) == 0 .or. size(b) == 0) return
        if (abs(a(1)) < tiny(1.0_dp)) return
        if (.not. all(ieee_is_finite(a)) .or. .not. all(ieee_is_finite(b))) return
        if (.not. all(ieee_is_finite(y))) return
        if (present(x)) then
            if (.not. all(ieee_is_finite(x))) return
        end if

        allocate(aa(nstate+1), bb(nstate+1))
        aa = 0.0_dp
        bb = 0.0_dp
        aa(1:size(a)) = a / a(1)
        bb(1:size(b)) = b / a(1)
        if (.not. all(ieee_is_finite(aa)) .or. .not. all(ieee_is_finite(bb))) return
        zi = 0.0_dp
        do i = 1, nstate
            do j = i+1, nstate+1
                lag = j-i
                history_y = 0.0_dp
                history_x = 0.0_dp
                if (lag <= size(y)) history_y = y(lag)
                if (present(x)) then
                    if (lag <= size(x)) history_x = x(lag)
                end if
                zi(i) = zi(i) + bb(j)*history_x - aa(j)*history_y
            end do
        end do
        if (.not. all(ieee_is_finite(zi))) then
            zi = nanval
            return
        end if
        if (present(status)) status = lfilter_status_success
    end function lfiltic_real

    function lfiltic_complex(b, a, y, x, status) result(zi)
        ! Complex-coefficient counterpart of lfiltic_real.
        complex(dp), intent(in) :: b(:) !! numerator coefficients in increasing lag order
        complex(dp), intent(in) :: a(:) !! denominator coefficients with nonzero a(1)
        complex(dp), intent(in) :: y(:) !! previous outputs, most recent first
        complex(dp), intent(in), optional :: x(:) !! previous inputs, most recent first; default zeros
        integer, intent(out), optional :: status !! zero on success, nonzero on invalid inputs
        complex(dp), allocatable :: zi(:)

        complex(dp), allocatable :: aa(:), bb(:)
        complex(dp) :: history_x, history_y
        real(dp) :: nanval
        integer :: nstate, i, j, lag

        nanval = ieee_value(0.0_dp, ieee_quiet_nan)
        if (present(status)) status = lfilter_status_invalid_input
        nstate = max(0, max(size(a), size(b)) - 1)
        allocate(zi(nstate))
        zi = cmplx(nanval, nanval, dp)
        if (size(a) == 0 .or. size(b) == 0) return
        if (abs(a(1)) < tiny(1.0_dp)) return
        if (.not. all(ieee_is_finite(real(a, dp))) .or. &
            .not. all(ieee_is_finite(aimag(a)))) return
        if (.not. all(ieee_is_finite(real(b, dp))) .or. &
            .not. all(ieee_is_finite(aimag(b)))) return
        if (.not. all(ieee_is_finite(real(y, dp))) .or. &
            .not. all(ieee_is_finite(aimag(y)))) return
        if (present(x)) then
            if (.not. all(ieee_is_finite(real(x, dp))) .or. &
                .not. all(ieee_is_finite(aimag(x)))) return
        end if

        allocate(aa(nstate+1), bb(nstate+1))
        aa = (0.0_dp, 0.0_dp)
        bb = (0.0_dp, 0.0_dp)
        aa(1:size(a)) = a / a(1)
        bb(1:size(b)) = b / a(1)
        if (.not. all(ieee_is_finite(real(aa, dp))) .or. &
            .not. all(ieee_is_finite(aimag(aa)))) return
        if (.not. all(ieee_is_finite(real(bb, dp))) .or. &
            .not. all(ieee_is_finite(aimag(bb)))) return
        zi = (0.0_dp, 0.0_dp)
        do i = 1, nstate
            do j = i+1, nstate+1
                lag = j-i
                history_y = (0.0_dp, 0.0_dp)
                history_x = (0.0_dp, 0.0_dp)
                if (lag <= size(y)) history_y = y(lag)
                if (present(x)) then
                    if (lag <= size(x)) history_x = x(lag)
                end if
                zi(i) = zi(i) + bb(j)*history_x - aa(j)*history_y
            end do
        end do
        if (.not. all(ieee_is_finite(real(zi, dp))) .or. &
            .not. all(ieee_is_finite(aimag(zi)))) then
            zi = cmplx(nanval, nanval, dp)
            return
        end if
        if (present(status)) status = lfilter_status_success
    end function lfiltic_complex


    ! The rank-two variants apply the one-dimensional transposed DF-II filter
    ! independently to each row or column. The axis is one-based, default 2.

    function lfilter_matrix_real(b,a,x,zi,zf,status,axis) result(y)
        real(dp), intent(in) :: b(:) !! numerator coefficients in increasing lag order
        real(dp), intent(in) :: a(:) !! denominator coefficients, with nonzero leading value
        real(dp), intent(in) :: x(:,:) !! input signal matrix
        real(dp), intent(in), optional :: zi(:,:) !! initial states, shaped as x with axis length nstate
        real(dp), allocatable, intent(out), optional :: zf(:,:) !! final states, shaped as zi
        integer, intent(out), optional :: status !! zero on success
        integer, intent(in), optional :: axis !! filter axis: 1 or 2; default 2
        real(dp), allocatable :: y(:,:)
        real(dp), allocatable :: row(:), last_state(:)
        real(dp) :: nanval
        integer :: along, nstate, i, code

        nanval = ieee_value(0.0_dp, ieee_quiet_nan)
        allocate(y(size(x,1),size(x,2)))
        y = nanval
        if (present(zf)) allocate(zf(0,0))
        if (present(status)) status = lfilter_status_invalid_input
        along = 2
        if (present(axis)) along = axis
        if (along /= 1 .and. along /= 2) return
        if (size(a) < 1 .or. size(b) < 1) return
        if (abs(a(1)) < tiny(1.0_dp)) return
        if (.not. (all(ieee_is_finite(a)) .and. all(ieee_is_finite(b)))) return
        nstate = max(size(b),size(a)) - 1
        if (present(zi)) then
            if (along == 1) then
                if (any(shape(zi) /= [nstate,size(x,2)])) return
            else
                if (any(shape(zi) /= [size(x,1),nstate])) return
            end if
            if (.not. (all(ieee_is_finite(zi)))) return
        end if
        if (present(zf)) then
            deallocate(zf)
            if (along == 1) then
                allocate(zf(nstate,size(x,2)))
            else
                allocate(zf(size(x,1),nstate))
            end if
        end if
        if (along == 1) then
            do i = 1, size(x,2)
                if (present(zi)) then
                    row = lfilter_real(b,a,x(:,i),zi=zi(:,i),zf=last_state,status=code)
                else
                    row = lfilter_real(b,a,x(:,i),zf=last_state,status=code)
                end if
                if (code /= lfilter_status_success) return
                y(:,i) = row
                if (present(zf)) zf(:,i) = last_state
            end do
        else
            do i = 1, size(x,1)
                if (present(zi)) then
                    row = lfilter_real(b,a,x(i,:),zi=zi(i,:),zf=last_state,status=code)
                else
                    row = lfilter_real(b,a,x(i,:),zf=last_state,status=code)
                end if
                if (code /= lfilter_status_success) return
                y(i,:) = row
                if (present(zf)) zf(i,:) = last_state
            end do
        end if
        if (present(status)) status = lfilter_status_success
    end function lfilter_matrix_real

    function lfilter_matrix_complex(b,a,x,zi,zf,status,axis) result(y)
        complex(dp), intent(in) :: b(:) !! numerator coefficients in increasing lag order
        complex(dp), intent(in) :: a(:) !! denominator coefficients, with nonzero leading value
        complex(dp), intent(in) :: x(:,:) !! input signal matrix
        complex(dp), intent(in), optional :: zi(:,:) !! initial states, shaped as x with axis length nstate
        complex(dp), allocatable, intent(out), optional :: zf(:,:) !! final states, shaped as zi
        integer, intent(out), optional :: status !! zero on success
        integer, intent(in), optional :: axis !! filter axis: 1 or 2; default 2
        complex(dp), allocatable :: y(:,:)
        complex(dp), allocatable :: row(:), last_state(:)
        real(dp) :: nanval
        integer :: along, nstate, i, code

        nanval = ieee_value(0.0_dp, ieee_quiet_nan)
        allocate(y(size(x,1),size(x,2)))
        y = cmplx(nanval,nanval,dp)
        if (present(zf)) allocate(zf(0,0))
        if (present(status)) status = lfilter_status_invalid_input
        along = 2
        if (present(axis)) along = axis
        if (along /= 1 .and. along /= 2) return
        if (size(a) < 1 .or. size(b) < 1) return
        if (abs(a(1)) < tiny(1.0_dp)) return
        if (.not. (all(ieee_is_finite(real(a,dp))) .and. all(ieee_is_finite(aimag(a))) .and. &
            all(ieee_is_finite(real(b,dp))) .and. all(ieee_is_finite(aimag(b))))) return
        nstate = max(size(b),size(a)) - 1
        if (present(zi)) then
            if (along == 1) then
                if (any(shape(zi) /= [nstate,size(x,2)])) return
            else
                if (any(shape(zi) /= [size(x,1),nstate])) return
            end if
            if (.not. (all(ieee_is_finite(real(zi,dp))) .and. all(ieee_is_finite(aimag(zi))))) return
        end if
        if (present(zf)) then
            deallocate(zf)
            if (along == 1) then
                allocate(zf(nstate,size(x,2)))
            else
                allocate(zf(size(x,1),nstate))
            end if
        end if
        if (along == 1) then
            do i = 1, size(x,2)
                if (present(zi)) then
                    row = lfilter_complex(b,a,x(:,i),zi=zi(:,i),zf=last_state,status=code)
                else
                    row = lfilter_complex(b,a,x(:,i),zf=last_state,status=code)
                end if
                if (code /= lfilter_status_success) return
                y(:,i) = row
                if (present(zf)) zf(:,i) = last_state
            end do
        else
            do i = 1, size(x,1)
                if (present(zi)) then
                    row = lfilter_complex(b,a,x(i,:),zi=zi(i,:),zf=last_state,status=code)
                else
                    row = lfilter_complex(b,a,x(i,:),zf=last_state,status=code)
                end if
                if (code /= lfilter_status_success) return
                y(i,:) = row
                if (present(zf)) zf(i,:) = last_state
            end do
        end if
        if (present(status)) status = lfilter_status_success
    end function lfilter_matrix_complex


    ! Rank-three filtering follows the same one-based axis convention as rank two.

    function lfilter_rank3_real(b,a,x,zi,zf,status,axis) result(y)
        real(dp), intent(in) :: b(:) !! numerator coefficients in increasing lag order
        real(dp), intent(in) :: a(:) !! denominator coefficients, nonzero leading value
        real(dp), intent(in) :: x(:,:,:) !! input rank-three array
        real(dp), intent(in), optional :: zi(:,:,:) !! initial DF-II states on each filtered lane
        real(dp), allocatable, intent(out), optional :: zf(:,:,:) !! final states
        integer, intent(out), optional :: status !! zero on success
        integer, intent(in), optional :: axis !! one-based filtering axis, default 3
        real(dp), allocatable :: y(:,:,:)
        real(dp), allocatable :: line(:), last_state(:)
        integer :: along, state_shape(3), i, j, code
        real(dp) :: nanval

        nanval = ieee_value(0.0_dp, ieee_quiet_nan)
        allocate(y(size(x,1),size(x,2),size(x,3)))
        y = nanval
        if (present(zf)) allocate(zf(0,0,0))
        if (present(status)) status = lfilter_status_invalid_input
        along = 3
        if (present(axis)) along = axis
        if (along < 1 .or. along > 3) return
        if (size(a) < 1 .or. size(b) < 1) return
        if (abs(a(1)) < tiny(1.0_dp)) return
        if (.not. (all(ieee_is_finite(a)) .and. all(ieee_is_finite(b)))) return
        state_shape = shape(x)
        state_shape(along) = max(size(a),size(b)) - 1
        if (present(zi)) then
            if (any(shape(zi) /= state_shape)) return
            if (.not. (all(ieee_is_finite(zi)))) return
        end if
        if (present(zf)) then
            deallocate(zf)
            allocate(zf(state_shape(1),state_shape(2),state_shape(3)))
        end if
        select case (along)
        case (1)
            do j = 1, size(x,3)
                do i = 1, size(x,2)
                    if (present(zi)) then
                        line = lfilter_real(b,a,x(:,i,j),zi=zi(:,i,j),zf=last_state,status=code)
                    else
                        line = lfilter_real(b,a,x(:,i,j),zf=last_state,status=code)
                    end if
                    if (code /= lfilter_status_success) return
                    y(:,i,j) = line
                    if (present(zf)) zf(:,i,j) = last_state
                end do
            end do
        case (2)
            do j = 1, size(x,3)
                do i = 1, size(x,1)
                    if (present(zi)) then
                        line = lfilter_real(b,a,x(i,:,j),zi=zi(i,:,j),zf=last_state,status=code)
                    else
                        line = lfilter_real(b,a,x(i,:,j),zf=last_state,status=code)
                    end if
                    if (code /= lfilter_status_success) return
                    y(i,:,j) = line
                    if (present(zf)) zf(i,:,j) = last_state
                end do
            end do
        case (3)
            do j = 1, size(x,2)
                do i = 1, size(x,1)
                    if (present(zi)) then
                        line = lfilter_real(b,a,x(i,j,:),zi=zi(i,j,:),zf=last_state,status=code)
                    else
                        line = lfilter_real(b,a,x(i,j,:),zf=last_state,status=code)
                    end if
                    if (code /= lfilter_status_success) return
                    y(i,j,:) = line
                    if (present(zf)) zf(i,j,:) = last_state
                end do
            end do
        end select
        if (present(status)) status = lfilter_status_success
    end function lfilter_rank3_real

    function lfilter_rank3_complex(b,a,x,zi,zf,status,axis) result(y)
        complex(dp), intent(in) :: b(:) !! numerator coefficients in increasing lag order
        complex(dp), intent(in) :: a(:) !! denominator coefficients, nonzero leading value
        complex(dp), intent(in) :: x(:,:,:) !! input rank-three array
        complex(dp), intent(in), optional :: zi(:,:,:) !! initial DF-II states on each filtered lane
        complex(dp), allocatable, intent(out), optional :: zf(:,:,:) !! final states
        integer, intent(out), optional :: status !! zero on success
        integer, intent(in), optional :: axis !! one-based filtering axis, default 3
        complex(dp), allocatable :: y(:,:,:)
        complex(dp), allocatable :: line(:), last_state(:)
        integer :: along, state_shape(3), i, j, code
        real(dp) :: nanval

        nanval = ieee_value(0.0_dp, ieee_quiet_nan)
        allocate(y(size(x,1),size(x,2),size(x,3)))
        y = cmplx(nanval,nanval,dp)
        if (present(zf)) allocate(zf(0,0,0))
        if (present(status)) status = lfilter_status_invalid_input
        along = 3
        if (present(axis)) along = axis
        if (along < 1 .or. along > 3) return
        if (size(a) < 1 .or. size(b) < 1) return
        if (abs(a(1)) < tiny(1.0_dp)) return
        if (.not. (all(ieee_is_finite(real(a,dp))) .and. all(ieee_is_finite(aimag(a))) .and. &
            all(ieee_is_finite(real(b,dp))) .and. all(ieee_is_finite(aimag(b))))) return
        state_shape = shape(x)
        state_shape(along) = max(size(a),size(b)) - 1
        if (present(zi)) then
            if (any(shape(zi) /= state_shape)) return
            if (.not. (all(ieee_is_finite(real(zi,dp))) .and. all(ieee_is_finite(aimag(zi))))) return
        end if
        if (present(zf)) then
            deallocate(zf)
            allocate(zf(state_shape(1),state_shape(2),state_shape(3)))
        end if
        select case (along)
        case (1)
            do j = 1, size(x,3)
                do i = 1, size(x,2)
                    if (present(zi)) then
                        line = lfilter_complex(b,a,x(:,i,j),zi=zi(:,i,j),zf=last_state,status=code)
                    else
                        line = lfilter_complex(b,a,x(:,i,j),zf=last_state,status=code)
                    end if
                    if (code /= lfilter_status_success) return
                    y(:,i,j) = line
                    if (present(zf)) zf(:,i,j) = last_state
                end do
            end do
        case (2)
            do j = 1, size(x,3)
                do i = 1, size(x,1)
                    if (present(zi)) then
                        line = lfilter_complex(b,a,x(i,:,j),zi=zi(i,:,j),zf=last_state,status=code)
                    else
                        line = lfilter_complex(b,a,x(i,:,j),zf=last_state,status=code)
                    end if
                    if (code /= lfilter_status_success) return
                    y(i,:,j) = line
                    if (present(zf)) zf(i,:,j) = last_state
                end do
            end do
        case (3)
            do j = 1, size(x,2)
                do i = 1, size(x,1)
                    if (present(zi)) then
                        line = lfilter_complex(b,a,x(i,j,:),zi=zi(i,j,:),zf=last_state,status=code)
                    else
                        line = lfilter_complex(b,a,x(i,j,:),zf=last_state,status=code)
                    end if
                    if (code /= lfilter_status_success) return
                    y(i,j,:) = line
                    if (present(zf)) zf(i,j,:) = last_state
                end do
            end do
        end select
        if (present(status)) status = lfilter_status_success
    end function lfilter_rank3_complex

    function lfilter_mixed_rrc_1d(b, a, x, zi, zf, status) result(y)
        real(dp), intent(in) :: b(:) !! numerator coefficients
        real(dp), intent(in) :: a(:) !! denominator coefficients
        complex(dp), intent(in) :: x(:) !! input samples
        complex(dp), intent(in), optional :: zi(:) !! initial complex delay states
        complex(dp), allocatable, intent(out), optional :: zf(:) !! final complex delay states
        integer, intent(out), optional :: status !! success or invalid-input status
        complex(dp), allocatable :: y(:)

        ! Promote mixed inputs to a common double-complex working type.
        y = lfilter_complex(cmplx(b,0.0_dp,dp),cmplx(a,0.0_dp,dp), &
            x,zi=zi,zf=zf,status=status)
    end function lfilter_mixed_rrc_1d

    function lfilter_mixed_rcr_1d(b, a, x, zi, zf, status) result(y)
        real(dp), intent(in) :: b(:) !! numerator coefficients
        complex(dp), intent(in) :: a(:) !! denominator coefficients
        real(dp), intent(in) :: x(:) !! input samples
        complex(dp), intent(in), optional :: zi(:) !! initial complex delay states
        complex(dp), allocatable, intent(out), optional :: zf(:) !! final complex delay states
        integer, intent(out), optional :: status !! success or invalid-input status
        complex(dp), allocatable :: y(:)

        ! Promote mixed inputs to a common double-complex working type.
        y = lfilter_complex(cmplx(b,0.0_dp,dp),a, &
            cmplx(x,0.0_dp,dp),zi=zi,zf=zf,status=status)
    end function lfilter_mixed_rcr_1d

    function lfilter_mixed_rcc_1d(b, a, x, zi, zf, status) result(y)
        real(dp), intent(in) :: b(:) !! numerator coefficients
        complex(dp), intent(in) :: a(:) !! denominator coefficients
        complex(dp), intent(in) :: x(:) !! input samples
        complex(dp), intent(in), optional :: zi(:) !! initial complex delay states
        complex(dp), allocatable, intent(out), optional :: zf(:) !! final complex delay states
        integer, intent(out), optional :: status !! success or invalid-input status
        complex(dp), allocatable :: y(:)

        ! Promote mixed inputs to a common double-complex working type.
        y = lfilter_complex(cmplx(b,0.0_dp,dp),a, &
            x,zi=zi,zf=zf,status=status)
    end function lfilter_mixed_rcc_1d

    function lfilter_mixed_crr_1d(b, a, x, zi, zf, status) result(y)
        complex(dp), intent(in) :: b(:) !! numerator coefficients
        real(dp), intent(in) :: a(:) !! denominator coefficients
        real(dp), intent(in) :: x(:) !! input samples
        complex(dp), intent(in), optional :: zi(:) !! initial complex delay states
        complex(dp), allocatable, intent(out), optional :: zf(:) !! final complex delay states
        integer, intent(out), optional :: status !! success or invalid-input status
        complex(dp), allocatable :: y(:)

        ! Promote mixed inputs to a common double-complex working type.
        y = lfilter_complex(b,cmplx(a,0.0_dp,dp), &
            cmplx(x,0.0_dp,dp),zi=zi,zf=zf,status=status)
    end function lfilter_mixed_crr_1d

    function lfilter_mixed_crc_1d(b, a, x, zi, zf, status) result(y)
        complex(dp), intent(in) :: b(:) !! numerator coefficients
        real(dp), intent(in) :: a(:) !! denominator coefficients
        complex(dp), intent(in) :: x(:) !! input samples
        complex(dp), intent(in), optional :: zi(:) !! initial complex delay states
        complex(dp), allocatable, intent(out), optional :: zf(:) !! final complex delay states
        integer, intent(out), optional :: status !! success or invalid-input status
        complex(dp), allocatable :: y(:)

        ! Promote mixed inputs to a common double-complex working type.
        y = lfilter_complex(b,cmplx(a,0.0_dp,dp), &
            x,zi=zi,zf=zf,status=status)
    end function lfilter_mixed_crc_1d

    function lfilter_mixed_ccr_1d(b, a, x, zi, zf, status) result(y)
        complex(dp), intent(in) :: b(:) !! numerator coefficients
        complex(dp), intent(in) :: a(:) !! denominator coefficients
        real(dp), intent(in) :: x(:) !! input samples
        complex(dp), intent(in), optional :: zi(:) !! initial complex delay states
        complex(dp), allocatable, intent(out), optional :: zf(:) !! final complex delay states
        integer, intent(out), optional :: status !! success or invalid-input status
        complex(dp), allocatable :: y(:)

        ! Promote mixed inputs to a common double-complex working type.
        y = lfilter_complex(b,a, &
            cmplx(x,0.0_dp,dp),zi=zi,zf=zf,status=status)
    end function lfilter_mixed_ccr_1d

    function lfilter_mixed_rrc_2d(b, a, x, zi, zf, status, axis) result(y)
        real(dp), intent(in) :: b(:) !! numerator coefficients
        real(dp), intent(in) :: a(:) !! denominator coefficients
        complex(dp), intent(in) :: x(:,:) !! input samples
        complex(dp), intent(in), optional :: zi(:,:) !! initial complex delay states
        complex(dp), allocatable, intent(out), optional :: zf(:,:) !! final complex delay states
        integer, intent(out), optional :: status !! success or invalid-input status
        integer, intent(in), optional :: axis !! filter axis, default 2
        complex(dp), allocatable :: y(:,:)

        ! Promote mixed inputs to a common double-complex working type.
        y = lfilter_matrix_complex(cmplx(b,0.0_dp,dp),cmplx(a,0.0_dp,dp), &
            x,zi=zi,zf=zf,status=status,axis=axis)
    end function lfilter_mixed_rrc_2d

    function lfilter_mixed_rcr_2d(b, a, x, zi, zf, status, axis) result(y)
        real(dp), intent(in) :: b(:) !! numerator coefficients
        complex(dp), intent(in) :: a(:) !! denominator coefficients
        real(dp), intent(in) :: x(:,:) !! input samples
        complex(dp), intent(in), optional :: zi(:,:) !! initial complex delay states
        complex(dp), allocatable, intent(out), optional :: zf(:,:) !! final complex delay states
        integer, intent(out), optional :: status !! success or invalid-input status
        integer, intent(in), optional :: axis !! filter axis, default 2
        complex(dp), allocatable :: y(:,:)

        ! Promote mixed inputs to a common double-complex working type.
        y = lfilter_matrix_complex(cmplx(b,0.0_dp,dp),a, &
            cmplx(x,0.0_dp,dp),zi=zi,zf=zf,status=status,axis=axis)
    end function lfilter_mixed_rcr_2d

    function lfilter_mixed_rcc_2d(b, a, x, zi, zf, status, axis) result(y)
        real(dp), intent(in) :: b(:) !! numerator coefficients
        complex(dp), intent(in) :: a(:) !! denominator coefficients
        complex(dp), intent(in) :: x(:,:) !! input samples
        complex(dp), intent(in), optional :: zi(:,:) !! initial complex delay states
        complex(dp), allocatable, intent(out), optional :: zf(:,:) !! final complex delay states
        integer, intent(out), optional :: status !! success or invalid-input status
        integer, intent(in), optional :: axis !! filter axis, default 2
        complex(dp), allocatable :: y(:,:)

        ! Promote mixed inputs to a common double-complex working type.
        y = lfilter_matrix_complex(cmplx(b,0.0_dp,dp),a, &
            x,zi=zi,zf=zf,status=status,axis=axis)
    end function lfilter_mixed_rcc_2d

    function lfilter_mixed_crr_2d(b, a, x, zi, zf, status, axis) result(y)
        complex(dp), intent(in) :: b(:) !! numerator coefficients
        real(dp), intent(in) :: a(:) !! denominator coefficients
        real(dp), intent(in) :: x(:,:) !! input samples
        complex(dp), intent(in), optional :: zi(:,:) !! initial complex delay states
        complex(dp), allocatable, intent(out), optional :: zf(:,:) !! final complex delay states
        integer, intent(out), optional :: status !! success or invalid-input status
        integer, intent(in), optional :: axis !! filter axis, default 2
        complex(dp), allocatable :: y(:,:)

        ! Promote mixed inputs to a common double-complex working type.
        y = lfilter_matrix_complex(b,cmplx(a,0.0_dp,dp), &
            cmplx(x,0.0_dp,dp),zi=zi,zf=zf,status=status,axis=axis)
    end function lfilter_mixed_crr_2d

    function lfilter_mixed_crc_2d(b, a, x, zi, zf, status, axis) result(y)
        complex(dp), intent(in) :: b(:) !! numerator coefficients
        real(dp), intent(in) :: a(:) !! denominator coefficients
        complex(dp), intent(in) :: x(:,:) !! input samples
        complex(dp), intent(in), optional :: zi(:,:) !! initial complex delay states
        complex(dp), allocatable, intent(out), optional :: zf(:,:) !! final complex delay states
        integer, intent(out), optional :: status !! success or invalid-input status
        integer, intent(in), optional :: axis !! filter axis, default 2
        complex(dp), allocatable :: y(:,:)

        ! Promote mixed inputs to a common double-complex working type.
        y = lfilter_matrix_complex(b,cmplx(a,0.0_dp,dp), &
            x,zi=zi,zf=zf,status=status,axis=axis)
    end function lfilter_mixed_crc_2d

    function lfilter_mixed_ccr_2d(b, a, x, zi, zf, status, axis) result(y)
        complex(dp), intent(in) :: b(:) !! numerator coefficients
        complex(dp), intent(in) :: a(:) !! denominator coefficients
        real(dp), intent(in) :: x(:,:) !! input samples
        complex(dp), intent(in), optional :: zi(:,:) !! initial complex delay states
        complex(dp), allocatable, intent(out), optional :: zf(:,:) !! final complex delay states
        integer, intent(out), optional :: status !! success or invalid-input status
        integer, intent(in), optional :: axis !! filter axis, default 2
        complex(dp), allocatable :: y(:,:)

        ! Promote mixed inputs to a common double-complex working type.
        y = lfilter_matrix_complex(b,a, &
            cmplx(x,0.0_dp,dp),zi=zi,zf=zf,status=status,axis=axis)
    end function lfilter_mixed_ccr_2d

    function lfilter_mixed_rrc_3d(b, a, x, zi, zf, status, axis) result(y)
        real(dp), intent(in) :: b(:) !! numerator coefficients
        real(dp), intent(in) :: a(:) !! denominator coefficients
        complex(dp), intent(in) :: x(:,:,:) !! input samples
        complex(dp), intent(in), optional :: zi(:,:,:) !! initial complex delay states
        complex(dp), allocatable, intent(out), optional :: zf(:,:,:) !! final complex delay states
        integer, intent(out), optional :: status !! success or invalid-input status
        integer, intent(in), optional :: axis !! filter axis, default 3
        complex(dp), allocatable :: y(:,:,:)

        ! Promote mixed inputs to a common double-complex working type.
        y = lfilter_rank3_complex(cmplx(b,0.0_dp,dp),cmplx(a,0.0_dp,dp), &
            x,zi=zi,zf=zf,status=status,axis=axis)
    end function lfilter_mixed_rrc_3d

    function lfilter_mixed_rcr_3d(b, a, x, zi, zf, status, axis) result(y)
        real(dp), intent(in) :: b(:) !! numerator coefficients
        complex(dp), intent(in) :: a(:) !! denominator coefficients
        real(dp), intent(in) :: x(:,:,:) !! input samples
        complex(dp), intent(in), optional :: zi(:,:,:) !! initial complex delay states
        complex(dp), allocatable, intent(out), optional :: zf(:,:,:) !! final complex delay states
        integer, intent(out), optional :: status !! success or invalid-input status
        integer, intent(in), optional :: axis !! filter axis, default 3
        complex(dp), allocatable :: y(:,:,:)

        ! Promote mixed inputs to a common double-complex working type.
        y = lfilter_rank3_complex(cmplx(b,0.0_dp,dp),a, &
            cmplx(x,0.0_dp,dp),zi=zi,zf=zf,status=status,axis=axis)
    end function lfilter_mixed_rcr_3d

    function lfilter_mixed_rcc_3d(b, a, x, zi, zf, status, axis) result(y)
        real(dp), intent(in) :: b(:) !! numerator coefficients
        complex(dp), intent(in) :: a(:) !! denominator coefficients
        complex(dp), intent(in) :: x(:,:,:) !! input samples
        complex(dp), intent(in), optional :: zi(:,:,:) !! initial complex delay states
        complex(dp), allocatable, intent(out), optional :: zf(:,:,:) !! final complex delay states
        integer, intent(out), optional :: status !! success or invalid-input status
        integer, intent(in), optional :: axis !! filter axis, default 3
        complex(dp), allocatable :: y(:,:,:)

        ! Promote mixed inputs to a common double-complex working type.
        y = lfilter_rank3_complex(cmplx(b,0.0_dp,dp),a, &
            x,zi=zi,zf=zf,status=status,axis=axis)
    end function lfilter_mixed_rcc_3d

    function lfilter_mixed_crr_3d(b, a, x, zi, zf, status, axis) result(y)
        complex(dp), intent(in) :: b(:) !! numerator coefficients
        real(dp), intent(in) :: a(:) !! denominator coefficients
        real(dp), intent(in) :: x(:,:,:) !! input samples
        complex(dp), intent(in), optional :: zi(:,:,:) !! initial complex delay states
        complex(dp), allocatable, intent(out), optional :: zf(:,:,:) !! final complex delay states
        integer, intent(out), optional :: status !! success or invalid-input status
        integer, intent(in), optional :: axis !! filter axis, default 3
        complex(dp), allocatable :: y(:,:,:)

        ! Promote mixed inputs to a common double-complex working type.
        y = lfilter_rank3_complex(b,cmplx(a,0.0_dp,dp), &
            cmplx(x,0.0_dp,dp),zi=zi,zf=zf,status=status,axis=axis)
    end function lfilter_mixed_crr_3d

    function lfilter_mixed_crc_3d(b, a, x, zi, zf, status, axis) result(y)
        complex(dp), intent(in) :: b(:) !! numerator coefficients
        real(dp), intent(in) :: a(:) !! denominator coefficients
        complex(dp), intent(in) :: x(:,:,:) !! input samples
        complex(dp), intent(in), optional :: zi(:,:,:) !! initial complex delay states
        complex(dp), allocatable, intent(out), optional :: zf(:,:,:) !! final complex delay states
        integer, intent(out), optional :: status !! success or invalid-input status
        integer, intent(in), optional :: axis !! filter axis, default 3
        complex(dp), allocatable :: y(:,:,:)

        ! Promote mixed inputs to a common double-complex working type.
        y = lfilter_rank3_complex(b,cmplx(a,0.0_dp,dp), &
            x,zi=zi,zf=zf,status=status,axis=axis)
    end function lfilter_mixed_crc_3d

    function lfilter_mixed_ccr_3d(b, a, x, zi, zf, status, axis) result(y)
        complex(dp), intent(in) :: b(:) !! numerator coefficients
        complex(dp), intent(in) :: a(:) !! denominator coefficients
        real(dp), intent(in) :: x(:,:,:) !! input samples
        complex(dp), intent(in), optional :: zi(:,:,:) !! initial complex delay states
        complex(dp), allocatable, intent(out), optional :: zf(:,:,:) !! final complex delay states
        integer, intent(out), optional :: status !! success or invalid-input status
        integer, intent(in), optional :: axis !! filter axis, default 3
        complex(dp), allocatable :: y(:,:,:)

        ! Promote mixed inputs to a common double-complex working type.
        y = lfilter_rank3_complex(b,a, &
            cmplx(x,0.0_dp,dp),zi=zi,zf=zf,status=status,axis=axis)
    end function lfilter_mixed_ccr_3d

    ! Each lane uses the same normalized coefficients and a private delay line.
    ! Strided array sections permit filtering without repacking the full array.
    subroutine filter_line_real(x, y, b, a, state)
        real(dp), intent(in) :: x(:) !! samples from one array channel
        real(dp), intent(out) :: y(:) !! filtered samples
        real(dp), intent(in) :: b(:) !! normalized padded numerator
        real(dp), intent(in) :: a(:) !! normalized padded denominator
        real(dp), intent(inout) :: state(:) !! transposed DF-II delay line
        real(dp) :: value
        integer :: k, j, ns

        ns = size(state)
        do k = 1, size(x)
            value = b(1)*x(k)
            if (ns > 0) value = value + state(1)
            y(k) = value
            do j = 1, ns - 1
                state(j) = state(j+1) + b(j+1)*x(k) - a(j+1)*value
            end do
            if (ns > 0) state(ns) = b(ns+1)*x(k) - a(ns+1)*value
        end do
    end subroutine filter_line_real

    subroutine filter_line_complex(x, y, b, a, state)
        complex(dp), intent(in) :: x(:) !! samples from one array channel
        complex(dp), intent(out) :: y(:) !! filtered samples
        complex(dp), intent(in) :: b(:) !! normalized padded numerator
        complex(dp), intent(in) :: a(:) !! normalized padded denominator
        complex(dp), intent(inout) :: state(:) !! transposed DF-II delay line
        complex(dp) :: value
        integer :: k, j, ns

        ns = size(state)
        do k = 1, size(x)
            value = b(1)*x(k)
            if (ns > 0) value = value + state(1)
            y(k) = value
            do j = 1, ns - 1
                state(j) = state(j+1) + b(j+1)*x(k) - a(j+1)*value
            end do
            if (ns > 0) state(ns) = b(ns+1)*x(k) - a(ns+1)*value
        end do
    end subroutine filter_line_complex

    ! CONTIGUOUS assumed-rank dummies are copied in/out by the compiler when
    ! the actual arguments are array sections.  The C intrinsics only create
    ! rank-one views of these buffers; they do not invoke external C code.
    subroutine lfilter_nd_real(b, a, x, y, zi, zf, status, axis)
        real(dp), intent(in) :: b(:) !! numerator coefficients, lowest lag first
        real(dp), intent(in) :: a(:) !! denominator coefficients, a(1) nonzero
        real(dp), contiguous, target, intent(in) :: x(..) !! input rank 1 through 15
        real(dp), contiguous, target, intent(out) :: y(..) !! caller-allocated output, shape(x)
        real(dp), contiguous, target, intent(in), optional :: zi(..) !! initial states, shape(x) with axis length ns
        real(dp), contiguous, target, intent(out), optional :: zf(..) !! final states, same shape as zi
        integer, intent(out), optional :: status !! success or invalid-input status
        integer, intent(in), optional :: axis !! one-based axis, defaults to last
        real(dp), pointer :: xin(:) => null(), yout(:) => null()
        real(dp), pointer :: state_in(:) => null(), state_out(:) => null()
        real(dp), allocatable :: aa(:), bb(:), state(:)
        real(dp) :: nanval
        integer :: ax, ns, dims(rank(x)), state_dims(rank(x)), k
        integer(int64) :: stride, block, state_block, outer, i, j, first, last, first_state

        nanval = ieee_value(0.0_dp, ieee_quiet_nan)
        if (present(status)) status = lfilter_status_invalid_input
        if (size(y,kind=int64) > 0_int64) then
            call c_f_pointer(c_loc(y), yout, [size(y,kind=int64)])
            yout = nanval
        end if
        if (present(zf)) then
            if (size(zf,kind=int64) > 0_int64) then
                call c_f_pointer(c_loc(zf), state_out, [size(zf,kind=int64)])
                state_out = nanval
            end if
        end if
        if (rank(x) < 1 .or. rank(x) /= rank(y)) return
        if (any(shape(x) /= shape(y))) return
        ax = rank(x)
        if (present(axis)) ax = axis
        if (ax < 1 .or. ax > rank(x)) return
        if (size(a) < 1 .or. size(b) < 1) return
        if (.not. ieee_is_finite(a(1))) return
        if (abs(a(1)) < tiny(1.0_dp)) return
        if (.not. (all(ieee_is_finite(a)) .and. all(ieee_is_finite(b)))) return

        ns = max(size(a),size(b))-1
        dims = shape(x)
        state_dims = dims
        state_dims(ax) = ns
        if (present(zi)) then
            if (rank(zi) /= rank(x)) return
            if (any(shape(zi) /= state_dims)) return
            if (size(zi,kind=int64) > 0_int64) then
                call c_f_pointer(c_loc(zi), state_in, [size(zi,kind=int64)])
                if (.not. all(ieee_is_finite(state_in))) return
            end if
        end if
        if (present(zf)) then
            if (rank(zf) /= rank(x)) return
            if (any(shape(zf) /= state_dims)) return
        end if
        allocate(aa(ns+1),bb(ns+1),state(ns))
        aa = 0.0_dp
        bb = 0.0_dp
        aa(1:size(a)) = a/a(1)
        bb(1:size(b)) = b/a(1)
        if (.not. (all(ieee_is_finite(aa)) .and. all(ieee_is_finite(bb)))) return

        if (size(x,kind=int64) > 0_int64) call c_f_pointer(c_loc(x), xin, [size(x,kind=int64)])
        stride = product(int(dims(:ax-1),int64))
        block = stride*int(dims(ax),int64)
        state_block = stride*int(ns,int64)
        outer = product(int(dims(ax+1:),int64))
        do j = 0_int64, outer-1_int64
            do i = 0_int64, stride-1_int64
                state = 0.0_dp
                first_state = j*state_block+i+1_int64
                if (present(zi)) then
                    do k = 1, ns
                        state(k) = state_in(first_state+(k-1)*stride)
                    end do
                end if
                if (dims(ax) > 0) then
                    first = j*block+i+1_int64
                    last = first+int(dims(ax)-1,int64)*stride
                    call filter_line_real(xin(first:last:stride),yout(first:last:stride),bb,aa,state)
                end if
                if (present(zf)) then
                    do k = 1, ns
                        state_out(first_state+(k-1)*stride) = state(k)
                    end do
                end if
            end do
        end do
        if (present(status)) status = lfilter_status_success
    end subroutine lfilter_nd_real

    subroutine lfilter_nd_complex(b, a, x, y, zi, zf, status, axis)
        complex(dp), intent(in) :: b(:) !! numerator coefficients, lowest lag first
        complex(dp), intent(in) :: a(:) !! denominator coefficients, a(1) nonzero
        complex(dp), contiguous, target, intent(in) :: x(..) !! input rank 1 through 15
        complex(dp), contiguous, target, intent(out) :: y(..) !! caller-allocated output, shape(x)
        complex(dp), contiguous, target, intent(in), optional :: zi(..) !! initial states, shape(x) with axis length ns
        complex(dp), contiguous, target, intent(out), optional :: zf(..) !! final states, same shape as zi
        integer, intent(out), optional :: status !! success or invalid-input status
        integer, intent(in), optional :: axis !! one-based axis, defaults to last
        complex(dp), pointer :: xin(:) => null(), yout(:) => null()
        complex(dp), pointer :: state_in(:) => null(), state_out(:) => null()
        complex(dp), allocatable :: aa(:), bb(:), state(:)
        real(dp) :: nanval
        integer :: ax, ns, dims(rank(x)), state_dims(rank(x)), k
        integer(int64) :: stride, block, state_block, outer, i, j, first, last, first_state

        nanval = ieee_value(0.0_dp, ieee_quiet_nan)
        if (present(status)) status = lfilter_status_invalid_input
        if (size(y,kind=int64) > 0_int64) then
            call c_f_pointer(c_loc(y), yout, [size(y,kind=int64)])
            yout = cmplx(nanval,nanval,dp)
        end if
        if (present(zf)) then
            if (size(zf,kind=int64) > 0_int64) then
                call c_f_pointer(c_loc(zf), state_out, [size(zf,kind=int64)])
                state_out = cmplx(nanval,nanval,dp)
            end if
        end if
        if (rank(x) < 1 .or. rank(x) /= rank(y)) return
        if (any(shape(x) /= shape(y))) return
        ax = rank(x)
        if (present(axis)) ax = axis
        if (ax < 1 .or. ax > rank(x)) return
        if (size(a) < 1 .or. size(b) < 1) return
        if (.not. (ieee_is_finite(real(a(1),dp)) .and. ieee_is_finite(aimag(a(1))))) return
        if (abs(a(1)) < tiny(1.0_dp)) return
        if (.not. (all(ieee_is_finite(real(a,dp))) .and. all(ieee_is_finite(aimag(a))) .and. &
            all(ieee_is_finite(real(b,dp))) .and. all(ieee_is_finite(aimag(b))))) return

        ns = max(size(a),size(b))-1
        dims = shape(x)
        state_dims = dims
        state_dims(ax) = ns
        if (present(zi)) then
            if (rank(zi) /= rank(x)) return
            if (any(shape(zi) /= state_dims)) return
            if (size(zi,kind=int64) > 0_int64) then
                call c_f_pointer(c_loc(zi), state_in, [size(zi,kind=int64)])
                if (.not. (all(ieee_is_finite(real(state_in,dp))) .and. &
                    all(ieee_is_finite(aimag(state_in))))) return
            end if
        end if
        if (present(zf)) then
            if (rank(zf) /= rank(x)) return
            if (any(shape(zf) /= state_dims)) return
        end if
        allocate(aa(ns+1),bb(ns+1),state(ns))
        aa = (0.0_dp,0.0_dp)
        bb = (0.0_dp,0.0_dp)
        aa(1:size(a)) = a/a(1)
        bb(1:size(b)) = b/a(1)
        if (.not. (all(ieee_is_finite(real(aa,dp))) .and. all(ieee_is_finite(aimag(aa))) .and. &
            all(ieee_is_finite(real(bb,dp))) .and. all(ieee_is_finite(aimag(bb))))) return

        if (size(x,kind=int64) > 0_int64) call c_f_pointer(c_loc(x), xin, [size(x,kind=int64)])
        stride = product(int(dims(:ax-1),int64))
        block = stride*int(dims(ax),int64)
        state_block = stride*int(ns,int64)
        outer = product(int(dims(ax+1:),int64))
        do j = 0_int64, outer-1_int64
            do i = 0_int64, stride-1_int64
                state = (0.0_dp,0.0_dp)
                first_state = j*state_block+i+1_int64
                if (present(zi)) then
                    do k = 1, ns
                        state(k) = state_in(first_state+(k-1)*stride)
                    end do
                end if
                if (dims(ax) > 0) then
                    first = j*block+i+1_int64
                    last = first+int(dims(ax)-1,int64)*stride
                    call filter_line_complex(xin(first:last:stride),yout(first:last:stride),bb,aa,state)
                end if
                if (present(zf)) then
                    do k = 1, ns
                        state_out(first_state+(k-1)*stride) = state(k)
                    end do
                end if
            end do
        end do
        if (present(status)) status = lfilter_status_success
    end subroutine lfilter_nd_complex


    ! A real input with complex coefficients produces a complex output and
    ! complex filter states.  One complex work vector is reused per channel.
    subroutine lfilter_nd_complex_real_input(b, a, x, y, zi, zf, status, axis)
        complex(dp), intent(in) :: b(:) !! complex numerator coefficients
        complex(dp), intent(in) :: a(:) !! complex denominator coefficients
        real(dp), contiguous, target, intent(in) :: x(..) !! input ranks 1 through 15
        complex(dp), contiguous, target, intent(out) :: y(..) !! complex output, same shape as x
        complex(dp), contiguous, target, intent(in), optional :: zi(..) !! initial states along axis
        complex(dp), contiguous, target, intent(out), optional :: zf(..) !! final states along axis
        integer, intent(out), optional :: status !! success or invalid-input status
        integer, intent(in), optional :: axis !! one-based filter axis, default last
        real(dp), pointer :: xin(:) => null()
        complex(dp), pointer :: yout(:) => null()
        complex(dp), pointer :: state_in(:) => null(), state_out(:) => null()
        complex(dp), allocatable :: aa(:), bb(:), state(:), line(:)
        real(dp) :: nanval
        integer :: ax, ns, dims(rank(x)), state_dims(rank(x)), k
        integer(int64) :: stride, block, state_block, outer, i, j, first, last, first_state

        nanval = ieee_value(0.0_dp, ieee_quiet_nan)
        if (present(status)) status = lfilter_status_invalid_input
        if (size(y,kind=int64) > 0_int64) then
            call c_f_pointer(c_loc(y), yout, [size(y,kind=int64)])
            yout = cmplx(nanval,nanval,dp)
        end if
        if (present(zf)) then
            if (size(zf,kind=int64) > 0_int64) then
                call c_f_pointer(c_loc(zf), state_out, [size(zf,kind=int64)])
                state_out = cmplx(nanval,nanval,dp)
            end if
        end if
        if (rank(x) < 1 .or. rank(x) /= rank(y)) return
        if (any(shape(x) /= shape(y))) return
        ax = rank(x)
        if (present(axis)) ax = axis
        if (ax < 1 .or. ax > rank(x)) return
        if (size(a) < 1 .or. size(b) < 1) return
        if (.not. (ieee_is_finite(real(a(1),dp)) .and. ieee_is_finite(aimag(a(1))))) return
        if (abs(a(1)) < tiny(1.0_dp)) return
        if (.not. (all(ieee_is_finite(real(a,dp))) .and. all(ieee_is_finite(aimag(a))) .and. &
            all(ieee_is_finite(real(b,dp))) .and. all(ieee_is_finite(aimag(b))))) return

        ns = max(size(a),size(b))-1
        dims = shape(x)
        state_dims = dims
        state_dims(ax) = ns
        if (present(zi)) then
            if (rank(zi) /= rank(x)) return
            if (any(shape(zi) /= state_dims)) return
            if (size(zi,kind=int64) > 0_int64) then
                call c_f_pointer(c_loc(zi), state_in, [size(zi,kind=int64)])
                if (.not. (all(ieee_is_finite(real(state_in,dp))) .and. &
                    all(ieee_is_finite(aimag(state_in))))) return
            end if
        end if
        if (present(zf)) then
            if (rank(zf) /= rank(x)) return
            if (any(shape(zf) /= state_dims)) return
        end if
        allocate(aa(ns+1),bb(ns+1),state(ns),line(dims(ax)))
        aa = (0.0_dp,0.0_dp)
        bb = (0.0_dp,0.0_dp)
        aa(1:size(a)) = a/a(1)
        bb(1:size(b)) = b/a(1)
        if (.not. (all(ieee_is_finite(real(aa,dp))) .and. all(ieee_is_finite(aimag(aa))) .and. &
            all(ieee_is_finite(real(bb,dp))) .and. all(ieee_is_finite(aimag(bb))))) return

        if (size(x,kind=int64) > 0_int64) call c_f_pointer(c_loc(x), xin, [size(x,kind=int64)])
        stride = product(int(dims(:ax-1),int64))
        block = stride*int(dims(ax),int64)
        state_block = stride*int(ns,int64)
        outer = product(int(dims(ax+1:),int64))
        do j = 0_int64, outer-1_int64
            do i = 0_int64, stride-1_int64
                state = (0.0_dp,0.0_dp)
                first_state = j*state_block+i+1_int64
                if (present(zi)) then
                    do k = 1, ns
                        state(k) = state_in(first_state+(k-1)*stride)
                    end do
                end if
                if (dims(ax) > 0) then
                    first = j*block+i+1_int64
                    last = first+int(dims(ax)-1,int64)*stride
                    line = cmplx(xin(first:last:stride),0.0_dp,dp)
                    call filter_line_complex(line,yout(first:last:stride),bb,aa,state)
                end if
                if (present(zf)) then
                    do k = 1, ns
                        state_out(first_state+(k-1)*stride) = state(k)
                    end do
                end if
            end do
        end do
        if (present(status)) status = lfilter_status_success
    end subroutine lfilter_nd_complex_real_input

    subroutine lfilter_nd_mixed_rrc(b, a, x, y, zi, zf, status, axis)
        real(dp), intent(in) :: b(:) !! numerator filter coefficients
        real(dp), intent(in) :: a(:) !! denominator filter coefficients
        complex(dp), contiguous, target, intent(in) :: x(..) !! input array up to rank 15
        complex(dp), contiguous, target, intent(out) :: y(..) !! complex output, shape(x)
        complex(dp), contiguous, target, intent(in), optional :: zi(..) !! initial states
        complex(dp), contiguous, target, intent(out), optional :: zf(..) !! final states
        integer, intent(out), optional :: status !! success or invalid-input status
        integer, intent(in), optional :: axis !! one-based filter axis, default last

        call lfilter_nd_complex(cmplx(b,0.0_dp,dp),cmplx(a,0.0_dp,dp), &
            x,y,zi=zi,zf=zf,status=status,axis=axis)
    end subroutine lfilter_nd_mixed_rrc


    subroutine lfilter_nd_mixed_rcr(b, a, x, y, zi, zf, status, axis)
        real(dp), intent(in) :: b(:) !! numerator filter coefficients
        complex(dp), intent(in) :: a(:) !! denominator filter coefficients
        real(dp), contiguous, target, intent(in) :: x(..) !! input array up to rank 15
        complex(dp), contiguous, target, intent(out) :: y(..) !! complex output, shape(x)
        complex(dp), contiguous, target, intent(in), optional :: zi(..) !! initial states
        complex(dp), contiguous, target, intent(out), optional :: zf(..) !! final states
        integer, intent(out), optional :: status !! success or invalid-input status
        integer, intent(in), optional :: axis !! one-based filter axis, default last

        call lfilter_nd_complex_real_input(cmplx(b,0.0_dp,dp),a, &
            x,y,zi=zi,zf=zf,status=status,axis=axis)
    end subroutine lfilter_nd_mixed_rcr


    subroutine lfilter_nd_mixed_rcc(b, a, x, y, zi, zf, status, axis)
        real(dp), intent(in) :: b(:) !! numerator filter coefficients
        complex(dp), intent(in) :: a(:) !! denominator filter coefficients
        complex(dp), contiguous, target, intent(in) :: x(..) !! input array up to rank 15
        complex(dp), contiguous, target, intent(out) :: y(..) !! complex output, shape(x)
        complex(dp), contiguous, target, intent(in), optional :: zi(..) !! initial states
        complex(dp), contiguous, target, intent(out), optional :: zf(..) !! final states
        integer, intent(out), optional :: status !! success or invalid-input status
        integer, intent(in), optional :: axis !! one-based filter axis, default last

        call lfilter_nd_complex(cmplx(b,0.0_dp,dp),a, &
            x,y,zi=zi,zf=zf,status=status,axis=axis)
    end subroutine lfilter_nd_mixed_rcc


    subroutine lfilter_nd_mixed_crr(b, a, x, y, zi, zf, status, axis)
        complex(dp), intent(in) :: b(:) !! numerator filter coefficients
        real(dp), intent(in) :: a(:) !! denominator filter coefficients
        real(dp), contiguous, target, intent(in) :: x(..) !! input array up to rank 15
        complex(dp), contiguous, target, intent(out) :: y(..) !! complex output, shape(x)
        complex(dp), contiguous, target, intent(in), optional :: zi(..) !! initial states
        complex(dp), contiguous, target, intent(out), optional :: zf(..) !! final states
        integer, intent(out), optional :: status !! success or invalid-input status
        integer, intent(in), optional :: axis !! one-based filter axis, default last

        call lfilter_nd_complex_real_input(b,cmplx(a,0.0_dp,dp), &
            x,y,zi=zi,zf=zf,status=status,axis=axis)
    end subroutine lfilter_nd_mixed_crr


    subroutine lfilter_nd_mixed_crc(b, a, x, y, zi, zf, status, axis)
        complex(dp), intent(in) :: b(:) !! numerator filter coefficients
        real(dp), intent(in) :: a(:) !! denominator filter coefficients
        complex(dp), contiguous, target, intent(in) :: x(..) !! input array up to rank 15
        complex(dp), contiguous, target, intent(out) :: y(..) !! complex output, shape(x)
        complex(dp), contiguous, target, intent(in), optional :: zi(..) !! initial states
        complex(dp), contiguous, target, intent(out), optional :: zf(..) !! final states
        integer, intent(out), optional :: status !! success or invalid-input status
        integer, intent(in), optional :: axis !! one-based filter axis, default last

        call lfilter_nd_complex(b,cmplx(a,0.0_dp,dp), &
            x,y,zi=zi,zf=zf,status=status,axis=axis)
    end subroutine lfilter_nd_mixed_crc


    subroutine lfilter_nd_mixed_ccr(b, a, x, y, zi, zf, status, axis)
        complex(dp), intent(in) :: b(:) !! numerator filter coefficients
        complex(dp), intent(in) :: a(:) !! denominator filter coefficients
        real(dp), contiguous, target, intent(in) :: x(..) !! input array up to rank 15
        complex(dp), contiguous, target, intent(out) :: y(..) !! complex output, shape(x)
        complex(dp), contiguous, target, intent(in), optional :: zi(..) !! initial states
        complex(dp), contiguous, target, intent(out), optional :: zf(..) !! final states
        integer, intent(out), optional :: status !! success or invalid-input status
        integer, intent(in), optional :: axis !! one-based filter axis, default last

        call lfilter_nd_complex_real_input(b,a, &
            x,y,zi=zi,zf=zf,status=status,axis=axis)
    end subroutine lfilter_nd_mixed_ccr

end module scifort_signal_lfilter
