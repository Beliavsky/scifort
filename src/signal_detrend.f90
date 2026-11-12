! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Independent least-squares trend removal on real and complex arrays.
module scifort_signal_detrend
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    implicit none
    private

    integer, parameter, public :: detrend_status_success = 0
    integer, parameter, public :: detrend_status_invalid_input = 1
    interface detrend
        module procedure detrend_real_1d, detrend_complex_1d
        module procedure detrend_real_2d, detrend_complex_2d
        module procedure detrend_real_3d, detrend_complex_3d
    end interface
    public :: detrend

contains

    function detrend_real_1d(x, kind, breaks, status) result(y)
        real(dp), intent(in) :: x(:) !! finite 1D observations
        character(len=*), intent(in), optional :: kind !! 'linear' (default) or 'constant'
        integer, intent(in), optional :: breaks(:) !! one-based segment starts in [1,n+1]
        integer, intent(out), optional :: status !! zero if valid, one otherwise
        real(dp), allocatable :: y(:)

        integer :: n, i, j, k, a, b, m, nbreak, item
        integer, allocatable :: cuts(:)
        real(dp) :: avg, slope, center, denom, numerator
        logical :: linear

        n = size(x)
        allocate(y(n))
        y = ieee_value(0.0_dp, ieee_quiet_nan)
        if (present(status)) status = detrend_status_invalid_input
        if (any(.not. ieee_is_finite(x))) return
        linear = .true.
        if (present(kind)) then
            select case (kind)
            case ('linear')
                linear = .true.
            case ('constant')
                linear = .false.
            case default
                return
            end select
        end if
        if (.not. linear) then
            y = x
            if (n > 0) y = x - sum(x) / real(n, dp)
            if (present(status)) status = detrend_status_success
            return
        end if
        nbreak = 0
        if (present(breaks)) nbreak = size(breaks)
        allocate(cuts(nbreak + 2))
        cuts(1) = 1
        cuts(2) = n + 1
        if (present(breaks)) then
            do i = 1, nbreak
                if (breaks(i) < 1 .or. breaks(i) > n + 1) return
                cuts(i + 2) = breaks(i)
            end do
        end if
        ! Stable insertion sort; breakpoints are normally a short list.
        do i = 2, size(cuts)
            item = cuts(i)
            j = i - 1
            do while (j >= 1)
                if (cuts(j) <= item) exit
                cuts(j + 1) = cuts(j)
                j = j - 1
            end do
            cuts(j + 1) = item
        end do
        do k = 1, size(cuts) - 1
            a = cuts(k)
            b = cuts(k + 1) - 1
            m = b - a + 1
            if (m < 1) cycle
            avg = sum(x(a:b)) / real(m, dp)
            center = 0.5_dp * real(m + 1, dp)
            numerator = 0.0_dp
            denom = 0.0_dp
            do i = 1, m
                slope = real(i, dp) - center
                numerator = numerator + slope * (x(a + i - 1) - avg)
                denom = denom + slope * slope
            end do
            slope = 0.0_dp
            if (denom > 0.0_dp) slope = numerator / denom
            do i = 1, m
                y(a + i - 1) = x(a + i - 1) - (avg + slope * (real(i, dp) - center))
            end do
        end do
        if (present(status)) status = detrend_status_success
    end function detrend_real_1d


    ! Breakpoint conventions are intentionally one-based, matching the
    ! original SciFort 1D interface. SciPy's zero-based bp needs +1.
    logical function valid_options(n,kind,breaks) result(ok)
        integer, intent(in) :: n !! length along the selected axis
        character(*), intent(in), optional :: kind !! linear or constant
        integer, intent(in), optional :: breaks(:) !! one-based segment starts
        logical :: linear

        ok = .true.
        linear = .true.
        if (present(kind)) then
            if (kind /= 'linear' .and. kind /= 'constant') ok = .false.
            linear = kind == 'linear'
        end if
        if (linear .and. present(breaks)) then
            if (any(breaks < 1) .or. any(breaks > n+1)) ok = .false.
        end if
    end function valid_options

    function detrend_complex_1d(x,kind,breaks,status) result(y)
        complex(dp), intent(in) :: x(:) !! observations
        character(*), intent(in), optional :: kind !! linear or constant
        integer, intent(in), optional :: breaks(:) !! one-based breakpoints
        integer, intent(out), optional :: status !! zero on success
        complex(dp), allocatable :: y(:)
        real(dp), allocatable :: yr(:), yi(:)
        integer :: ierr

        allocate(y(size(x)))
        y = cmplx(ieee_value(0.0_dp,ieee_quiet_nan), &
            ieee_value(0.0_dp,ieee_quiet_nan),dp)
        if (present(status)) status = detrend_status_invalid_input
        if (.not. all(ieee_is_finite(real(x,dp)))) return
        if (.not. all(ieee_is_finite(aimag(x)))) return
        yr = detrend_real_1d(real(x,dp),kind,breaks,ierr)
        if (ierr /= detrend_status_success) return
        yi = detrend_real_1d(aimag(x),kind,breaks,ierr)
        if (ierr /= detrend_status_success) return
        y = cmplx(yr,yi,dp)
        if (present(status)) status = detrend_status_success
    end function detrend_complex_1d


    function detrend_real_2d(x,kind,breaks,axis,status) result(y)
        real(dp), intent(in) :: x(:,:) !! data array
        character(*), intent(in), optional :: kind !! linear (default) or constant
        integer, intent(in), optional :: breaks(:) !! one-based segment starts
        integer, intent(in), optional :: axis !! one-based axis, default 2
        integer, intent(out), optional :: status !! zero on success
        real(dp), allocatable :: y(:,:)
        integer :: dim, i,j

        allocate(y, mold=x)
        y = ieee_value(0.0_dp,ieee_quiet_nan)
        if (present(status)) status = detrend_status_invalid_input
        dim = 2
        if (present(axis)) dim = axis
        if (dim < 1 .or. dim > 2) return
        if (.not. valid_options(size(x,dim),kind,breaks)) return
        if (.not. all(ieee_is_finite(x))) return
        select case (dim)
        case (1)
            do j = 1, size(x,2)
                y(:,j) = detrend_real_1d(x(:,j),kind,breaks)
            end do
        case (2)
            do i = 1, size(x,1)
                y(i,:) = detrend_real_1d(x(i,:),kind,breaks)
            end do
        end select
        if (present(status)) status = detrend_status_success
    end function detrend_real_2d

    function detrend_complex_2d(x,kind,breaks,axis,status) result(y)
        complex(dp), intent(in) :: x(:,:) !! data array
        character(*), intent(in), optional :: kind !! linear (default) or constant
        integer, intent(in), optional :: breaks(:) !! one-based segment starts
        integer, intent(in), optional :: axis !! one-based axis, default 2
        integer, intent(out), optional :: status !! zero on success
        complex(dp), allocatable :: y(:,:)
        integer :: dim, i,j

        allocate(y, mold=x)
        y = cmplx(ieee_value(0.0_dp,ieee_quiet_nan), &
            ieee_value(0.0_dp,ieee_quiet_nan),dp)
        if (present(status)) status = detrend_status_invalid_input
        dim = 2
        if (present(axis)) dim = axis
        if (dim < 1 .or. dim > 2) return
        if (.not. valid_options(size(x,dim),kind,breaks)) return
        if (.not. all(ieee_is_finite(real(x,dp)))) return
        if (.not. all(ieee_is_finite(aimag(x)))) return
        select case (dim)
        case (1)
            do j = 1, size(x,2)
                y(:,j) = detrend_complex_1d(x(:,j),kind,breaks)
            end do
        case (2)
            do i = 1, size(x,1)
                y(i,:) = detrend_complex_1d(x(i,:),kind,breaks)
            end do
        end select
        if (present(status)) status = detrend_status_success
    end function detrend_complex_2d

    function detrend_real_3d(x,kind,breaks,axis,status) result(y)
        real(dp), intent(in) :: x(:,:,:) !! data array
        character(*), intent(in), optional :: kind !! linear (default) or constant
        integer, intent(in), optional :: breaks(:) !! one-based segment starts
        integer, intent(in), optional :: axis !! one-based axis, default 3
        integer, intent(out), optional :: status !! zero on success
        real(dp), allocatable :: y(:,:,:)
        integer :: dim, i,j,k

        allocate(y, mold=x)
        y = ieee_value(0.0_dp,ieee_quiet_nan)
        if (present(status)) status = detrend_status_invalid_input
        dim = 3
        if (present(axis)) dim = axis
        if (dim < 1 .or. dim > 3) return
        if (.not. valid_options(size(x,dim),kind,breaks)) return
        if (.not. all(ieee_is_finite(x))) return
        select case (dim)
        case (1)
            do k = 1, size(x,3)
                do j = 1, size(x,2)
                    y(:,j,k) = detrend_real_1d(x(:,j,k),kind,breaks)
                end do
            end do
        case (2)
            do k = 1, size(x,3)
                do i = 1, size(x,1)
                    y(i,:,k) = detrend_real_1d(x(i,:,k),kind,breaks)
                end do
            end do
        case (3)
            do j = 1, size(x,2)
                do i = 1, size(x,1)
                    y(i,j,:) = detrend_real_1d(x(i,j,:),kind,breaks)
                end do
            end do
        end select
        if (present(status)) status = detrend_status_success
    end function detrend_real_3d

    function detrend_complex_3d(x,kind,breaks,axis,status) result(y)
        complex(dp), intent(in) :: x(:,:,:) !! data array
        character(*), intent(in), optional :: kind !! linear (default) or constant
        integer, intent(in), optional :: breaks(:) !! one-based segment starts
        integer, intent(in), optional :: axis !! one-based axis, default 3
        integer, intent(out), optional :: status !! zero on success
        complex(dp), allocatable :: y(:,:,:)
        integer :: dim, i,j,k

        allocate(y, mold=x)
        y = cmplx(ieee_value(0.0_dp,ieee_quiet_nan), &
            ieee_value(0.0_dp,ieee_quiet_nan),dp)
        if (present(status)) status = detrend_status_invalid_input
        dim = 3
        if (present(axis)) dim = axis
        if (dim < 1 .or. dim > 3) return
        if (.not. valid_options(size(x,dim),kind,breaks)) return
        if (.not. all(ieee_is_finite(real(x,dp)))) return
        if (.not. all(ieee_is_finite(aimag(x)))) return
        select case (dim)
        case (1)
            do k = 1, size(x,3)
                do j = 1, size(x,2)
                    y(:,j,k) = detrend_complex_1d(x(:,j,k),kind,breaks)
                end do
            end do
        case (2)
            do k = 1, size(x,3)
                do i = 1, size(x,1)
                    y(i,:,k) = detrend_complex_1d(x(i,:,k),kind,breaks)
                end do
            end do
        case (3)
            do j = 1, size(x,2)
                do i = 1, size(x,1)
                    y(i,j,:) = detrend_complex_1d(x(i,j,:),kind,breaks)
                end do
            end do
        end select
        if (present(status)) status = detrend_status_success
    end function detrend_complex_3d

end module scifort_signal_detrend
