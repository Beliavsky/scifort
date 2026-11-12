! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Public SciPy-style mirror-symmetric first- and second-order recursive
! smoothers. Reuses the independently implemented recurrences for splines.
module scifort_signal_symiir
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_signal_spline2d, only : symmetric_iir1, symmetric_iir2
    implicit none
    private

    integer, parameter, public :: symiir_status_success = 0
    integer, parameter, public :: symiir_status_invalid_input = 1
    integer, parameter, public :: symiir_status_no_convergence = 2
    interface symiirorder1
        module procedure symiirorder1_real
        module procedure symiirorder1_complex
        module procedure symiirorder1_real2d
        module procedure symiirorder1_complex2d
    end interface
    interface symiirorder2
        module procedure symiirorder2_real
        module procedure symiirorder2_complex
        module procedure symiirorder2_real2d
        module procedure symiirorder2_complex2d
    end interface
    public :: symiirorder1, symiirorder2

contains

    function symiir_kernel(signal, order, first, second, precision, status) result(y)
        complex(dp), intent(in) :: signal(:) !! samples with mirror-symmetric endpoints
        integer, intent(in) :: order !! 1 for real pole, 2 for complex pole pair
        real(dp), intent(in) :: first !! first-order numerator or second-order pole radius
        real(dp), intent(in) :: second !! first-order pole or second-order pole angle
        real(dp), intent(in), optional :: precision !! boundary tolerance, default 1e-11
        integer, intent(out) :: status !! success, invalid input or insufficient tail
        complex(dp), allocatable :: y(:)
        real(dp) :: pole, eps
        integer :: n

        allocate(y(0))
        status = symiir_status_invalid_input
        n = size(signal)
        if (n < 1 .or. .not. ieee_is_finite(first) .or. .not. ieee_is_finite(second)) return
        if (.not. all(ieee_is_finite(real(signal,dp))) .or. &
            .not. all(ieee_is_finite(aimag(signal)))) return
        eps = 1.0e-11_dp
        if (present(precision)) then
            if (.not. ieee_is_finite(precision)) return
            if (precision > 0.0_dp) eps = precision
        end if
        if (eps <= 0.0_dp .or. eps >= 1.0_dp) return
        select case (order)
        case (1)
            pole = second
            if (abs(pole) >= 1.0_dp) return
        case (2)
            pole = first
            if (abs(pole) >= 1.0_dp) return
            if (pole /= 0.0_dp) then
                if (n < 2 .or. abs(sin(second)) < 1.0e-12_dp) return
            end if
        case default
            return
        end select
        ! The underlying boundary series must have decayed before the last
        ! sample; otherwise truncating to the finite signal is unreliable.
        if (abs(pole) > 0.0_dp) then
            if (real(n,dp)*log(abs(pole)) > log(eps)) then
                status = symiir_status_no_convergence
                return
            end if
        end if
        if (order == 1) then
            y = symmetric_iir1(signal,pole,first,eps)
        else if (pole == 0.0_dp) then
            y = signal
        else
            y = symmetric_iir2(signal,pole,second,eps)
        end if
        if (.not. all(ieee_is_finite(real(y,dp))) .or. &
            .not. all(ieee_is_finite(aimag(y)))) then
            deallocate(y)
            allocate(y(0))
            status = symiir_status_no_convergence
            return
        end if
        status = symiir_status_success
    end function symiir_kernel

    function symiirorder1_real(signal, c0, z1, precision, status) result(y)
        real(dp), intent(in) :: signal(:) !! real input samples
        real(dp), intent(in) :: c0 !! symmetric transfer-function numerator
        real(dp), intent(in) :: z1 !! real stable pole
        real(dp), optional, intent(in) :: precision !! truncation tolerance
        integer, optional, intent(out) :: status !! zero on success
        real(dp), allocatable :: y(:)
        complex(dp), allocatable :: z(:)
        integer :: code
        z = symiir_kernel(cmplx(signal,0.0_dp,dp),1,c0,z1,precision,code)
        y = real(z,dp)
        if (present(status)) status = code
    end function symiirorder1_real

    function symiirorder1_complex(signal, c0, z1, precision, status) result(y)
        complex(dp), intent(in) :: signal(:) !! complex input samples
        real(dp), intent(in) :: c0 !! symmetric transfer-function numerator
        real(dp), intent(in) :: z1 !! stable real pole
        real(dp), optional, intent(in) :: precision !! truncation tolerance
        integer, optional, intent(out) :: status !! zero on success
        complex(dp), allocatable :: y(:)
        integer :: code
        y = symiir_kernel(signal,1,c0,z1,precision,code)
        if (present(status)) status = code
    end function symiirorder1_complex

    function symiirorder2_real(signal, r, omega, precision, status) result(y)
        real(dp), intent(in) :: signal(:) !! real samples
        real(dp), intent(in) :: r !! pole radius, abs(r) less than one
        real(dp), intent(in) :: omega !! pole angle, radians
        real(dp), optional, intent(in) :: precision !! boundary tolerance
        integer, optional, intent(out) :: status !! zero on success
        real(dp), allocatable :: y(:)
        complex(dp), allocatable :: z(:)
        integer :: code
        z = symiir_kernel(cmplx(signal,0.0_dp,dp),2,r,omega,precision,code)
        y = real(z,dp)
        if (present(status)) status = code
    end function symiirorder2_real

    function symiirorder2_complex(signal, r, omega, precision, status) result(y)
        complex(dp), intent(in) :: signal(:) !! complex samples
        real(dp), intent(in) :: r !! pole radius, abs(r) less than one
        real(dp), intent(in) :: omega !! pole angle, radians
        real(dp), optional, intent(in) :: precision !! boundary tolerance
        integer, optional, intent(out) :: status !! zero on success
        complex(dp), allocatable :: y(:)
        integer :: code
        y = symiir_kernel(signal,2,r,omega,precision,code)
        if (present(status)) status = code
    end function symiirorder2_complex

    function symiirorder1_real2d(signal, c0, z1, precision, status) result(y)
        real(dp), intent(in) :: signal(:,:) !! row-major collection of sample sequences
        real(dp), intent(in) :: c0 !! symmetric transfer-function numerator
        real(dp), intent(in) :: z1 !! stable real pole
        real(dp), optional, intent(in) :: precision !! boundary tolerance
        integer, optional, intent(out) :: status !! zero on success
        real(dp), allocatable :: y(:,:), row(:)
        integer :: i, code
        allocate(y(0,0))
        if (present(status)) status = symiir_status_invalid_input
        if (size(signal,1) < 1 .or. size(signal,2) < 1) return
        deallocate(y)
        allocate(y(size(signal,1),size(signal,2)))
        do i = 1,size(signal,1)
            row = symiirorder1_real(signal(i,:),c0,z1,precision,code)
            if (code /= 0) then
                deallocate(y)
                allocate(y(0,0))
                if (present(status)) status = code
                return
            end if
            y(i,:) = row
        end do
        if (present(status)) status = symiir_status_success
    end function symiirorder1_real2d

    function symiirorder1_complex2d(signal, c0, z1, precision, status) result(y)
        complex(dp), intent(in) :: signal(:,:) !! complex sequences along the second axis
        real(dp), intent(in) :: c0 !! symmetric transfer-function numerator
        real(dp), intent(in) :: z1 !! stable real pole
        real(dp), optional, intent(in) :: precision !! boundary tolerance
        integer, optional, intent(out) :: status !! zero on success
        complex(dp), allocatable :: y(:,:), row(:)
        integer :: i, code
        allocate(y(0,0))
        if (present(status)) status = symiir_status_invalid_input
        if (size(signal,1) < 1 .or. size(signal,2) < 1) return
        deallocate(y)
        allocate(y(size(signal,1),size(signal,2)))
        do i = 1,size(signal,1)
            row = symiirorder1_complex(signal(i,:),c0,z1,precision,code)
            if (code /= 0) then
                deallocate(y)
                allocate(y(0,0))
                if (present(status)) status = code
                return
            end if
            y(i,:) = row
        end do
        if (present(status)) status = symiir_status_success
    end function symiirorder1_complex2d

    function symiirorder2_real2d(signal, r, omega, precision, status) result(y)
        real(dp), intent(in) :: signal(:,:) !! real sequences along the second axis
        real(dp), intent(in) :: r !! pole radius
        real(dp), intent(in) :: omega !! pole angle in radians
        real(dp), optional, intent(in) :: precision !! boundary tolerance
        integer, optional, intent(out) :: status !! zero on success
        real(dp), allocatable :: y(:,:), row(:)
        integer :: i, code
        allocate(y(0,0))
        if (present(status)) status = symiir_status_invalid_input
        if (size(signal,1) < 1 .or. size(signal,2) < 1) return
        deallocate(y)
        allocate(y(size(signal,1),size(signal,2)))
        do i = 1,size(signal,1)
            row = symiirorder2_real(signal(i,:),r,omega,precision,code)
            if (code /= 0) then
                deallocate(y)
                allocate(y(0,0))
                if (present(status)) status = code
                return
            end if
            y(i,:) = row
        end do
        if (present(status)) status = symiir_status_success
    end function symiirorder2_real2d

    function symiirorder2_complex2d(signal, r, omega, precision, status) result(y)
        complex(dp), intent(in) :: signal(:,:) !! complex sequences along the second axis
        real(dp), intent(in) :: r !! pole radius
        real(dp), intent(in) :: omega !! pole angle in radians
        real(dp), optional, intent(in) :: precision !! boundary tolerance
        integer, optional, intent(out) :: status !! zero on success
        complex(dp), allocatable :: y(:,:), row(:)
        integer :: i, code
        allocate(y(0,0))
        if (present(status)) status = symiir_status_invalid_input
        if (size(signal,1) < 1 .or. size(signal,2) < 1) return
        deallocate(y)
        allocate(y(size(signal,1),size(signal,2)))
        do i = 1,size(signal,1)
            row = symiirorder2_complex(signal(i,:),r,omega,precision,code)
            if (code /= 0) then
                deallocate(y)
                allocate(y(0,0))
                if (present(status)) status = code
                return
            end if
            y(i,:) = row
        end do
        if (present(status)) status = symiir_status_success
    end function symiirorder2_complex2d
end module scifort_signal_symiir
