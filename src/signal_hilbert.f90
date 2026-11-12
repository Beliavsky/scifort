! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! FFT-based one- and two-dimensional analytic signals (SciPy signal.hilbert,
! signal.hilbert2). The two-dimensional transform uses the single-orthant
! convention, not a transform applied independently to a real part.
module scifort_signal_hilbert
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_signal_spectral, only : forward_fft
    implicit none
    private

    integer, parameter, public :: hilbert_status_success = 0
    integer, parameter, public :: hilbert_status_invalid_input = 1

    interface hilbert2
        module procedure hilbert2_shape
        module procedure hilbert2_square
    end interface hilbert2

    public :: hilbert, hilbert2

contains

    ! Frequencies are in FFT order, starting with the DC component.
    pure function positive_frequency_mask(n) result(mask)
        integer, intent(in) :: n !! positive transform length
        real(dp) :: mask(n)

        mask = 0.0_dp
        mask(1) = 1.0_dp
        if (mod(n, 2) == 0) then
            mask(2:n/2) = 2.0_dp
            mask(n/2+1) = 1.0_dp
        else
            mask(2:(n+1)/2) = 2.0_dp
        end if
    end function positive_frequency_mask

    function hilbert(x, n, status) result(y)
        real(dp), intent(in) :: x(:) !! real-valued input signal
        integer, intent(in), optional :: n !! FFT length, truncates or zero-pads input
        integer, intent(out), optional :: status !! success or invalid input
        complex(dp), allocatable :: y(:)
        complex(dp), allocatable :: work(:)
        integer :: nf, ncopy

        if (present(status)) status = hilbert_status_invalid_input
        allocate(y(0))
        nf = size(x)
        if (present(n)) nf = n
        if (nf < 1) return
        if (.not. all(ieee_is_finite(x))) return

        allocate(work(nf))
        work = cmplx(0.0_dp, 0.0_dp, dp)
        ncopy = min(size(x), nf)
        work(:ncopy) = cmplx(x(:ncopy), 0.0_dp, dp)
        call forward_fft(work)
        work = work * positive_frequency_mask(nf)
        work = conjg(work)
        call forward_fft(work)
        deallocate(y)
        y = conjg(work)/real(nf,dp)
        if (present(status)) status = hilbert_status_success
    end function hilbert

    function hilbert2_shape(x, n, status) result(y)
        real(dp), intent(in) :: x(:,:) !! two-dimensional real-valued samples
        integer, intent(in), optional :: n(:) !! two transform lengths, optional
        integer, intent(out), optional :: status !! success or invalid input
        complex(dp), allocatable :: y(:,:)
        complex(dp), allocatable :: work(:,:)
        real(dp), allocatable :: mask1(:), mask2(:)
        integer :: nr, nc, i, j

        if (present(status)) status = hilbert_status_invalid_input
        allocate(y(0,0))
        nr = size(x,1)
        nc = size(x,2)
        if (present(n)) then
            if (size(n) /= 2) return
            nr = n(1)
            nc = n(2)
        end if
        if (nr < 1 .or. nc < 1) return
        if (.not. all(ieee_is_finite(x))) return

        allocate(work(nr,nc))
        work = cmplx(0.0_dp, 0.0_dp, dp)
        work(:min(size(x,1),nr),:min(size(x,2),nc)) = &
            cmplx(x(:min(size(x,1),nr),:min(size(x,2),nc)),0.0_dp,dp)

        do j = 1, nc
            call forward_fft(work(:,j))
        end do
        do i = 1, nr
            call forward_fft(work(i,:))
        end do

        mask1 = positive_frequency_mask(nr)
        mask2 = positive_frequency_mask(nc)
        ! SciPy hilbert2 (1.17) excludes the unpaired Nyquist bins when
        ! either transformed axis has even length, unlike one-dimensional hilbert.
        if (mod(nr, 2) == 0) mask1(nr/2+1) = 0.0_dp
        if (mod(nc, 2) == 0) mask2(nc/2+1) = 0.0_dp
        do j = 1, nc
            do i = 1, nr
                work(i,j) = work(i,j) * mask1(i) * mask2(j)
            end do
        end do

        work = conjg(work)
        do j = 1, nc
            call forward_fft(work(:,j))
        end do
        do i = 1, nr
            call forward_fft(work(i,:))
        end do
        deallocate(y)
        y = conjg(work)/(real(nr,dp)*real(nc,dp))
        if (present(status)) status = hilbert_status_success
    end function hilbert2_shape

    function hilbert2_square(x, n, status) result(y)
        real(dp), intent(in) :: x(:,:) !! real-valued two-dimensional samples
        integer, intent(in) :: n !! common FFT length for both axes
        integer, intent(out), optional :: status !! success or invalid input
        complex(dp), allocatable :: y(:,:)

        y = hilbert2_shape(x, [n,n], status)
    end function hilbert2_square

end module scifort_signal_hilbert
