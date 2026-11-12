! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Exact finite-window population cross-spectral matrix for a stationary
! real VAR(1). This SciFort extension is not part of scipy.signal.
module scifort_signal_var1_taper_spectrum
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    implicit none
    private

    integer, parameter, public :: var1_spectrum_success = 0
    integer, parameter, public :: var1_spectrum_invalid_input = 1
    integer, parameter, public :: var1_spectrum_nonstationary = 2
    real(dp), parameter :: twopi = 2.0_dp*acos(-1.0_dp)

    type, public :: var1_taper_spectrum_result
        real(dp), allocatable :: frequencies(:) !! cycles per sample, supplied by caller
        real(dp), allocatable :: lag0_covariance(:,:) !! stationary VAR(1) covariance
        complex(dp), allocatable :: power(:,:,:) !! expected tapered outer product (frequency,channel,channel)
        integer :: status = var1_spectrum_invalid_input !! 0 success, 1 invalid, 2 unstable or no convergence
    end type var1_taper_spectrum_result

    public :: var1_expected_taper_spectrum

contains

    ! Positive definiteness check without external LAPACK dependencies.
    pure logical function positive_definite(s) result(ok)
        real(dp), intent(in) :: s(:,:) !! symmetric innovation covariance
        real(dp) :: lower(size(s,1),size(s,1)), pivot, scale
        integer :: i, j, k, n

        n = size(s,1)
        lower = 0.0_dp
        scale = maxval(abs(s))
        ok = .false.
        if (scale <= 0.0_dp) return
        do i = 1, n
            pivot = s(i,i)-sum(lower(i,1:i-1)**2)
            if (pivot <= 64.0_dp*epsilon(scale)*scale) return
            lower(i,i) = sqrt(pivot)
            do j = i+1, n
                lower(j,i) = s(j,i)
                do k = 1, i-1
                    lower(j,i) = lower(j,i)-lower(j,k)*lower(i,k)
                end do
                lower(j,i) = lower(j,i)/lower(i,i)
            end do
        end do
        ok = .true.
    end function positive_definite

    ! Smith doubling: Gamma = Sigma + A Gamma A^T, B <- B^2.
    ! The increment shrinks quadratically for a stable state matrix.
    subroutine stationary_covariance(a,sigma,gamma,success)
        real(dp), intent(in) :: a(:,:) !! real VAR(1) transition, spectral radius below one
        real(dp), intent(in) :: sigma(:,:) !! positive-definite innovation covariance
        real(dp), intent(out) :: gamma(:,:) !! stationary lag-zero covariance if successful
        logical, intent(out) :: success !! convergence and finite arithmetic
        real(dp) :: b(size(a,1),size(a,2)), delta(size(a,1),size(a,2))
        integer :: iteration

        gamma = sigma
        b = a
        success = .false.
        do iteration = 1, 80
            if (.not. all(ieee_is_finite(b))) return
            delta = matmul(b,matmul(gamma,transpose(b)))
            if (.not. all(ieee_is_finite(delta))) return
            if (maxval(abs(delta)) <= 32.0_dp*epsilon(1.0_dp)*maxval(abs(gamma))) then
                gamma = gamma+delta
                gamma = 0.5_dp*(gamma+transpose(gamma))
                success = all(ieee_is_finite(gamma))
                return
            end if
            gamma = gamma+delta
            if (.not. all(ieee_is_finite(gamma))) return
            b = matmul(b,b)
        end do
    end subroutine stationary_covariance

    function var1_expected_taper_spectrum(a,innovation_covariance,taper,frequencies, &
        normalized) result(result)
        real(dp), intent(in) :: a(:,:) !! stable VAR(1) coefficient matrix (channels,channels)
        real(dp), intent(in) :: innovation_covariance(:,:) !! symmetric positive-definite innovation covariance
        real(dp), intent(in) :: taper(:) !! finite window weights, nonzero squared norm
        real(dp), intent(in) :: frequencies(:) !! finite Fourier frequencies in cycles/sample
        logical, intent(in), optional :: normalized !! divide by sum(taper**2), default false
        type(var1_taper_spectrum_result) :: result
        real(dp), allocatable :: sigma(:,:), gamma(:,:), lag_cov(:,:), power_a(:,:)
        real(dp) :: weight, norm_weight, angle
        complex(dp) :: phase
        integer :: n, l, nf, h, t, f, i, j
        logical :: converged

        n = size(a,1)
        l = size(taper)
        nf = size(frequencies)
        if (n < 1 .or. l < 1 .or. nf < 1) return
        if (size(a,2) /= n) return
        if (size(innovation_covariance,1) /= n .or. &
            size(innovation_covariance,2) /= n) return
        if (.not. all(ieee_is_finite(a))) return
        if (.not. all(ieee_is_finite(innovation_covariance))) return
        if (.not. all(ieee_is_finite(taper))) return
        if (.not. all(ieee_is_finite(frequencies))) return
        norm_weight = sum(taper**2)
        if (.not. ieee_is_finite(norm_weight) .or. norm_weight <= 0.0_dp) return
        if (maxval(abs(innovation_covariance-transpose(innovation_covariance))) > &
            64.0_dp*epsilon(1.0_dp)*maxval(abs(innovation_covariance))) return
        sigma = 0.5_dp*(innovation_covariance+transpose(innovation_covariance))
        if (.not. positive_definite(sigma)) return
        allocate(gamma(n,n),lag_cov(n,n),power_a(n,n))
        call stationary_covariance(a,sigma,gamma,converged)
        if (.not. converged) then
            result%status = var1_spectrum_nonstationary
            return
        end if
        if (.not. positive_definite(gamma)) then
            result%status = var1_spectrum_nonstationary
            return
        end if
        allocate(result%frequencies(nf),result%lag0_covariance(n,n),result%power(nf,n,n))
        result%frequencies = frequencies
        result%lag0_covariance = gamma
        do f = 1, nf
            result%power(f,:,:) = cmplx(norm_weight*gamma,0.0_dp,dp)
        end do
        power_a = a
        do h = 1, l-1
            weight = 0.0_dp
            do t = 1, l-h
                weight = weight+taper(t+h)*taper(t)
            end do
            lag_cov = matmul(power_a,gamma)
            do f = 1, nf
                angle = twopi*modulo(modulo(frequencies(f),1.0_dp)*real(h,dp),1.0_dp)
                phase = cmplx(cos(angle),-sin(angle),dp)
                do j = 1, n
                    do i = 1, n
                        result%power(f,i,j) = result%power(f,i,j)+weight* &
                            (phase*lag_cov(i,j)+conjg(phase)*lag_cov(j,i))
                    end do
                end do
            end do
            power_a = matmul(power_a,a)
        end do
        if (present(normalized)) then
            if (normalized) result%power = result%power/norm_weight
        end if
        if (.not. all(ieee_is_finite(real(result%power,dp))) .or. &
            .not. all(ieee_is_finite(aimag(result%power)))) then
            deallocate(result%frequencies,result%lag0_covariance,result%power)
            return
        end if
        result%status = var1_spectrum_success
    end function var1_expected_taper_spectrum

end module scifort_signal_var1_taper_spectrum
