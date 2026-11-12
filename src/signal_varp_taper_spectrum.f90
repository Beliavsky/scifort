! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Exact finite-window population cross-spectra for stable real VAR(p).
! Band integrals are analytic, not quadrature over a sampled frequency grid.
! Original SciFort extension; these functions are not scipy.signal APIs.
module scifort_signal_varp_taper_spectrum
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    implicit none
    private

    integer, parameter, public :: varp_spectrum_success = 0
    integer, parameter, public :: varp_spectrum_invalid_input = 1
    integer, parameter, public :: varp_spectrum_nonstationary = 2
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: twopi = 2.0_dp*pi

    type, public :: varp_taper_spectrum_result
        real(dp), allocatable :: frequencies(:) !! cycles/sample; band centers for band results
        real(dp), allocatable :: bands(:,:) !! exact integration boundaries (2,nband), if bands
        real(dp), allocatable :: lag0_covariance(:,:) !! stationary covariance of the observed channels
        complex(dp), allocatable :: power(:,:,:) !! expected spectral matrix (frequency or band, channel, channel)
        integer :: status = varp_spectrum_invalid_input !! 0 success, 1 invalid, 2 unstable/no convergence
    end type varp_taper_spectrum_result

    public :: varp_expected_taper_spectrum, varp_expected_band_taper_spectrum

contains

    pure logical function positive_definite(a) result(ok)
        real(dp), intent(in) :: a(:,:) !! symmetric candidate innovation or state covariance
        real(dp) :: lower(size(a,1),size(a,1)), pivot, scale
        integer :: i, j, k, n

        n = size(a,1)
        ok = .false.
        if (n < 1 .or. size(a,2) /= n) return
        scale = maxval(abs(a))
        if (.not. ieee_is_finite(scale) .or. scale <= 0.0_dp) return
        lower = 0.0_dp
        do i = 1, n
            pivot = a(i,i)-sum(lower(i,1:i-1)**2)
            if (.not. ieee_is_finite(pivot)) return
            if (pivot <= 64.0_dp*epsilon(1.0_dp)*scale) return
            lower(i,i) = sqrt(pivot)
            do j = i+1, n
                lower(j,i) = a(j,i)
                do k = 1, i-1
                    lower(j,i) = lower(j,i)-lower(j,k)*lower(i,k)
                end do
                lower(j,i) = lower(j,i)/lower(i,i)
            end do
        end do
        ok = .true.
    end function positive_definite

    ! Smith doubling works for the singular innovation covariance of the
    ! companion process; only the observed innovation block must be PD.
    subroutine state_covariance(f,q,gamma,ok)
        real(dp), intent(in) :: f(:,:) !! companion transition matrix (np,np)
        real(dp), intent(in) :: q(:,:) !! state innovation covariance, positive semidefinite
        real(dp), intent(out) :: gamma(:,:) !! stationary state covariance
        logical, intent(out) :: ok !! finite and converged Lyapunov series
        real(dp) :: transition(size(f,1),size(f,2)), delta(size(f,1),size(f,2))
        integer :: iteration

        ok = .false.
        gamma = q
        transition = f
        do iteration = 1, 80
            if (.not. all(ieee_is_finite(transition))) return
            delta = matmul(transition,matmul(gamma,transpose(transition)))
            if (.not. all(ieee_is_finite(delta))) return
            if (maxval(abs(delta)) <= 32.0_dp*epsilon(1.0_dp)*maxval(abs(gamma))) then
                gamma = gamma+delta
                gamma = 0.5_dp*(gamma+transpose(gamma))
                ok = all(ieee_is_finite(gamma))
                return
            end if
            gamma = gamma+delta
            if (.not. all(ieee_is_finite(gamma))) return
            transition = matmul(transition,transition)
        end do
    end subroutine state_covariance

    subroutine window_covariances(a,sigma,taper,covariance,weights,status)
        real(dp), intent(in) :: a(:,:,:) !! VAR coefficients (channel,channel,lag), lag 1..p
        real(dp), intent(in) :: sigma(:,:) !! symmetric positive-definite innovation covariance
        real(dp), intent(in) :: taper(:) !! finite window, nonzero sum of squares
        real(dp), allocatable, intent(out) :: covariance(:,:,:) !! C(h)=E[x(t+h)x(t)^T], h=0..length-1
        real(dp), allocatable, intent(out) :: weights(:) !! taper lag autocorrelations, h=0..length-1
        integer, intent(out) :: status !! 0 success; 1 invalid input; 2 nonstationary
        real(dp), allocatable :: f(:,:), q(:,:), gamma(:,:), state_cross(:,:)
        integer :: n, p, l, m, i, j, h
        logical :: ok

        status = varp_spectrum_invalid_input
        n = size(a,1)
        p = size(a,3)
        l = size(taper)
        if (n < 1 .or. p < 1 .or. l < 1) return
        if (size(a,2) /= n) return
        if (size(sigma,1) /= n .or. size(sigma,2) /= n) return
        if (.not. all(ieee_is_finite(a))) return
        if (.not. all(ieee_is_finite(sigma))) return
        if (.not. all(ieee_is_finite(taper))) return
        if (maxval(abs(sigma-transpose(sigma))) > &
            64.0_dp*epsilon(1.0_dp)*maxval(abs(sigma))) return
        if (.not. positive_definite(0.5_dp*(sigma+transpose(sigma)))) return
        if (p > huge(p)/n) return
        m = p*n
        allocate(weights(0:l-1))
        weights(0) = sum(taper**2)
        if (.not. ieee_is_finite(weights(0)) .or. weights(0) <= 0.0_dp) return
        do h = 1, l-1
            weights(h) = sum(taper(1:l-h)*taper(h+1:l))
        end do
        allocate(f(m,m),q(m,m),gamma(m,m),state_cross(m,n))
        f = 0.0_dp
        q = 0.0_dp
        q(1:n,1:n) = sigma
        do j = 1, p
            f(1:n,(j-1)*n+1:j*n) = a(:,:,j)
        end do
        do j = 2, p
            do i = 1, n
                f((j-1)*n+i,(j-2)*n+i) = 1.0_dp
            end do
        end do
        call state_covariance(f,q,gamma,ok)
        if (.not. ok) then
            status = varp_spectrum_nonstationary
            return
        end if
        if (.not. positive_definite(gamma(1:n,1:n))) then
            status = varp_spectrum_nonstationary
            return
        end if
        allocate(covariance(n,n,0:l-1))
        state_cross = gamma(:,1:n)
        covariance(:,:,0) = gamma(1:n,1:n)
        do h = 1, l-1
            state_cross = matmul(f,state_cross)
            if (.not. all(ieee_is_finite(state_cross))) then
                deallocate(covariance)
                status = varp_spectrum_nonstationary
                return
            end if
            covariance(:,:,h) = state_cross(1:n,:)
        end do
        status = varp_spectrum_success
    end subroutine window_covariances

    function varp_expected_taper_spectrum(a,innovation_covariance,taper,frequencies, &
            normalized) result(out)
        real(dp), intent(in) :: a(:,:,:) !! stable real VAR(p) coefficients (n,n,p)
        real(dp), intent(in) :: innovation_covariance(:,:) !! positive-definite real innovation covariance
        real(dp), intent(in) :: taper(:) !! sample-by-sample window weights
        real(dp), intent(in) :: frequencies(:) !! cycles/sample; may be negative, wrapped, or unsorted
        logical, intent(in), optional :: normalized !! divide by sum(taper**2); default false
        type(varp_taper_spectrum_result) :: out
        real(dp), allocatable :: cov(:,:,:), weight(:)
        real(dp) :: angle, normalizer
        complex(dp) :: phase
        integer :: n, nf, f, h, status

        nf = size(frequencies)
        if (nf < 1) return
        if (.not. all(ieee_is_finite(frequencies))) return
        call window_covariances(a,innovation_covariance,taper,cov,weight,status)
        out%status = status
        if (status /= varp_spectrum_success) return
        n = size(a,1)
        allocate(out%frequencies(nf),out%lag0_covariance(n,n),out%power(nf,n,n))
        out%frequencies = frequencies
        out%lag0_covariance = cov(:,:,0)
        do f = 1, nf
            out%power(f,:,:) = cmplx(weight(0)*cov(:,:,0),0.0_dp,dp)
            do h = 1, size(taper)-1
                angle = twopi*modulo(modulo(frequencies(f),1.0_dp)*real(h,dp),1.0_dp)
                phase = cmplx(cos(angle),-sin(angle),dp)
                out%power(f,:,:) = out%power(f,:,:)+weight(h)* &
                    (phase*cov(:,:,h)+conjg(phase)*transpose(cov(:,:,h)))
            end do
        end do
        normalizer = 1.0_dp
        if (present(normalized)) then
            if (normalized) normalizer = weight(0)
        end if
        out%power = out%power/normalizer
        if (.not. all(ieee_is_finite(real(out%power,dp))) .or. &
            .not. all(ieee_is_finite(aimag(out%power)))) then
            deallocate(out%frequencies,out%lag0_covariance,out%power)
            out%status = varp_spectrum_invalid_input
        end if
    end function varp_expected_taper_spectrum

    ! Analytic integral of exp(-2*pi*i*f*h) on [left,right], avoiding
    ! cancellation when the frequency band is very narrow.
    pure function lag_band_integral(h,left,right) result(value)
        integer, intent(in) :: h !! positive lag in samples
        real(dp), intent(in) :: left !! lower integration frequency, cycles/sample
        real(dp), intent(in) :: right !! upper integration frequency, cycles/sample
        complex(dp) :: value
        real(dp) :: width, x, sinc_value, mid, angle

        width = right-left
        x = pi*real(h,dp)*width
        if (abs(x) < 1.0e-3_dp) then
            sinc_value = 1.0_dp-x*x/6.0_dp+x**4/120.0_dp-x**6/5040.0_dp
        else
            sinc_value = sin(x)/x
        end if
        mid = 0.5_dp*left+0.5_dp*right
        angle = twopi*modulo(modulo(mid,1.0_dp)*real(h,dp),1.0_dp)
        value = width*sinc_value*cmplx(cos(angle),-sin(angle),dp)
    end function lag_band_integral

    function varp_expected_band_taper_spectrum(a,innovation_covariance,taper,bands, &
            average,normalized) result(out)
        real(dp), intent(in) :: a(:,:,:) !! stable real VAR(p) coefficients (n,n,p)
        real(dp), intent(in) :: innovation_covariance(:,:) !! positive-definite innovation covariance
        real(dp), intent(in) :: taper(:) !! sample-by-sample window weights
        real(dp), intent(in) :: bands(:,:) !! lower/upper edges in cycles/sample, (2,number of bands)
        logical, intent(in), optional :: average !! return band-mean spectral density; default false
        logical, intent(in), optional :: normalized !! divide by sum(taper**2); default false
        type(varp_taper_spectrum_result) :: out
        real(dp), allocatable :: cov(:,:,:), weight(:)
        real(dp) :: width, normalizer
        complex(dp) :: factor
        integer :: n, nb, b, h, status

        if (size(bands,1) /= 2) return
        nb = size(bands,2)
        if (nb < 1) return
        if (.not. all(ieee_is_finite(bands))) return
        if (any(bands(2,:) <= bands(1,:))) return
        if (.not. all(ieee_is_finite(bands(2,:)-bands(1,:)))) return
        call window_covariances(a,innovation_covariance,taper,cov,weight,status)
        out%status = status
        if (status /= varp_spectrum_success) return
        n = size(a,1)
        allocate(out%frequencies(nb),out%bands(2,nb),out%lag0_covariance(n,n),out%power(nb,n,n))
        out%bands = bands
        out%frequencies = 0.5_dp*bands(1,:)+0.5_dp*bands(2,:)
        out%lag0_covariance = cov(:,:,0)
        do b = 1, nb
            width = bands(2,b)-bands(1,b)
            out%power(b,:,:) = cmplx(width*weight(0)*cov(:,:,0),0.0_dp,dp)
            do h = 1, size(taper)-1
                factor = lag_band_integral(h,bands(1,b),bands(2,b))
                out%power(b,:,:) = out%power(b,:,:)+weight(h)* &
                    (factor*cov(:,:,h)+conjg(factor)*transpose(cov(:,:,h)))
            end do
            normalizer = 1.0_dp
            if (present(normalized)) then
                if (normalized) normalizer = weight(0)
            end if
            if (present(average)) then
                if (average) normalizer = normalizer*width
            end if
            out%power(b,:,:) = out%power(b,:,:)/normalizer
        end do
        if (.not. all(ieee_is_finite(real(out%power,dp))) .or. &
            .not. all(ieee_is_finite(aimag(out%power)))) then
            deallocate(out%frequencies,out%bands,out%lag0_covariance,out%power)
            out%status = varp_spectrum_invalid_input
        end if
    end function varp_expected_band_taper_spectrum

end module scifort_signal_varp_taper_spectrum
