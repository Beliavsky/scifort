! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Bilinear s-to-z transformations in transfer-function and ZPK forms.
! Frequencies are in hertz; no frequency prewarping is performed.
module scifort_signal_bilinear
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_signal_frequency_transforms, only : analog_tf_transform_result, &
        analog_zpk_transform_result, frequency_transform_success, &
        frequency_transform_invalid_input, frequency_transform_too_expensive
    implicit none
    private

    integer, parameter :: max_bilinear_order = 32
    integer, parameter :: max_bilinear_roots = 128
    public :: bilinear, bilinear_zpk

    interface bilinear
        module procedure bilinear_real, bilinear_complex
    end interface bilinear
    interface bilinear_zpk
        module procedure bilinear_zpk_real_gain, bilinear_zpk_complex_gain
    end interface bilinear_zpk
contains
    pure logical function finite_complex(z) result(ok)
        complex(dp), intent(in) :: z(:) !! complex array to validate
        ok = all(ieee_is_finite(real(z, dp))) .and. all(ieee_is_finite(aimag(z)))
    end function finite_complex

    pure function linear_product(poly, c0, c1) result(out)
        complex(dp), intent(in) :: poly(:) !! ascending polynomial coefficients
        real(dp), intent(in) :: c0 !! constant factor of the linear polynomial
        real(dp), intent(in) :: c1 !! linear factor of the linear polynomial
        complex(dp) :: out(size(poly)+1)
        integer :: j
        out = (0.0_dp, 0.0_dp)
        do j = 1, size(poly)
            out(j) = out(j) + c0*poly(j)
            out(j+1) = out(j+1) + c1*poly(j)
        end do
    end function linear_product

    function substituted_poly(coeff, degree, q) result(out)
        complex(dp), intent(in) :: coeff(:) !! descending analog coefficients
        integer, intent(in) :: degree !! max analog polynomial degree
        real(dp), intent(in) :: q !! twice the sampling rate
        complex(dp) :: out(degree+1)
        complex(dp), allocatable :: term(:)
        integer :: j, k, power

        out = (0.0_dp, 0.0_dp)
        do j = 1, size(coeff)
            if (abs(coeff(j)) == 0.0_dp) cycle
            power = size(coeff)-j
            term = [cmplx(1.0_dp,0.0_dp,dp)]
            do k = 1, power
                term = linear_product(term,-q,q)
            end do
            do k = power+1, degree
                term = linear_product(term,1.0_dp,1.0_dp)
            end do
            out = out + coeff(j)*term
        end do
    end function substituted_poly

    function bilinear_complex(b, a, fs) result(out)
        complex(dp), intent(in) :: b(:) !! analog numerator, descending powers of s
        complex(dp), intent(in) :: a(:) !! analog denominator, descending powers of s
        real(dp), intent(in), optional :: fs !! sampling rate in hertz (default 1)
        type(analog_tf_transform_result) :: out
        complex(dp), allocatable :: pb(:), pa(:)
        real(dp) :: rate
        integer :: degree

        rate = 1.0_dp
        if (present(fs)) rate = fs
        if (.not. ieee_is_finite(rate)) return
        if (rate <= 0.0_dp .or. rate >= huge(rate)/2.0_dp) return
        if (size(a) == 0 .or. size(b) == 0) return
        if (.not. finite_complex(a) .or. .not. finite_complex(b)) return
        if (all(abs(a) == 0.0_dp)) return
        degree = max(size(a),size(b))-1
        if (degree > max_bilinear_order) then
            out%status = frequency_transform_too_expensive
            return
        end if
        pb = substituted_poly(b,degree,2.0_dp*rate)
        pa = substituted_poly(a,degree,2.0_dp*rate)
        if (.not. finite_complex(pa) .or. .not. finite_complex(pb)) return
        if (abs(pa(degree+1)) <= tiny(rate)) return
        out%b = pb(degree+1:1:-1)/pa(degree+1)
        out%a = pa(degree+1:1:-1)/pa(degree+1)
        if (.not. finite_complex(out%b) .or. .not. finite_complex(out%a)) then
            deallocate(out%b,out%a)
            return
        end if
        out%status = frequency_transform_success
    end function bilinear_complex

    function bilinear_real(b, a, fs) result(out)
        real(dp), intent(in) :: b(:) !! real analog numerator, descending powers
        real(dp), intent(in) :: a(:) !! real analog denominator, descending powers
        real(dp), intent(in), optional :: fs !! sampling rate in hertz (default 1)
        type(analog_tf_transform_result) :: out
        out = bilinear_complex(cmplx(b,0.0_dp,dp),cmplx(a,0.0_dp,dp),fs)
    end function bilinear_real

    function bilinear_zpk_complex_gain(z, p, k, fs) result(out)
        complex(dp), intent(in) :: z(:) !! analog zeros
        complex(dp), intent(in) :: p(:) !! analog poles
        complex(dp), intent(in) :: k !! analog gain
        real(dp), intent(in) :: fs !! sampling rate in hertz
        type(analog_zpk_transform_result) :: out
        real(dp) :: fs2
        integer :: degree
        complex(dp) :: gain_ratio

        if (.not. ieee_is_finite(fs)) return
        if (fs <= 0.0_dp .or. fs >= huge(fs)/2.0_dp) return
        if (size(z) > size(p)) return
        if (size(p) > max_bilinear_roots) then
            out%status = frequency_transform_too_expensive
            return
        end if
        if (.not. finite_complex(z) .or. .not. finite_complex(p)) return
        if (.not. ieee_is_finite(real(k,dp)) .or. .not. ieee_is_finite(aimag(k))) return
        fs2 = 2.0_dp*fs
        if (any(abs(fs2-z) == 0.0_dp)) return
        if (any(abs(fs2-p) == 0.0_dp)) return
        degree = size(p)-size(z)
        allocate(out%zeros(size(p)),out%poles(size(p)))
        if (size(z) > 0) out%zeros(:size(z)) = (fs2+z)/(fs2-z)
        if (degree > 0) out%zeros(size(z)+1:) = cmplx(-1.0_dp,0.0_dp,dp)
        if (size(p) > 0) out%poles = (fs2+p)/(fs2-p)
        gain_ratio = cmplx(1.0_dp,0.0_dp,dp)
        if (size(z) > 0) gain_ratio = product(fs2-z)
        if (size(p) > 0) gain_ratio = gain_ratio/product(fs2-p)
        out%gain = k*gain_ratio
        if (.not. finite_complex(out%zeros) .or. .not. finite_complex(out%poles)) return
        if (.not. ieee_is_finite(real(out%gain,dp)) .or. &
            .not. ieee_is_finite(aimag(out%gain))) return
        out%status = frequency_transform_success
    end function bilinear_zpk_complex_gain

    function bilinear_zpk_real_gain(z, p, k, fs) result(out)
        complex(dp), intent(in) :: z(:) !! analog zeros
        complex(dp), intent(in) :: p(:) !! analog poles
        real(dp), intent(in) :: k !! real analog gain
        real(dp), intent(in) :: fs !! sampling rate in hertz
        type(analog_zpk_transform_result) :: out
        out = bilinear_zpk_complex_gain(z,p,cmplx(k,0.0_dp,dp),fs)
    end function bilinear_zpk_real_gain
end module scifort_signal_bilinear
