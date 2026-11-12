! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Analog IIR ZPK and second-order-section representations.
! Polynomial coefficients in SOS rows are descending powers of s, unlike
! the digital SOS delay-polynomial convention. See CODE_PROVENANCE.md.
module scifort_signal_analog_sos
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal_analog_design, only : analog_iir, analog_tf_result, &
        analog_design_success
    use scifort_signal_elliptic, only : elliptic_analog_zpk, elliptic_zpk_result
    use scifort_signal_analog_prototypes, only : analog_prototype_result, &
        besselap, prototype_success
    use scifort_signal_frequency_transforms, only : analog_zpk_transform_result, &
        frequency_transform_success, lp2lp_zpk, lp2hp_zpk, lp2bp_zpk, lp2bs_zpk
    use scifort_signal_design, only : frequency_response, design_status_success, &
        design_status_singular
    implicit none
    private

    integer, parameter, public :: analog_output_status_success = 0
    integer, parameter, public :: analog_output_status_invalid_input = 1

    type, public :: analog_zpk_result
        complex(dp), allocatable :: zeros(:) !! finite s-plane transmission zeros
        complex(dp), allocatable :: poles(:) !! s-plane poles
        real(dp) :: gain = 0.0_dp !! gain times product(s-z)/product(s-p)
        integer :: status = analog_output_status_invalid_input
    end type analog_zpk_result

    type, public :: analog_sos_result
        real(dp), allocatable :: sos(:,:) !! (nsections,6), descending powers of s
        integer :: status = analog_output_status_invalid_input
    end type analog_sos_result

    interface iirfilter_zpk
        module procedure iirfilter_zpk_scalar, iirfilter_zpk_vector
    end interface iirfilter_zpk
    interface iirfilter_sos
        module procedure iirfilter_sos_scalar, iirfilter_sos_vector
    end interface iirfilter_sos
    public :: iirfilter_zpk, iirfilter_sos, zpk2sos_analog, sosfreqs

contains

    function iirfilter_zpk_scalar(order, cutoff, btype, ftype, rp, rs) result(out)
        integer, intent(in) :: order !! analog prototype order, 1 through 20
        real(dp), intent(in) :: cutoff !! critical angular frequency in rad/s
        character(*), intent(in), optional :: btype !! lowpass or highpass
        character(*), intent(in), optional :: ftype !! butter, cheby1, cheby2, ellip, bessel
        real(dp), intent(in), optional :: rp !! passband ripple in dB
        real(dp), intent(in), optional :: rs !! stopband attenuation in dB
        type(analog_zpk_result) :: out
        out = iirfilter_zpk_vector(order, [cutoff], btype, ftype, rp, rs)
    end function iirfilter_zpk_scalar

    function iirfilter_zpk_vector(order, cutoff, btype, ftype, rp, rs) result(out)
        integer, intent(in) :: order !! analog prototype order, 1 through 20
        real(dp), intent(in) :: cutoff(:) !! one or two angular frequencies in rad/s
        character(*), intent(in), optional :: btype !! lowpass/highpass/bandpass/bandstop
        character(*), intent(in), optional :: ftype !! butter, cheby1, cheby2, ellip, bessel
        real(dp), intent(in), optional :: rp !! passband ripple in dB
        real(dp), intent(in), optional :: rs !! stopband attenuation in dB
        type(analog_zpk_result) :: out
        type(analog_tf_result) :: tf
        type(elliptic_zpk_result) :: e
        type(analog_prototype_result) :: proto
        type(analog_zpk_transform_result) :: mapped
        real(dp) :: center, bandwidth
        character(:), allocatable :: family, mode

        family = 'butter'
        if (present(ftype)) family = trim(ftype)
        mode = 'lowpass'
        if (present(btype)) mode = trim(btype)
        select case (family)
        case ('butter', 'butterworth')
            tf = analog_iir(order, cutoff, mode, 'butter', roots_only=.true.)
        case ('cheby1', 'chebyshev1')
            if (.not. present(rp)) return
            tf = analog_iir(order, cutoff, mode, 'cheby1', rp, roots_only=.true.)
        case ('cheby2', 'chebyshev2')
            if (.not. present(rs)) return
            tf = analog_iir(order, cutoff, mode, 'cheby2', rs, roots_only=.true.)
        case ('bessel', 'thomson')
            if (present(rp) .or. present(rs)) return
            if (size(cutoff) < 1 .or. size(cutoff) > 2) return
            if (.not. all(ieee_is_finite(cutoff))) return
            if (any(cutoff <= 0.0_dp)) return
            select case (mode)
            case ('low', 'lowpass', 'high', 'highpass')
                if (size(cutoff) /= 1) return
            case ('band', 'bandpass', 'stop', 'bandstop')
                if (size(cutoff) /= 2) return
                if (cutoff(1) >= cutoff(2)) return
            case default
                return
            end select
            proto = besselap(order)
            if (proto%status /= prototype_success) return
            select case (mode)
            case ('low', 'lowpass')
                mapped = lp2lp_zpk(proto%zeros,proto%poles,proto%gain,cutoff(1))
            case ('high', 'highpass')
                mapped = lp2hp_zpk(proto%zeros,proto%poles,proto%gain,cutoff(1))
            case ('band', 'bandpass', 'stop', 'bandstop')
                center = sqrt(cutoff(1))*sqrt(cutoff(2))
                bandwidth = cutoff(2)-cutoff(1)
                if (.not. ieee_is_finite(center) .or. .not. ieee_is_finite(bandwidth)) return
                if (center <= 0.0_dp .or. bandwidth <= 0.0_dp) return
                if (mode == 'band' .or. mode == 'bandpass') then
                    mapped = lp2bp_zpk(proto%zeros,proto%poles,proto%gain,center,bandwidth)
                else
                    mapped = lp2bs_zpk(proto%zeros,proto%poles,proto%gain,center,bandwidth)
                end if
            end select
            if (mapped%status /= frequency_transform_success) return
            if (abs(aimag(mapped%gain)) > &
                1.0e-8_dp*max(1.0_dp,abs(real(mapped%gain,dp)))) return
            out%zeros = mapped%zeros
            out%poles = mapped%poles
            out%gain = real(mapped%gain,dp)
            out%status = analog_output_status_success
            return
        case ('ellip', 'elliptic', 'cauer')
            if (.not. present(rp) .or. .not. present(rs)) return
            e = elliptic_analog_zpk(order, rp, rs, cutoff, mode)
            if (e%status /= 0) return
            out%zeros = e%zeros
            out%poles = e%poles
            out%gain = e%gain
            out%status = analog_output_status_success
            return
        case default
            return
        end select
        if (tf%status /= analog_design_success) return
        out%zeros = tf%zeros
        out%poles = tf%poles
        out%gain = tf%gain
        out%status = analog_output_status_success
    end function iirfilter_zpk_vector

    function iirfilter_sos_scalar(order, cutoff, btype, ftype, rp, rs) result(out)
        integer, intent(in) :: order !! analog prototype order
        real(dp), intent(in) :: cutoff !! critical angular frequency
        character(*), intent(in), optional :: btype !! lowpass or highpass
        character(*), intent(in), optional :: ftype !! butter, cheby1, cheby2, ellip, bessel
        real(dp), intent(in), optional :: rp !! passband ripple in dB
        real(dp), intent(in), optional :: rs !! stopband attenuation in dB
        type(analog_sos_result) :: out
        out = iirfilter_sos_vector(order, [cutoff], btype, ftype, rp, rs)
    end function iirfilter_sos_scalar

    function iirfilter_sos_vector(order, cutoff, btype, ftype, rp, rs) result(out)
        integer, intent(in) :: order !! analog prototype order
        real(dp), intent(in) :: cutoff(:) !! critical angular frequencies
        character(*), intent(in), optional :: btype !! lowpass/highpass/bandpass/bandstop
        character(*), intent(in), optional :: ftype !! butter, cheby1, cheby2, ellip, bessel
        real(dp), intent(in), optional :: rp !! passband ripple in dB
        real(dp), intent(in), optional :: rs !! stopband attenuation in dB
        type(analog_sos_result) :: out
        type(analog_zpk_result) :: zpk
        zpk = iirfilter_zpk_vector(order, cutoff, btype, ftype, rp, rs)
        if (zpk%status /= analog_output_status_success) return
        out = zpk2sos_analog(zpk%zeros, zpk%poles, zpk%gain)
    end function iirfilter_sos_vector

    function zpk2sos_analog(z, p, k) result(out)
        complex(dp), intent(in) :: z(:) !! finite transmission zeros
        complex(dp), intent(in) :: p(:) !! finite poles
        real(dp), intent(in) :: k !! transfer function gain
        type(analog_sos_result) :: out
        real(dp), allocatable :: zpoly(:,:), ppoly(:,:)
        integer, allocatable :: zdegree(:), pdegree(:)
        integer :: np, nz, ns, j, nfull
        logical :: ok

        np = size(p)
        nz = size(z)
        if (nz > np .or. np > 40) return
        if (.not. ieee_is_finite(k)) return
        if (.not. finite_roots(z) .or. .not. finite_roots(p)) return
        if (np == 0) then
            allocate(out%sos(1,6))
            out%sos(1,:) = [0.0_dp, 0.0_dp, k, 0.0_dp, 0.0_dp, 1.0_dp]
            out%status = analog_output_status_success
            return
        end if
        call factor_real_polynomial(p, ppoly, pdegree, ok)
        if (.not. ok) return
        call factor_real_polynomial(z, zpoly, zdegree, ok)
        if (.not. ok) return
        ns = (np+1)/2
        allocate(out%sos(ns,6))
        out%sos(:,1) = 0.0_dp
        out%sos(:,2) = 0.0_dp
        out%sos(:,3) = 1.0_dp
        out%sos(:,4:6) = ppoly
        ! Each quadratic numerator has its own section; any remaining
        ! first-order numerator goes in the last available section.
        nfull = nz/2
        do j = 1, nfull
            out%sos(j,1:3) = zpoly(j,:)
        end do
        if (mod(nz,2) /= 0) out%sos(ns,1:3) = zpoly(size(zpoly,1),:)
        ! Restrict each numerator's polynomial degree to its denominator's.
        do j = 1, ns
            if (out%sos(j,1) /= 0.0_dp .and. pdegree(j) < 2) return
            if (out%sos(j,2) /= 0.0_dp .and. pdegree(j) < 1) return
        end do
        out%sos(1,1:3) = k*out%sos(1,1:3)
        if (.not. all(ieee_is_finite(out%sos))) return
        out%status = analog_output_status_success
    end function zpk2sos_analog

    pure logical function finite_roots(roots) result(ok)
        complex(dp), intent(in) :: roots(:) !! finite complex roots
        ok = all(ieee_is_finite(real(roots,dp))) .and. all(ieee_is_finite(aimag(roots)))
    end function finite_roots

    subroutine factor_real_polynomial(roots, factors, degrees, ok)
        complex(dp), intent(in) :: roots(:) !! roots in conjugate pairs
        real(dp), allocatable, intent(out) :: factors(:,:) !! descending quadratic factors
        integer, allocatable, intent(out) :: degrees(:) !! degree of each factor
        logical, intent(out) :: ok !! true for conjugate-symmetric roots
        logical :: used(size(roots))
        complex(dp) :: r1, r2, c1, c2
        integer :: n, ns, keep_real, first, second, i, j
        real(dp), parameter :: root_tol = 1.0e-8_dp

        ok = .false.
        n = size(roots)
        ns = (n+1)/2
        allocate(factors(ns,3), degrees(ns))
        factors = 0.0_dp
        degrees = 2
        used = .false.
        keep_real = 0
        if (mod(n,2) /= 0) then
            do i = 1, n
                if (abs(aimag(roots(i))) <= root_tol*max(1.0_dp,abs(roots(i)))) then
                    keep_real = i
                    exit
                end if
            end do
            if (keep_real == 0) return
        end if
        do i = 1, n/2
            first = 0
            ! Prefer complex pairs, preserving a real root for odd degree.
            do j = 1, n
                if (used(j) .or. j == keep_real) cycle
                if (abs(aimag(roots(j))) > root_tol*max(1.0_dp,abs(roots(j)))) then
                    first = j
                    exit
                end if
            end do
            if (first == 0) then
                do j = 1, n
                    if (.not. used(j) .and. j /= keep_real) then
                        first = j
                        exit
                    end if
                end do
            end if
            if (first == 0) return
            r1 = roots(first)
            second = 0
            do j = 1, n
                if (used(j) .or. j == first .or. j == keep_real) cycle
                if (abs(aimag(r1)) <= root_tol*max(1.0_dp,abs(r1))) then
                    if (abs(aimag(roots(j))) > root_tol*max(1.0_dp,abs(roots(j)))) cycle
                else
                    if (abs(roots(j)-conjg(r1)) > root_tol*max(1.0_dp,abs(r1))) cycle
                end if
                second = j
                exit
            end do
            if (second == 0) return
            r2 = roots(second)
            used(first) = .true.
            used(second) = .true.
            c1 = -(r1+r2)
            c2 = r1*r2
            if (abs(aimag(c1)) > root_tol*max(1.0_dp,abs(c1))) return
            if (abs(aimag(c2)) > root_tol*max(1.0_dp,abs(c2))) return
            factors(i,:) = [1.0_dp, real(c1,dp), real(c2,dp)]
        end do
        if (keep_real /= 0) then
            factors(ns,:) = [0.0_dp, 1.0_dp, -real(roots(keep_real),dp)]
            degrees(ns) = 1
            used(keep_real) = .true.
        end if
        ok = all(used) .and. all(ieee_is_finite(factors))
    end subroutine factor_real_polynomial

    function sosfreqs(sos, w) result(out)
        real(dp), intent(in) :: sos(:,:) !! analog sections in descending powers of s
        real(dp), intent(in) :: w(:) !! angular frequencies in rad/s
        type(frequency_response) :: out
        complex(dp) :: s, h, num, den
        integer :: i, j

        if (size(sos,2) /= 6 .or. size(sos,1) < 1) return
        if (.not. all(ieee_is_finite(sos))) return
        if (.not. all(ieee_is_finite(w))) return
        if (any(all(sos(:,4:6) == 0.0_dp,dim=2))) return
        out%frequencies = w
        allocate(out%response(size(w)))
        out%status = design_status_success
        do i = 1, size(w)
            s = cmplx(0.0_dp,w(i),dp)
            h = cmplx(1.0_dp,0.0_dp,dp)
            do j = 1, size(sos,1)
                num = (sos(j,1)*s+sos(j,2))*s+sos(j,3)
                den = (sos(j,4)*s+sos(j,5))*s+sos(j,6)
                if (abs(den) <= tiny(1.0_dp)) then
                    out%status = design_status_singular
                    h = cmplx(ieee_value(0.0_dp,ieee_quiet_nan), &
                        ieee_value(0.0_dp,ieee_quiet_nan), dp)
                    exit
                end if
                h = h*(num/den)
            end do
            out%response(i) = h
        end do
    end function sosfreqs

end module scifort_signal_analog_sos
