! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Chebyshev-II inverse-ripple, analog/digital IIR selection, and notch/peak design.
! Analog prototype equations, bilinear mapping and biquad formulas follow
! standard signal-processing texts; see CODE_PROVENANCE.md.
module scifort_signal_iir_extra
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_signal_analog_design, only : analog_iir, analog_tf_result, analog_design_success
    use scifort_signal_elliptic, only : elliptic_analog_zpk, elliptic_zpk_result, &
        elliptic_k, elliptic_k_complement
    use scifort_signal_bessel_design, only : bessel
    use scifort_signal_design, only : butter_result, butter, cheby1, &
        design_status_success, design_status_invalid_input
    implicit none
    private
    real(dp), parameter :: pi = acos(-1.0_dp)

    type, public :: iir_order_result
        integer :: order = 0
        real(dp) :: cutoff = 0.0_dp
        integer :: status = design_status_invalid_input
        real(dp), allocatable :: cutoffs(:) !! one or two critical frequencies
    end type iir_order_result

    interface cheby2
        module procedure cheby2_scalar
        module procedure cheby2_vector
    end interface
    interface ellip
        module procedure ellip_scalar, ellip_vector
    end interface ellip
    interface iirfilter
        module procedure iirfilter_scalar
        module procedure iirfilter_vector
    end interface
    interface buttord
        module procedure buttord_scalar, buttord_vector
    end interface buttord
    interface cheb1ord
        module procedure cheb1ord_scalar, cheb1ord_vector
    end interface cheb1ord
    interface cheb2ord
        module procedure cheb2ord_scalar, cheb2ord_vector
    end interface cheb2ord
    interface ellipord
        module procedure ellipord_scalar, ellipord_vector
    end interface ellipord
    interface iirdesign
        module procedure iirdesign_scalar, iirdesign_vector
    end interface iirdesign
    public :: cheby2, ellip, iirfilter, iirdesign, buttord, cheb1ord, cheb2ord, &
        ellipord, iirnotch, iirpeak, iircomb, gammatone

contains

    function cheby2_scalar(order, rs, cutoff, btype, fs, analog) result(out)
        integer, intent(in) :: order !! lowpass prototype order, 1 through 20
        real(dp), intent(in) :: rs !! minimum stopband attenuation in dB
        real(dp), intent(in) :: cutoff !! digital critical frequency
        character(*), intent(in), optional :: btype !! lowpass (default), highpass
        real(dp), intent(in), optional :: fs !! sampling rate, default 2
        logical, intent(in), optional :: analog !! analog s-domain design, default false
        type(butter_result) :: out
        out = cheby2_vector(order, rs, [cutoff], btype, fs, analog)
    end function cheby2_scalar

    function cheby2_vector(order, rs, cutoff, btype, fs, analog) result(out)
        integer, intent(in) :: order !! lowpass prototype order, 1 through 20
        real(dp), intent(in) :: rs !! minimum stopband attenuation in dB
        real(dp), intent(in) :: cutoff(:) !! one or two digital critical frequencies
        character(*), intent(in), optional :: btype !! lowpass, highpass, bandpass, bandstop
        real(dp), intent(in), optional :: fs !! sampling rate, default 2
        logical, intent(in), optional :: analog !! analog s-domain design, default false
        type(butter_result) :: out
        type(analog_tf_result) :: tf
        complex(dp), allocatable :: p0(:), z0(:), poles(:), zeros(:)
        complex(dp) :: root, disc, v, scale
        logical, allocatable :: at_infinity(:)
        real(dp) :: rate, omega1, omega2, bw, center, eps, mu, theta, m
        integer :: i, nz, n, np, j
        character(:), allocatable :: mode
        if (present(analog)) then
            if (analog) then
                if (present(fs)) return
                tf = analog_iir(order, cutoff, btype, 'cheby2', rs)
                if (tf%status /= analog_design_success) return
                out%b = tf%b
                out%a = tf%a
                out%status = design_status_success
                return
            end if
        end if
        ! Validate before allocating arrays and avoid unbounded polynomial orders.
        rate = 2.0_dp
        if (present(fs)) rate = fs
        if (order < 1 .or. order > 20) return
        if (.not. ieee_is_finite(rate) .or. .not. ieee_is_finite(rs)) return
        if (rate <= 0.0_dp .or. rs <= 0.0_dp .or. rs > 100.0_dp) return
        mode = 'lowpass'
        if (present(btype)) mode = trim(btype)
        select case (mode)
        case ('lowpass', 'highpass')
            if (size(cutoff) /= 1) return
        case ('bandpass', 'bandstop')
            if (size(cutoff) /= 2) return
            if (cutoff(1) >= cutoff(2)) return
        case default
            return
        end select
        if (.not. all(ieee_is_finite(cutoff))) return
        if (any(cutoff <= 0.0_dp) .or. any(cutoff >= rate/2.0_dp)) return
        omega1 = 2.0_dp*rate*tan(pi*cutoff(1)/rate)
        omega2 = 0.0_dp
        if (size(cutoff) == 2) omega2 = 2.0_dp*rate*tan(pi*cutoff(2)/rate)
        if (.not. ieee_is_finite(omega1) .or. .not. ieee_is_finite(omega2)) return
        bw = omega2 - omega1
        center = sqrt(omega1)*sqrt(omega2)
        eps = 1.0_dp / sqrt(10.0_dp**(rs/10.0_dp) - 1.0_dp)
        mu = asinh(1.0_dp/eps)/real(order,dp)
        allocate(p0(order), z0(order))
        nz = 0
        do i = 1, order
            m = real(-order+1 + 2*(i-1), dp)
            theta = pi*m/(2.0_dp*real(order,dp))
            ! Inverse Chebyshev analog low-pass prototype, cutoff at -rs dB.
            root = cmplx(-sinh(mu)*cos(theta), -cosh(mu)*sin(theta), dp)
            p0(i) = 1.0_dp/root
            if (abs(m) > 0.5_dp) then
                nz = nz + 1
                z0(nz) = cmplx(0.0_dp, 1.0_dp/sin(theta), dp)
            end if
        end do
        if (nz /= order-mod(order,2)) return
        n = order
        if (size(cutoff) == 2) n = 2*order
        allocate(poles(n), zeros(n), at_infinity(n))
        at_infinity = .false.
        np = 0
        j = 0
        do i = 1, order
            select case (mode)
            case ('lowpass')
                np = np + 1
                poles(np) = omega1*p0(i)
            case ('highpass')
                np = np + 1
                poles(np) = omega1/p0(i)
            case ('bandpass')
                v = bw*p0(i)
                disc = sqrt(v*v-cmplx(4.0_dp*center**2,0.0_dp,dp))
                poles(np+1) = (v+disc)/2.0_dp
                poles(np+2) = (v-disc)/2.0_dp
                np = np+2
            case ('bandstop')
                v = bw/p0(i)
                disc = sqrt(v*v-cmplx(4.0_dp*center**2,0.0_dp,dp))
                poles(np+1) = (v+disc)/2.0_dp
                poles(np+2) = (v-disc)/2.0_dp
                np = np+2
            end select
        end do
        ! Transform each finite prototype transmission zero.
        do i = 1, nz
            select case (mode)
            case ('lowpass')
                j = j+1
                zeros(j) = omega1*z0(i)
            case ('highpass')
                j = j+1
                zeros(j) = omega1/z0(i)
            case ('bandpass')
                v = bw*z0(i)
                disc = sqrt(v*v-cmplx(4.0_dp*center**2,0.0_dp,dp))
                zeros(j+1) = (v+disc)/2.0_dp
                zeros(j+2) = (v-disc)/2.0_dp
                j = j+2
            case ('bandstop')
                v = bw/z0(i)
                disc = sqrt(v*v-cmplx(4.0_dp*center**2,0.0_dp,dp))
                zeros(j+1) = (v+disc)/2.0_dp
                zeros(j+2) = (v-disc)/2.0_dp
                j = j+2
            end select
        end do
        ! The analog zeros at infinity map to z=-1; bandpass adds
        ! one zero at z=+1 and one at z=-1, while bandstop gets +/-i*w0.
        if (mod(order,2) == 1) then
            select case (mode)
            case ('lowpass', 'highpass')
                j = j+1
                at_infinity(j) = .true.
            case ('bandpass')
                j = j+2
                at_infinity(j-1:j) = .true.
            case ('bandstop')
                j = j+2
                zeros(j-1) = cmplx(0.0_dp,center,dp)
                zeros(j) = conjg(zeros(j-1))
            end select
        end if
        if (j /= n) return
        scale = cmplx(2.0_dp*rate,0.0_dp,dp)
        poles = (scale+poles)/(scale-poles)
        do i = 1,n
            if (at_infinity(i)) then
                if (mode == 'highpass') then
                    zeros(i) = cmplx(1.0_dp,0.0_dp,dp)
                else if (mode == 'bandpass') then
                    if (i == n-1) then
                        zeros(i) = cmplx(1.0_dp,0.0_dp,dp)
                    else
                        zeros(i) = cmplx(-1.0_dp,0.0_dp,dp)
                    end if
                else
                    zeros(i) = cmplx(-1.0_dp,0.0_dp,dp)
                end if
            else
                zeros(i) = (scale+zeros(i))/(scale-zeros(i))
            end if
        end do
        if (.not. all(ieee_is_finite(real(poles,dp)))) return
        if (.not. all(ieee_is_finite(aimag(poles)))) return
        if (.not. all(ieee_is_finite(real(zeros,dp)))) return
        if (.not. all(ieee_is_finite(aimag(zeros)))) return
        if (any(abs(poles) >= 1.0_dp+1.0e-10_dp)) return
        call form_sections(zeros, poles, mode, rate, center, out)
    end function cheby2_vector

    function ellip_scalar(order, rp, rs, cutoff, btype, fs, analog) result(out)
        integer, intent(in) :: order !! elliptic filter prototype order (1 through 20)
        real(dp), intent(in) :: rp !! maximum passband ripple in dB
        real(dp), intent(in) :: rs !! minimum stopband attenuation in dB
        real(dp), intent(in) :: cutoff !! critical frequency
        character(*), intent(in), optional :: btype !! lowpass or highpass
        real(dp), intent(in), optional :: fs !! sampling rate (digital only), default 2
        logical, intent(in), optional :: analog !! true for analog filter
        type(butter_result) :: out
        out = ellip_vector(order, rp, rs, [cutoff], btype, fs, analog)
    end function ellip_scalar

    function ellip_vector(order, rp, rs, cutoff, btype, fs, analog) result(out)
        integer, intent(in) :: order !! elliptic filter prototype order (1 through 20)
        real(dp), intent(in) :: rp !! maximum passband ripple in dB
        real(dp), intent(in) :: rs !! minimum stopband attenuation in dB
        real(dp), intent(in) :: cutoff(:) !! one or two critical frequencies
        character(*), intent(in), optional :: btype !! lowpass, highpass, bandpass, bandstop
        real(dp), intent(in), optional :: fs !! sampling rate (digital only), default 2
        logical, intent(in), optional :: analog !! true for analog filter
        type(butter_result) :: out
        type(elliptic_zpk_result) :: analog_roots
        complex(dp), allocatable :: zd(:), pd(:)
        real(dp), allocatable :: omega(:)
        real(dp), allocatable :: coeff_b(:), coeff_a(:)
        real(dp) :: rate, warped_center, pass_gain
        integer :: n, nz
        logical :: is_analog
        character(:), allocatable :: mode

        is_analog = .false.
        if (present(analog)) is_analog = analog
        if (is_analog .and. present(fs)) return
        mode = 'lowpass'
        if (present(btype)) mode = trim(btype)
        if (is_analog) then
            analog_roots = elliptic_analog_zpk(order,rp,rs,cutoff,mode)
            if (analog_roots%status /= 0) return
            call polynomial_from_roots(analog_roots%zeros,coeff_b)
            call polynomial_from_roots(analog_roots%poles,coeff_a)
            out%b = analog_roots%gain*coeff_b
            out%a = coeff_a
            if (.not. all(ieee_is_finite(out%b))) return
            if (.not. all(ieee_is_finite(out%a))) return
            out%status = design_status_success
            return
        end if
        rate = 2.0_dp
        if (present(fs)) rate = fs
        if (.not. ieee_is_finite(rate) .or. rate <= 0.0_dp) return
        if (.not. all(ieee_is_finite(cutoff))) return
        if (any(cutoff <= 0.0_dp) .or. any(cutoff >= rate/2.0_dp)) return
        omega = 2.0_dp*rate*tan(pi*cutoff/rate)
        if (.not. all(ieee_is_finite(omega))) return
        analog_roots = elliptic_analog_zpk(order,rp,rs,omega,mode)
        if (analog_roots%status /= 0) return
        n = size(analog_roots%poles)
        nz = size(analog_roots%zeros)
        allocate(pd(n),zd(n))
        pd = (2.0_dp*rate+analog_roots%poles)/(2.0_dp*rate-analog_roots%poles)
        zd(:nz) = (2.0_dp*rate+analog_roots%zeros)/(2.0_dp*rate-analog_roots%zeros)
        zd(nz+1:) = cmplx(-1.0_dp,0.0_dp,dp)
        if (.not. all(ieee_is_finite(real(pd,dp)))) return
        if (.not. all(ieee_is_finite(aimag(pd)))) return
        if (.not. all(ieee_is_finite(real(zd,dp)))) return
        if (.not. all(ieee_is_finite(aimag(zd)))) return
        if (any(abs(pd) >= 1.0_dp+1.0e-10_dp)) return
        warped_center = 0.0_dp
        if (size(omega) == 2) warped_center = sqrt(omega(1))*sqrt(omega(2))
        pass_gain = 1.0_dp
        if (mod(order,2) == 0) pass_gain = 10.0_dp**(-rp/20.0_dp)
        call form_sections(zd,pd,mode,rate,warped_center,out,pass_gain)
    end function ellip_vector

    pure subroutine polynomial_from_roots(roots, coeff)
        complex(dp), intent(in) :: roots(:) !! roots of real polynomial
        real(dp), allocatable, intent(out) :: coeff(:) !! descending real coefficients
        complex(dp) :: c(size(roots)+1)
        integer :: i, j
        c = cmplx(0.0_dp,0.0_dp,dp)
        c(1) = cmplx(1.0_dp,0.0_dp,dp)
        do i=1,size(roots)
            do j=i+1,2,-1
                c(j) = c(j)-roots(i)*c(j-1)
            end do
        end do
        coeff = real(c,dp)
    end subroutine polynomial_from_roots

    subroutine form_sections(zeros, poles, mode, fs, center, out, target_gain)
        complex(dp), intent(in) :: zeros(:) !! digital zeros
        complex(dp), intent(in) :: poles(:) !! digital poles
        character(*), intent(in) :: mode !! filter passband selection
        real(dp), intent(in) :: fs !! sampling frequency
        real(dp), intent(in) :: center !! warped center frequency for bandpass
        type(butter_result), intent(inout) :: out !! resulting transfer function and SOS
        real(dp), intent(in), optional :: target_gain !! required gain at reference passband frequency
        real(dp), allocatable :: br(:), ar(:)
        complex(dp) :: q, hv_num, hv_den
        real(dp) :: ref, gain, pcoef(3), zcoef(3)
        integer :: i, n, ns, len
        logical :: ok1, ok2
        n = size(poles)
        if (n < 1 .or. size(zeros) /= n) return
        ns = (n+1)/2
        allocate(out%sos(ns,6))
        out%sos = 0.0_dp
        do i = 1,ns
            call next_factor(poles,i,n,pcoef,ok1)
            call next_factor(zeros,i,n,zcoef,ok2)
            if (.not. ok1 .or. .not. ok2) return
            out%sos(i,4:6) = pcoef
            out%sos(i,1:3) = zcoef
        end do
        ref = 0.0_dp
        if (mode == 'highpass') ref = fs/2.0_dp
        if (mode == 'bandpass') ref = fs/pi * atan(center/(2.0_dp*fs))
        q = cmplx(cos(-2.0_dp*pi*ref/fs),sin(-2.0_dp*pi*ref/fs),dp)
        hv_num = cmplx(1.0_dp,0.0_dp,dp)
        hv_den = cmplx(1.0_dp,0.0_dp,dp)
        do i = 1,ns
            hv_num = hv_num*poly_at(out%sos(i,1:3),q)
            hv_den = hv_den*poly_at(out%sos(i,4:6),q)
        end do
        if (abs(hv_num) <= tiny(1.0_dp)) return
        gain = abs(hv_den/hv_num)
        if (present(target_gain)) gain = gain*target_gain
        if (.not. ieee_is_finite(gain)) return
        out%sos(1,1:3) = gain*out%sos(1,1:3)
        out%b = [1.0_dp]
        out%a = [1.0_dp]
        do i = 1,ns
            len = 3
            if (mod(n,2) == 1 .and. i == ns) len = 2
            br = convolve_coefficients(out%b,out%sos(i,1:len))
            ar = convolve_coefficients(out%a,out%sos(i,4:3+len))
            call move_alloc(br,out%b)
            call move_alloc(ar,out%a)
        end do
        if (.not. all(ieee_is_finite(out%b))) return
        if (.not. all(ieee_is_finite(out%a))) return
        out%status = design_status_success
    contains
        subroutine next_factor(roots, section, ntotal, coeff, ok)
            complex(dp), intent(in) :: roots(:)
            integer, intent(in) :: section !! SOS index
            integer, intent(in) :: ntotal !! number of roots
            real(dp), intent(out) :: coeff(3) !! ascending real polynomial coefficients
            logical, intent(out) :: ok !! root grouping succeeded
            ! Compute the unmarked roots for each section deterministically.
            ! Calls per SOS repeat the preceding selections; n<=40.
            logical :: chosen(size(roots))
            integer :: t, x, y, j
            complex(dp) :: z1, z2
            chosen = .false.
            ok = .false.
            coeff = 0.0_dp
            do t=1,section
                x = 0
                y = 0
                if (t == (ntotal+1)/2 .and. mod(ntotal,2) == 1) then
                    do j=1,ntotal
                        if (chosen(j)) cycle
                        if (abs(aimag(roots(j))) <= 1.0e-7_dp) then
                            x = j
                            exit
                        end if
                    end do
                    if (x == 0) return
                    chosen(x) = .true.
                    if (t == section) then
                        coeff = [1.0_dp, -real(roots(x),dp), 0.0_dp]
                        ok = .true.
                    end if
                else
                    do j=1,ntotal
                        if (chosen(j)) cycle
                        ! Preserve the unmatched real root for the final
                        ! first-order section when the degree is odd.
                        if (mod(ntotal,2)==1 .and. abs(aimag(roots(j)))<=1.0e-7_dp) cycle
                        x = j
                        exit
                    end do
                    if (x == 0) return
                    do j=1,ntotal
                        if (chosen(j) .or. j == x) cycle
                        if (abs(aimag(roots(x))) < 1.0e-7_dp) then
                            if (abs(aimag(roots(j))) > 1.0e-7_dp) cycle
                        else
                            if (abs(roots(j)-conjg(roots(x))) > 1.0e-6_dp) cycle
                        end if
                        y = j
                        exit
                    end do
                    if (y == 0) return
                    chosen(x) = .true.
                    chosen(y) = .true.
                    if (t == section) then
                        z1 = roots(x)
                        z2 = roots(y)
                        coeff = [1.0_dp, -real(z1+z2,dp), real(z1*z2,dp)]
                        ok = .true.
                    end if
                end if
            end do
        end subroutine next_factor
    end subroutine form_sections

    pure complex(dp) function poly_at(a, q) result(y)
        real(dp), intent(in) :: a(:) !! ascending polynomial coefficients
        complex(dp), intent(in) :: q !! unit-circle evaluation point
        integer :: j
        y = (0.0_dp,0.0_dp)
        do j=size(a),1,-1
            y = y*q+a(j)
        end do
    end function poly_at

    pure function convolve_coefficients(a,b) result(c)
        real(dp), intent(in) :: a(:) !! left polynomial
        real(dp), intent(in) :: b(:) !! right polynomial
        real(dp) :: c(size(a)+size(b)-1)
        integer :: i,j
        c = 0.0_dp
        do j=1,size(b)
            do i=1,size(a)
                c(i+j-1) = c(i+j-1)+a(i)*b(j)
            end do
        end do
    end function convolve_coefficients

    function iirfilter_scalar(order, cutoff, btype, ftype, rp, rs, fs, analog) result(out)
        integer, intent(in) :: order !! prototype order
        real(dp), intent(in) :: cutoff !! digital critical frequency
        character(*), intent(in), optional :: btype !! lowpass/highpass
        character(*), intent(in), optional :: ftype !! butter, cheby1, cheby2, ellip, bessel
        real(dp), intent(in), optional :: rp !! passband ripple for cheby1 or ellip
        real(dp), intent(in), optional :: rs !! stopband attenuation for cheby2 or ellip
        real(dp), intent(in), optional :: fs !! sampling rate
        logical, intent(in), optional :: analog !! analog s-domain design, default false
        type(butter_result) :: out
        out = iirfilter_vector(order,[cutoff],btype,ftype,rp,rs,fs,analog)
    end function iirfilter_scalar

    function iirfilter_vector(order, cutoff, btype, ftype, rp, rs, fs, analog) result(out)
        integer, intent(in) :: order !! prototype order
        real(dp), intent(in) :: cutoff(:) !! critical frequency array
        character(*), intent(in), optional :: btype !! filter type
        character(*), intent(in), optional :: ftype !! butter, cheby1, cheby2, ellip, bessel
        real(dp), intent(in), optional :: rp !! Chebyshev-I ripple in dB
        real(dp), intent(in), optional :: rs !! Chebyshev-II stopband attenuation in dB
        real(dp), intent(in), optional :: fs !! sampling frequency
        logical, intent(in), optional :: analog !! analog s-domain design, default false
        type(butter_result) :: out
        character(:), allocatable :: family
        family = 'butter'
        if (present(ftype)) family = trim(ftype)
        select case(family)
        case('butter', 'butterworth')
            out = butter(order,cutoff,btype=btype,fs=fs,analog=analog)
        case('cheby1', 'chebyshev1')
            if (.not. present(rp)) return
            out = cheby1(order,rp,cutoff,btype=btype,fs=fs,analog=analog)
        case('cheby2', 'chebyshev2')
            if (.not. present(rs)) return
            out = cheby2(order,rs,cutoff,btype=btype,fs=fs,analog=analog)
        case('bessel', 'thomson')
            if (present(rp) .or. present(rs)) return
            out = bessel(order,cutoff,btype=btype,fs=fs,analog=analog)
        case('ellip', 'elliptic', 'cauer')
            if (.not. present(rp) .or. .not. present(rs)) return
            out = ellip(order,rp,rs,cutoff,btype=btype,fs=fs,analog=analog)
        end select
    end function iirfilter_vector


    function buttord_scalar(wp, ws, gpass, gstop, fs, analog) result(out)
        real(dp), intent(in) :: wp !! passband edge
        real(dp), intent(in) :: ws !! stopband edge
        real(dp), intent(in) :: gpass !! maximum passband loss in dB
        real(dp), intent(in) :: gstop !! minimum stopband loss in dB
        real(dp), intent(in), optional :: fs !! digital sample rate, default 2
        logical, intent(in), optional :: analog !! analog rad/s frequencies
        type(iir_order_result) :: out
        out = select_iir_order([wp], [ws], gpass, gstop, fs, analog, 'butter')
    end function buttord_scalar

    function buttord_vector(wp, ws, gpass, gstop, fs, analog) result(out)
        real(dp), intent(in) :: wp(:) !! two band passband edges
        real(dp), intent(in) :: ws(:) !! two band stopband edges
        real(dp), intent(in) :: gpass !! maximum passband loss in dB
        real(dp), intent(in) :: gstop !! minimum stopband loss in dB
        real(dp), intent(in), optional :: fs !! digital sample rate, default 2
        logical, intent(in), optional :: analog !! analog rad/s frequencies
        type(iir_order_result) :: out
        out = select_iir_order(wp, ws, gpass, gstop, fs, analog, 'butter')
    end function buttord_vector

    function cheb1ord_scalar(wp, ws, gpass, gstop, fs, analog) result(out)
        real(dp), intent(in) :: wp !! passband edge
        real(dp), intent(in) :: ws !! stopband edge
        real(dp), intent(in) :: gpass !! maximum passband loss in dB
        real(dp), intent(in) :: gstop !! minimum stopband loss in dB
        real(dp), intent(in), optional :: fs !! digital sample rate, default 2
        logical, intent(in), optional :: analog !! analog rad/s frequencies
        type(iir_order_result) :: out
        out = select_iir_order([wp], [ws], gpass, gstop, fs, analog, 'cheby1')
    end function cheb1ord_scalar

    function cheb1ord_vector(wp, ws, gpass, gstop, fs, analog) result(out)
        real(dp), intent(in) :: wp(:) !! two band passband edges
        real(dp), intent(in) :: ws(:) !! two band stopband edges
        real(dp), intent(in) :: gpass !! maximum passband loss in dB
        real(dp), intent(in) :: gstop !! minimum stopband loss in dB
        real(dp), intent(in), optional :: fs !! digital sample rate, default 2
        logical, intent(in), optional :: analog !! analog rad/s frequencies
        type(iir_order_result) :: out
        out = select_iir_order(wp, ws, gpass, gstop, fs, analog, 'cheby1')
    end function cheb1ord_vector

    function cheb2ord_scalar(wp, ws, gpass, gstop, fs, analog) result(out)
        real(dp), intent(in) :: wp !! passband edge
        real(dp), intent(in) :: ws !! stopband edge
        real(dp), intent(in) :: gpass !! maximum passband loss in dB
        real(dp), intent(in) :: gstop !! minimum stopband loss in dB
        real(dp), intent(in), optional :: fs !! digital sample rate, default 2
        logical, intent(in), optional :: analog !! analog rad/s frequencies
        type(iir_order_result) :: out
        out = select_iir_order([wp], [ws], gpass, gstop, fs, analog, 'cheby2')
    end function cheb2ord_scalar

    function cheb2ord_vector(wp, ws, gpass, gstop, fs, analog) result(out)
        real(dp), intent(in) :: wp(:) !! two band passband edges
        real(dp), intent(in) :: ws(:) !! two band stopband edges
        real(dp), intent(in) :: gpass !! maximum passband loss in dB
        real(dp), intent(in) :: gstop !! minimum stopband loss in dB
        real(dp), intent(in), optional :: fs !! digital sample rate, default 2
        logical, intent(in), optional :: analog !! analog rad/s frequencies
        type(iir_order_result) :: out
        out = select_iir_order(wp, ws, gpass, gstop, fs, analog, 'cheby2')
    end function cheb2ord_vector

    function ellipord_scalar(wp, ws, gpass, gstop, fs, analog) result(out)
        real(dp), intent(in) :: wp !! passband edge
        real(dp), intent(in) :: ws !! stopband edge
        real(dp), intent(in) :: gpass !! maximum passband ripple in dB
        real(dp), intent(in) :: gstop !! minimum stopband attenuation in dB
        real(dp), intent(in), optional :: fs !! sampling rate (digital only)
        logical, intent(in), optional :: analog !! analog frequency units, rad/s
        type(iir_order_result) :: out
        out = select_iir_order([wp],[ws],gpass,gstop,fs,analog,'ellip')
    end function ellipord_scalar

    function ellipord_vector(wp, ws, gpass, gstop, fs, analog) result(out)
        real(dp), intent(in) :: wp(:) !! one or two passband edges
        real(dp), intent(in) :: ws(:) !! one or two stopband edges
        real(dp), intent(in) :: gpass !! maximum passband ripple in dB
        real(dp), intent(in) :: gstop !! minimum stopband attenuation in dB
        real(dp), intent(in), optional :: fs !! sampling rate (digital only)
        logical, intent(in), optional :: analog !! analog frequency units, rad/s
        type(iir_order_result) :: out
        out = select_iir_order(wp,ws,gpass,gstop,fs,analog,'ellip')
    end function ellipord_vector

    function select_iir_order(wp, ws, gpass, gstop, fs, analog, family) result(out)
        real(dp), intent(in) :: wp(:) !! scalar or two-element passband array
        real(dp), intent(in) :: ws(:) !! scalar or two-element stopband array
        real(dp), intent(in) :: gpass !! maximum passband attenuation, dB
        real(dp), intent(in) :: gstop !! minimum stopband attenuation, dB
        real(dp), intent(in), optional :: fs !! digital sample rate
        logical, intent(in), optional :: analog !! true for analog rad/s
        character(*), intent(in) :: family !! butter, cheby1, cheby2
        type(iir_order_result) :: out
        real(dp) :: rate, passb(2), stopb(2), nat, gp, gs
        real(dp) :: degree, w0, alpha, delta, product, root
        real(dp) :: bounds(2), critical(2), temp, m0, m1
        integer :: nedge, mode, j
        logical :: is_analog

        nedge = size(wp)
        if (nedge /= 1 .and. nedge /= 2) return
        if (size(ws) /= nedge) return
        if (.not. all(ieee_is_finite(wp)) .or. .not. all(ieee_is_finite(ws))) return
        if (.not. ieee_is_finite(gpass) .or. .not. ieee_is_finite(gstop)) return
        if (gpass <= 0.0_dp .or. gstop <= gpass .or. gstop > 100.0_dp) return
        is_analog = .false.
        if (present(analog)) is_analog = analog
        if (is_analog .and. present(fs)) return
        rate = 2.0_dp
        if (present(fs)) rate = fs
        if (.not. ieee_is_finite(rate) .or. rate <= 0.0_dp) return
        if (any(wp <= 0.0_dp) .or. any(ws <= 0.0_dp)) return
        if (nedge == 2) then
            if (wp(1) >= wp(2) .or. ws(1) >= ws(2)) return
        end if
        if (.not. is_analog) then
            if (any(wp >= rate/2.0_dp) .or. any(ws >= rate/2.0_dp)) return
        end if
        if (nedge == 1) then
            mode = 1
            if (wp(1) > ws(1)) mode = 2
            if (wp(1) == ws(1)) return
        else if (wp(1) > ws(1)) then
            mode = 4 ! bandpass: stop1 < pass1 < pass2 < stop2
            if (.not. (ws(1) < wp(1) .and. wp(2) < ws(2))) return
        else
            mode = 3 ! bandstop: pass1 < stop1 < stop2 < pass2
            if (.not. (wp(1) < ws(1) .and. ws(2) < wp(2))) return
        end if
        passb = 0.0_dp
        stopb = 0.0_dp
        passb(:nedge) = wp
        stopb(:nedge) = ws
        if (.not. is_analog) then
            passb(:nedge) = tan(pi*wp/rate)
            stopb(:nedge) = tan(pi*ws/rate)
            if (.not. all(ieee_is_finite(passb(:nedge)))) return
            if (.not. all(ieee_is_finite(stopb(:nedge)))) return
        end if
        gp = 10.0_dp**(gpass/10.0_dp)-1.0_dp
        gs = 10.0_dp**(gstop/10.0_dp)-1.0_dp
        if (.not. ieee_is_finite(gp) .or. .not. ieee_is_finite(gs)) return
        if (gp <= 0.0_dp .or. gs <= gp) return

        select case (mode)
        case (1)
            nat = stopb(1)/passb(1)
        case (2)
            nat = passb(1)/stopb(1)
        case (3)
            ! The two passband boundaries may move inward to minimize the
            ! prototype order, as in SciPy's bandstop order selection.
            bounds = passb
            do j = 1, 2
                if (j == 1) then
                    passb(1) = optimize_bandstop_edge(1, passb, stopb, bounds(1), stopb(1))
                else
                    passb(2) = optimize_bandstop_edge(2, passb, stopb, stopb(2), bounds(2))
                end if
            end do
            nat = bandstop_nat(passb, stopb)
        case (4)
            nat = minval(abs((stopb**2-passb(1)*passb(2)) / &
                (stopb*(passb(1)-passb(2)))))
        end select
        if (.not. ieee_is_finite(nat)) return
        if (nat <= 1.0_dp) return
        select case (family)
        case ('butter')
            degree = log(gs/gp)/(2.0_dp*log(nat))
        case ('cheby1', 'cheby2')
            degree = acosh(sqrt(gs/gp))/acosh(nat)
        case ('ellip')
            m0 = 1.0_dp/(nat*nat)
            m1 = gp/gs
            degree = elliptic_k(m0)*elliptic_k_complement(m1) / &
                (elliptic_k_complement(m0)*elliptic_k(m1))
        case default
            return
        end select
        if (.not. ieee_is_finite(degree)) return
        if (degree > 20.0_dp) return
        out%order = max(1, ceiling(degree))
        product = passb(1)*passb(2)
        delta = passb(2)-passb(1)
        select case (family)
        case ('butter')
            w0 = gp**(-1.0_dp/(2.0_dp*real(out%order,dp)))
            select case (mode)
            case (1)
                critical(1) = w0*passb(1)
            case (2)
                critical(1) = passb(1)/w0
            case (3)
                root = sqrt(delta**2+4.0_dp*w0**2*product)
                critical(2) = (delta+root)/(2.0_dp*w0)
                critical(1) = product/critical(2)
            case (4)
                temp = w0*delta/2.0_dp
                root = sqrt(temp**2+product)
                critical(1) = product/(root+temp)
                critical(2) = root+temp
            end select
        case ('cheby1', 'ellip')
            critical(:nedge) = passb(:nedge)
        case ('cheby2')
            alpha = acosh(sqrt(gs/gp))/real(out%order,dp)
            w0 = 1.0_dp/cosh(alpha)
            select case (mode)
            case (1)
                critical(1) = passb(1)/w0
            case (2)
                critical(1) = passb(1)*w0
            case (3)
                temp = w0*delta/2.0_dp
                root = sqrt(temp**2+product)
                critical(1) = product/(root+temp)
                critical(2) = root+temp
            case (4)
                temp = delta/(2.0_dp*w0)
                root = sqrt(temp**2+product)
                critical(1) = product/(root+temp)
                critical(2) = root+temp
            end select
        end select
        if (.not. all(ieee_is_finite(critical(:nedge)))) return
        if (any(critical(:nedge) <= 0.0_dp)) return
        if (nedge == 2) then
            if (critical(1) >= critical(2)) return
        end if
        if (.not. is_analog) critical(:nedge) = rate/pi*atan(critical(:nedge))
        if (.not. is_analog) then
            if (any(critical(:nedge) >= rate/2.0_dp)) return
        end if
        out%cutoff = critical(1)
        out%cutoffs = critical(:nedge)
        out%status = design_status_success
    end function select_iir_order

    pure real(dp) function bandstop_nat(passb, stopb) result(nat)
        real(dp), intent(in) :: passb(2) !! warped passband boundaries
        real(dp), intent(in) :: stopb(2) !! warped stopband boundaries
        nat = minval(abs(stopb*(passb(1)-passb(2))/ &
            (stopb**2-passb(1)*passb(2))))
    end function bandstop_nat

    pure real(dp) function bandstop_cost(edge, index, passb, stopb) result(cost)
        real(dp), intent(in) :: edge !! candidate passband boundary
        integer, intent(in) :: index !! which boundary to adjust
        real(dp), intent(in) :: passb(2) !! current passband boundaries
        real(dp), intent(in) :: stopb(2) !! stopband boundaries
        real(dp) :: candidate(2), natural
        candidate = passb
        candidate(index) = edge
        natural = bandstop_nat(candidate, stopb)
        cost = huge(1.0_dp)
        if (.not. ieee_is_finite(natural)) return
        if (natural > 1.0_dp) cost = 1.0_dp/log(natural)
    end function bandstop_cost

    real(dp) function optimize_bandstop_edge(index, passb, stopb, lower, upper) result(edge)
        integer, intent(in) :: index !! index of optimizable passband boundary
        real(dp), intent(in) :: passb(2) !! current passband boundaries
        real(dp), intent(in) :: stopb(2) !! stopband boundaries
        real(dp), intent(in) :: lower !! strict lower search boundary
        real(dp), intent(in) :: upper !! strict upper search boundary
        real(dp), parameter :: golden = 0.3819660112501051518_dp
        real(dp) :: a, b, c, d, fc, fd, epsilon
        integer :: it
        a = lower
        b = upper
        c = a+golden*(b-a)
        d = b-golden*(b-a)
        fc = bandstop_cost(c,index,passb,stopb)
        fd = bandstop_cost(d,index,passb,stopb)
        ! One-dimensional bounded minimization; leave a small safety margin
        ! around stopband boundaries and converge on a scale-relative tolerance.
        epsilon = 2.0e-12_dp*max(1.0_dp,abs(lower),abs(upper))
        do it = 1, 110
            if (b-a <= epsilon) exit
            if (fc <= fd) then
                b = d
                d = c
                fd = fc
                c = a+golden*(b-a)
                fc = bandstop_cost(c,index,passb,stopb)
            else
                a = c
                c = d
                fc = fd
                d = b-golden*(b-a)
                fd = bandstop_cost(d,index,passb,stopb)
            end if
        end do
        edge = 0.5_dp*(a+b)
    end function optimize_bandstop_edge

    function iirdesign_scalar(wp, ws, gpass, gstop, ftype, fs, analog) result(out)
        real(dp), intent(in) :: wp !! passband edge
        real(dp), intent(in) :: ws !! stopband edge
        real(dp), intent(in) :: gpass !! maximum passband loss in dB
        real(dp), intent(in) :: gstop !! minimum stopband loss in dB
        character(*), intent(in), optional :: ftype !! butter, cheby1, cheby2
        real(dp), intent(in), optional :: fs !! digital sample rate
        logical, intent(in), optional :: analog !! analog rad/s specification
        type(butter_result) :: out
        out = iirdesign_vector([wp],[ws],gpass,gstop,ftype,fs,analog)
    end function iirdesign_scalar

    function iirdesign_vector(wp, ws, gpass, gstop, ftype, fs, analog) result(out)
        real(dp), intent(in) :: wp(:) !! one or two passband boundaries
        real(dp), intent(in) :: ws(:) !! one or two stopband boundaries
        real(dp), intent(in) :: gpass !! maximum passband loss in dB
        real(dp), intent(in) :: gstop !! minimum stopband loss in dB
        character(*), intent(in), optional :: ftype !! butter, cheby1, cheby2
        real(dp), intent(in), optional :: fs !! digital sample rate
        logical, intent(in), optional :: analog !! analog rad/s specification
        type(butter_result) :: out
        type(iir_order_result) :: selection
        character(:), allocatable :: family, mode
        family = 'butter'
        if (present(ftype)) family = trim(ftype)
        select case (family)
        case ('butter', 'butterworth')
            selection = buttord(wp,ws,gpass,gstop,fs,analog)
        case ('cheby1', 'chebyshev1')
            selection = cheb1ord(wp,ws,gpass,gstop,fs,analog)
        case ('cheby2', 'chebyshev2')
            selection = cheb2ord(wp,ws,gpass,gstop,fs,analog)
        case ('ellip', 'elliptic', 'cauer')
            selection = ellipord(wp,ws,gpass,gstop,fs,analog)
        case default
            return
        end select
        if (selection%status /= design_status_success) return
        mode = 'lowpass'
        if (size(wp) == 1) then
            if (wp(1) > ws(1)) mode = 'highpass'
        else if (wp(1) > ws(1)) then
            mode = 'bandpass'
        else
            mode = 'bandstop'
        end if
        out = iirfilter(selection%order, selection%cutoffs, btype=mode, ftype=family, &
            rp=gpass, rs=gstop, fs=fs, analog=analog)
    end function iirdesign_vector

    pure function ascii_lower(word) result(lower)
        character(*), intent(in) :: word !! ASCII text to fold to lower case
        character(len(word)) :: lower
        integer :: i, c
        lower = word
        do i=1,len(lower)
            c = iachar(lower(i:i))
            if (c >= iachar('A') .and. c <= iachar('Z')) lower(i:i) = achar(c+32)
        end do
    end function ascii_lower

    function iircomb(w0, qfactor, ftype, fs, pass_zero) result(out)
        real(dp), intent(in) :: w0 !! fundamental frequency, strictly below Nyquist
        real(dp), intent(in) :: qfactor !! ratio of fundamental frequency to bandwidth
        character(*), intent(in), optional :: ftype !! notch (default) or peak
        real(dp), intent(in), optional :: fs !! sampling frequency, default 2
        logical, intent(in), optional :: pass_zero !! whether the response passes DC
        type(butter_result) :: out
        real(dp) :: rate, beta, ax, bx, cx, omega, ratio, sign_coef
        integer :: n
        logical :: peak, dc
        character(:), allocatable :: mode

        rate = 2.0_dp
        if (present(fs)) rate = fs
        if (.not. ieee_is_finite(rate) .or. .not. ieee_is_finite(w0)) return
        if (.not. ieee_is_finite(qfactor)) return
        if (rate <= 0.0_dp .or. w0 <= 0.0_dp .or. w0 >= rate / 2.0_dp) return
        if (qfactor <= 0.0_dp) return
        mode = 'notch'
        if (present(ftype)) mode = ascii_lower(trim(ftype))
        select case (mode)
        case ('notch')
            peak = .false.
        case ('peak')
            peak = .true.
        case default
            return
        end select
        ratio = rate / w0
        ! Cap allocation and avoid overflow in nint for extremely small w0.
        if (.not. ieee_is_finite(ratio)) return
        if (ratio > 1000000.0_dp) return
        n = nint(ratio)
        if (n <= 2) return
        if (abs(w0 - rate / real(n,dp)) / rate > 1.0e-14_dp) return
        omega = 2.0_dp * pi * w0 / rate
        beta = tan(real(n,dp) * omega / (4.0_dp * qfactor))
        if (.not. ieee_is_finite(beta)) return
        if (abs(1.0_dp + beta) < 100.0_dp * epsilon(beta)) return
        ! A positive beta gives the stable, physical passband bandwidth.
        if (beta <= 0.0_dp) return
        ax = (1.0_dp - beta) / (1.0_dp + beta)
        if (peak) then
            bx = beta / (1.0_dp + beta)
            cx = -bx
        else
            bx = 1.0_dp / (1.0_dp + beta)
            cx = bx
        end if
        dc = .false.
        if (present(pass_zero)) dc = pass_zero
        sign_coef = 1.0_dp
        if ((peak .and. dc) .or. (.not. peak .and. .not. dc)) sign_coef = -1.0_dp
        allocate(out%b(n+1), out%a(n+1))
        out%b = 0.0_dp
        out%a = 0.0_dp
        out%b(1) = bx
        out%b(n+1) = sign_coef * cx
        out%a(1) = 1.0_dp
        out%a(n+1) = sign_coef * ax
        out%status = design_status_success
    end function iircomb

    function gammatone(freq, ftype, order, numtaps, fs) result(out)
        real(dp), intent(in) :: freq !! center frequency in units of fs
        character(*), intent(in) :: ftype !! fir or iir
        integer, intent(in), optional :: order !! FIR order, default 4
        integer, intent(in), optional :: numtaps !! FIR length, default max(int(fs*0.015),15)
        real(dp), intent(in), optional :: fs !! sample rate, default 2
        type(butter_result) :: out
        real(dp) :: rate, bw, t, scale, twopi, period, phase, damp, gain
        real(dp) :: g3, g4, cr
        complex(dp) :: g, g1, g2, g5, denom
        integer :: taps, norder, i
        character(:), allocatable :: mode

        rate = 2.0_dp
        if (present(fs)) rate = fs
        if (.not. ieee_is_finite(rate) .or. .not. ieee_is_finite(freq)) return
        if (rate <= 0.0_dp .or. freq <= 0.0_dp .or. freq >= rate / 2.0_dp) return
        mode = ascii_lower(trim(ftype))
        if (mode /= 'fir' .and. mode /= 'iir') return
        twopi = 2.0_dp * pi
        ! Glasberg/Moore ERB fit, used in Slaney's auditory filter bank.
        bw = 1.019_dp * (freq / 9.26449_dp + 24.7_dp)
        if (mode == 'fir') then
            norder = 4
            if (present(order)) norder = order
            if (norder < 1 .or. norder > 24) return
            ! SciPy uses max(int(fs*0.015),15); cap runaway allocations.
            if (rate > 6.66666666e7_dp) return
            taps = max(int(rate * 0.015_dp), 15)
            if (present(numtaps)) taps = numtaps
            if (taps < 1 .or. taps > 1000000) return
            allocate(out%b(taps), out%a(1))
            out%a = 1.0_dp
            scale = 2.0_dp * (twopi * bw)**norder / gamma(real(norder,dp)) / rate
            do i=1,taps
                t = real(i-1,dp) / rate
                out%b(i) = t**(norder-1) * exp(-twopi*bw*t) * cos(twopi*freq*t) * scale
            end do
            if (.not. all(ieee_is_finite(out%b))) then
                deallocate(out%b, out%a)
                return
            end if
        else
            ! Four gammatone stages, 8th-order real IIR transfer function.
            period = 1.0_dp / rate
            bw = twopi * bw
            phase = twopi * freq * period
            damp = exp(-bw * period)
            cr = cos(phase)
            g3 = sqrt(3.0_dp + 2.0_dp**1.5_dp) * sin(phase)
            g4 = sqrt(3.0_dp - 2.0_dp**1.5_dp) * sin(phase)
            g1 = -2.0_dp * period * exp(cmplx(0.0_dp,2.0_dp*phase,dp))
            g2 = 2.0_dp * period * damp * exp(cmplx(0.0_dp,phase,dp))
            g5 = exp(cmplx(0.0_dp,2.0_dp*phase,dp))
            denom = (-2.0_dp*damp*damp - 2.0_dp*g5 + 2.0_dp*(1.0_dp+g5)*damp)**4
            if (abs(denom) <= tiny(1.0_dp)) return
            g = (g1+g2*(cr-g4)) * (g1+g2*(cr+g4)) &
                * (g1+g2*(cr-g3)) * (g1+g2*(cr+g3)) / denom
            gain = abs(g)
            if (.not. ieee_is_finite(gain)) return
            if (gain <= tiny(1.0_dp)) return
            allocate(out%b(5), out%a(9))
            out%b = period**4 / gain * [ &
                1.0_dp, -4.0_dp*cr*damp, 6.0_dp*cos(2.0_dp*phase)*damp**2, &
                -4.0_dp*cos(3.0_dp*phase)*damp**3, cos(4.0_dp*phase)*damp**4]
            out%a = [ &
                1.0_dp, -8.0_dp*cr*damp, &
                4.0_dp*(4.0_dp+3.0_dp*cos(2.0_dp*phase))*damp**2, &
                -8.0_dp*(6.0_dp*cr+cos(3.0_dp*phase))*damp**3, &
                2.0_dp*(18.0_dp+16.0_dp*cos(2.0_dp*phase)+cos(4.0_dp*phase))*damp**4, &
                -8.0_dp*(6.0_dp*cr+cos(3.0_dp*phase))*damp**5, &
                4.0_dp*(4.0_dp+3.0_dp*cos(2.0_dp*phase))*damp**6, &
                -8.0_dp*cr*damp**7, damp**8]
            if (.not. all(ieee_is_finite(out%b)) .or. .not. all(ieee_is_finite(out%a))) then
                deallocate(out%b, out%a)
                return
            end if
        end if
        out%status = design_status_success
    end function gammatone

    function iirnotch(w0, qfactor, fs) result(out)
        real(dp), intent(in) :: w0 !! notch center frequency
        real(dp), intent(in) :: qfactor !! center frequency divided by -3dB bandwidth
        real(dp), intent(in), optional :: fs !! sampling rate, default 2
        type(butter_result) :: out
        out = resonator(w0,qfactor,fs,.false.)
    end function iirnotch

    function iirpeak(w0, qfactor, fs) result(out)
        real(dp), intent(in) :: w0 !! resonant center frequency
        real(dp), intent(in) :: qfactor !! center frequency divided by -3dB bandwidth
        real(dp), intent(in), optional :: fs !! sampling rate, default 2
        type(butter_result) :: out
        out = resonator(w0,qfactor,fs,.true.)
    end function iirpeak

    function resonator(w0,qfactor,fs,peak) result(out)
        real(dp), intent(in) :: w0 !! notch or peak center
        real(dp), intent(in) :: qfactor !! quality factor
        real(dp), intent(in), optional :: fs !! sampling rate
        logical, intent(in) :: peak !! true for peak, false for notch
        type(butter_result) :: out
        real(dp) :: rate, wn, bw, beta, gain, cosine
        rate = 2.0_dp
        if (present(fs)) rate = fs
        if (.not. ieee_is_finite(w0)) return
        if (.not. ieee_is_finite(qfactor)) return
        if (.not. ieee_is_finite(rate)) return
        if (rate <= 0.0_dp .or. qfactor <= 0.0_dp) return
        if (w0 <= 0.0_dp .or. w0 >= rate/2.0_dp) return
        wn = 2.0_dp*pi*w0/rate
        bw = wn/qfactor
        ! Reject widths whose tangent would give unstable digital biquads.
        if (bw <= 0.0_dp .or. bw >= pi) return
        beta = tan(bw/2.0_dp)
        gain = 1.0_dp/(1.0_dp+beta)
        cosine = cos(wn)
        if (peak) then
            out%b = (1.0_dp-gain)*[1.0_dp, 0.0_dp, -1.0_dp]
        else
            out%b = gain*[1.0_dp, -2.0_dp*cosine, 1.0_dp]
        end if
        out%a = [1.0_dp,-2.0_dp*gain*cosine,2.0_dp*gain-1.0_dp]
        allocate(out%sos(1,6))
        out%sos(1,1:3) = out%b
        out%sos(1,4:6) = out%a
        out%status = design_status_success
    end function resonator
end module scifort_signal_iir_extra
