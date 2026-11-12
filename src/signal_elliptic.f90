! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Elliptic (Cauer) analog IIR prototypes and frequency transformations.
! Real Jacobi elliptic functions use the descending arithmetic-geometric
! mean; complete integrals use the AGM, and the prototype modulus is found
! from the elliptic degree equation. See CODE_PROVENANCE.md.
module scifort_signal_elliptic
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    implicit none
    private
    real(dp), parameter :: pi = acos(-1.0_dp)

    type, public :: elliptic_zpk_result
        complex(dp), allocatable :: zeros(:) !! finite analog transmission zeros
        complex(dp), allocatable :: poles(:) !! analog filter poles
        real(dp) :: gain = 0.0_dp !! analog gain multiplying product(s-z)/product(s-p)
        integer :: status = 1 !! zero on success
    end type elliptic_zpk_result

    public :: elliptic_analog_zpk, elliptic_k, elliptic_k_complement

contains

    pure real(dp) function elliptic_k(m) result(value)
        real(dp), intent(in) :: m !! parameter (modulus squared), 0 <= m < 1
        real(dp) :: a, b, anew
        integer :: i
        value = huge(1.0_dp)
        if (.not. ieee_is_finite(m)) return
        if (m < 0.0_dp .or. m >= 1.0_dp) return
        a = 1.0_dp
        b = sqrt(1.0_dp-m)
        do i = 1, 80
            anew = (a+b)/2.0_dp
            b = sqrt(a*b)
            a = anew
            if (abs(a-b) <= 2.0_dp*epsilon(a)*a) exit
        end do
        value = pi/(2.0_dp*a)
    end function elliptic_k

    pure real(dp) function elliptic_k_complement(m) result(value)
        real(dp), intent(in) :: m !! complementary parameter, 0 < m <= 1
        real(dp) :: a, b, anew
        integer :: i
        value = huge(1.0_dp)
        if (.not. ieee_is_finite(m)) return
        if (m <= 0.0_dp .or. m > 1.0_dp) return
        a = 1.0_dp
        b = sqrt(m)
        do i = 1, 80
            anew = (a+b)/2.0_dp
            b = sqrt(a*b)
            a = anew
            if (abs(a-b) <= 2.0_dp*epsilon(a)*a) exit
        end do
        value = pi/(2.0_dp*a)
    end function elliptic_k_complement

    pure subroutine real_jacobi(u, m_complement, sn, cn, dn)
        real(dp), intent(in) :: u !! real Jacobi elliptic argument
        real(dp), intent(in) :: m_complement !! 1 - elliptic parameter, without cancellation
        real(dp), intent(out) :: sn !! sn(u|m)
        real(dp), intent(out) :: cn !! cn(u|m)
        real(dp), intent(out) :: dn !! dn(u|m)
        real(dp) :: a(0:32), c(0:32), b, phi, arg
        integer :: i, j, n
        a(0) = 1.0_dp
        b = sqrt(m_complement)
        c(0) = sqrt(max(0.0_dp, 1.0_dp-m_complement))
        n = 0
        do i = 1, 32
            c(i) = (a(i-1)-b)/2.0_dp
            a(i) = (a(i-1)+b)/2.0_dp
            b = sqrt(a(i-1)*b)
            n = i
            if (abs(c(i)) <= 2.0_dp*epsilon(b)*a(i)) exit
        end do
        phi = scale(a(n)*u,n)
        do j = n, 1, -1
            arg = c(j)*sin(phi)/a(j)
            phi = 0.5_dp*(phi+asin(max(-1.0_dp,min(1.0_dp,arg))))
        end do
        sn = sin(phi)
        cn = cos(phi)
        dn = sqrt(cn*cn+m_complement*sn*sn)
    end subroutine real_jacobi

    ! Carlson's duplication formula for RF(x,y,z), x,y,z >= 0.
    pure real(dp) function carlson_rf(x, y, z) result(value)
        real(dp), intent(in) :: x !! nonnegative x
        real(dp), intent(in) :: y !! nonnegative y
        real(dp), intent(in) :: z !! positive z
        real(dp) :: a, b, c, mean, sx, sy, sz, e2, e3, lam
        integer :: i
        a = x
        b = y
        c = z
        do i = 1, 80
            mean = (a+b+c)/3.0_dp
            sx = (mean-a)/mean
            sy = (mean-b)/mean
            sz = (mean-c)/mean
            if (max(abs(sx), abs(sy), abs(sz)) < 0.001_dp) exit
            lam = sqrt(a*b)+sqrt(a*c)+sqrt(b*c)
            a = (a+lam)/4.0_dp
            b = (b+lam)/4.0_dp
            c = (c+lam)/4.0_dp
        end do
        e2 = sx*sy+sy*sz+sz*sx
        e3 = sx*sy*sz
        value = (1.0_dp-e2/10.0_dp+e3/14.0_dp+e2*e2/24.0_dp &
            -3.0_dp*e2*e3/44.0_dp-5.0_dp*e2**3/208.0_dp) / sqrt(mean)
    end function carlson_rf

    ! Compute m and 1-m separately to retain precision when m approaches 1.
    pure subroutine elliptic_degree(order, m1, m, mc, ok)
        integer, intent(in) :: order !! lowpass prototype order, >= 2
        real(dp), intent(in) :: m1 !! ripple ratio in (0,1)
        real(dp), intent(out) :: m !! prototype elliptic parameter
        real(dp), intent(out) :: mc !! complementary parameter 1-m
        logical, intent(out) :: ok !! success
        real(dp) :: k1, k1p, q, numerator, denominator, loglo, loghi, logmid, cand
        real(dp) :: target, ratio, step
        integer :: j, i
        ok = .false.
        m = 0.0_dp
        mc = 1.0_dp
        if (m1 <= 0.0_dp .or. m1 >= 1.0_dp) return
        k1 = elliptic_k(m1)
        k1p = elliptic_k_complement(m1)
        q = exp(-pi*k1p/(real(order,dp)*k1))
        numerator = 0.0_dp
        denominator = 1.0_dp
        do j = 0, 99
            step = q**(j*(j+1))
            numerator = numerator+step
            if (j > 0) denominator = denominator+2.0_dp*q**(j*j)
            if (step < 1.0e-24_dp .and. j > 4) exit
        end do
        m = 16.0_dp*q*(numerator/denominator)**4
        if (m < 0.98_dp .and. m > 0.0_dp) then
            mc = 1.0_dp-m
            ok = .true.
            return
        end if
        ! Close to the singularity, the nome formula loses digits in 1-m.
        ! Solve the same degree equation using log(1-m) as the variable.
        target = k1p/(real(order,dp)*k1)
        loglo = log(tiny(1.0_dp))
        loghi = log(0.5_dp)
        do i = 1, 120
            logmid = (loglo+loghi)/2.0_dp
            cand = exp(logmid)
            ratio = elliptic_k(cand)/elliptic_k_complement(cand)
            if (ratio > target) then
                loghi = logmid
            else
                loglo = logmid
            end if
        end do
        mc = exp((loglo+loghi)/2.0_dp)
        m = 1.0_dp-mc
        ok = mc > 0.0_dp .and. mc < 1.0_dp
    end subroutine elliptic_degree

    function elliptic_analog_zpk(order, rp, rs, cutoff, btype) result(out)
        integer, intent(in) :: order !! prototype order, 1 through 20
        real(dp), intent(in) :: rp !! passband ripple in dB
        real(dp), intent(in) :: rs !! stopband attenuation in dB
        real(dp), intent(in) :: cutoff(:) !! one or two critical angular frequencies
        character(*), intent(in), optional :: btype !! lowpass (default), highpass, bandpass, bandstop
        type(elliptic_zpk_result) :: out
        complex(dp) :: prototype_p(20), prototype_z(20), p, v, disc
        complex(dp) :: numerator, denominator
        real(dp) :: m, mc, k, k1, m1, eps, epssq, r, v0, sv, cv, dv
        real(dp) :: sn, cn, dn, bw, center, phi, gain
        integer :: i, j, nz, np, degree, nfinal, kidx
        logical :: ok
        character(:), allocatable :: mode

        if (order < 1 .or. order > 20) return
        if (.not. ieee_is_finite(rp) .or. .not. ieee_is_finite(rs)) return
        if (rp <= 0.0_dp .or. rs <= rp .or. rs > 100.0_dp) return
        if (size(cutoff) < 1 .or. size(cutoff) > 2) return
        if (.not. all(ieee_is_finite(cutoff))) return
        if (any(cutoff <= 0.0_dp)) return
        mode = 'lowpass'
        if (present(btype)) mode = trim(btype)
        select case (mode)
        case ('lowpass', 'low', 'highpass', 'high')
            if (size(cutoff) /= 1) return
        case ('bandpass', 'band', 'bandstop', 'stop')
            if (size(cutoff) /= 2) return
            if (cutoff(1) >= cutoff(2)) return
        case default
            return
        end select
        epssq = expm1_local(log(10.0_dp)*rp/10.0_dp)
        eps = sqrt(epssq)
        m1 = epssq/expm1_local(log(10.0_dp)*rs/10.0_dp)
        nz = 0
        np = 0
        if (order == 1) then
            prototype_p(1) = cmplx(-1.0_dp/eps,0.0_dp,dp)
            np = 1
        else
            call elliptic_degree(order,m1,m,mc,ok)
            if (.not. ok) return
            k = elliptic_k_complement(mc)
            k1 = elliptic_k(m1)
            ! sc(v|1-m1)=1/eps, with v=F(atan(1/eps)|1-m1).
            phi = atan2(1.0_dp,eps)
            r = sin(phi)*carlson_rf(cos(phi)**2,cos(phi)**2+m1*sin(phi)**2,1.0_dp)
            v0 = k*r/(real(order,dp)*k1)
            call real_jacobi(v0, m, sv, cv, dv)
            do j = 1-mod(order,2), order-1, 2
                call real_jacobi(real(j,dp)*k/real(order,dp),mc,sn,cn,dn)
                if (j > 0) then
                    nz = nz+1
                    prototype_z(nz) = cmplx(0.0_dp,1.0_dp/(sqrt(m)*sn),dp)
                    nz = nz+1
                    prototype_z(nz) = conjg(prototype_z(nz-1))
                end if
                p = -(cn*dn*sv*cv+cmplx(0.0_dp,sn*dv,dp)) / &
                    (cv*cv+m*sn*sn*sv*sv)
                np = np+1
                prototype_p(np) = p
                if (j > 0) then
                    np = np+1
                    prototype_p(np) = conjg(p)
                end if
            end do
        end if
        if (np /= order .or. nz /= order-mod(order,2)) return
        numerator = cmplx(1.0_dp,0.0_dp,dp)
        denominator = cmplx(1.0_dp,0.0_dp,dp)
        do i = 1, order
            numerator = numerator*(-prototype_p(i))
        end do
        do i = 1, nz
            denominator = denominator*(-prototype_z(i))
        end do
        gain = real(numerator/denominator,dp)
        if (mod(order,2) == 0) gain = gain/sqrt(1.0_dp+epssq)
        degree = order-nz
        nfinal = order
        if (size(cutoff) == 2) nfinal = 2*order
        allocate(out%poles(nfinal))
        bw = 0.0_dp
        center = 0.0_dp
        if (size(cutoff) == 2) then
            bw = cutoff(2)-cutoff(1)
            center = sqrt(cutoff(1))*sqrt(cutoff(2))
        end if
        select case (mode)
        case ('low', 'lowpass')
            out%poles = cutoff(1)*prototype_p(:order)
            allocate(out%zeros(nz))
            out%zeros = cutoff(1)*prototype_z(:nz)
            gain = gain*cutoff(1)**degree
        case ('high', 'highpass')
            out%poles = cutoff(1)/prototype_p(:order)
            allocate(out%zeros(order))
            out%zeros(:nz) = cutoff(1)/prototype_z(:nz)
            out%zeros(nz+1:) = cmplx(0.0_dp,0.0_dp,dp)
            gain = gain*real(denominator/numerator,dp)
        case ('band', 'bandpass')
            do i = 1, order
                v = bw*prototype_p(i)/2.0_dp
                disc = sqrt(v*v-cmplx(center**2,0.0_dp,dp))
                out%poles(2*i-1) = v+disc
                out%poles(2*i) = v-disc
            end do
            allocate(out%zeros(2*nz+degree))
            do i = 1, nz
                v = bw*prototype_z(i)/2.0_dp
                disc = sqrt(v*v-cmplx(center**2,0.0_dp,dp))
                out%zeros(2*i-1) = v+disc
                out%zeros(2*i) = v-disc
            end do
            out%zeros(2*nz+1:) = cmplx(0.0_dp,0.0_dp,dp)
            gain = gain*bw**degree
        case ('stop', 'bandstop')
            do i = 1, order
                v = bw/(2.0_dp*prototype_p(i))
                disc = sqrt(v*v-cmplx(center**2,0.0_dp,dp))
                out%poles(2*i-1) = v+disc
                out%poles(2*i) = v-disc
            end do
            allocate(out%zeros(2*nz+2*degree))
            do i = 1, nz
                v = bw/(2.0_dp*prototype_z(i))
                disc = sqrt(v*v-cmplx(center**2,0.0_dp,dp))
                out%zeros(2*i-1) = v+disc
                out%zeros(2*i) = v-disc
            end do
            do i = 1, degree
                kidx = 2*nz+2*i-1
                out%zeros(kidx) = cmplx(0.0_dp,center,dp)
                out%zeros(kidx+1) = conjg(out%zeros(kidx))
            end do
            gain = gain*real(denominator/numerator,dp)
        end select
        if (.not. ieee_is_finite(gain)) return
        if (.not. all(ieee_is_finite(real(out%poles,dp)))) return
        if (.not. all(ieee_is_finite(aimag(out%poles)))) return
        if (.not. all(ieee_is_finite(real(out%zeros,dp)))) return
        if (.not. all(ieee_is_finite(aimag(out%zeros)))) return
        out%gain = gain
        out%status = 0
    end function elliptic_analog_zpk

    pure real(dp) function expm1_local(x) result(y)
        real(dp), intent(in) :: x !! argument of exp(x)-1
        if (abs(x) < 0.0001_dp) then
            y = x*(1.0_dp+x*(0.5_dp+x*(1.0_dp/6.0_dp+x/24.0_dp)))
        else
            y = exp(x)-1.0_dp
        end if
    end function expm1_local

end module scifort_signal_elliptic
