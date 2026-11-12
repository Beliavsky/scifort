! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Real Kelvin functions of order zero. This is original SciFort code based on
! the defining Bessel identities and large-argument asymptotic expansions in
! NIST DLMF sections 10.61 and 10.67. No third-party implementation source is
! copied or adapted.

module scifort_kelvin
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_constants, only : scifort_pi
    use scifort_kinds, only : dp
    use scifort_math, only : negative_infinity, positive_infinity, quiet_nan
    implicit none
    private

    integer, parameter, public :: kelvin_status_success = 0
    integer, parameter, public :: kelvin_status_invalid_argument = 1
    integer, parameter, public :: kelvin_status_numerical_failure = 2
    integer, parameter, public :: kelvin_max_zeros = 200

    integer, parameter :: root_ber = 1
    integer, parameter :: root_bei = 2
    integer, parameter :: root_ker = 3
    integer, parameter :: root_kei = 4
    integer, parameter :: root_berp = 5
    integer, parameter :: root_beip = 6
    integer, parameter :: root_kerp = 7
    integer, parameter :: root_keip = 8

    real(dp), parameter :: euler_gamma = &
        0.577215664901532860606512090082402431042_dp
    real(dp), parameter :: inv_sqrt_two = &
        0.707106781186547524400844362104849039285_dp

    public :: bei
    public :: bei_zeros
    public :: beip
    public :: beip_zeros
    public :: ber
    public :: ber_zeros
    public :: berp
    public :: berp_zeros
    public :: kei
    public :: kei_zeros
    public :: keip
    public :: keip_zeros
    public :: kelvin
    public :: kelvin_zeros
    public :: ker
    public :: ker_zeros
    public :: kerp
    public :: kerp_zeros

contains

    pure elemental function ber(x) result(y)
        real(dp), intent(in) :: x !! real Kelvin argument
        real(dp) :: y

        complex(dp) :: be
        complex(dp) :: bep
        complex(dp) :: ke
        complex(dp) :: kep

        call kelvin(x, be, ke, bep, kep)
        y = real(be, dp)
    end function ber

    pure elemental function bei(x) result(y)
        real(dp), intent(in) :: x !! real Kelvin argument
        real(dp) :: y

        complex(dp) :: be
        complex(dp) :: bep
        complex(dp) :: ke
        complex(dp) :: kep

        call kelvin(x, be, ke, bep, kep)
        y = aimag(be)
    end function bei

    pure elemental function berp(x) result(y)
        real(dp), intent(in) :: x !! real Kelvin argument
        real(dp) :: y

        complex(dp) :: be
        complex(dp) :: bep
        complex(dp) :: ke
        complex(dp) :: kep

        call kelvin(x, be, ke, bep, kep)
        y = real(bep, dp)
    end function berp

    pure elemental function beip(x) result(y)
        real(dp), intent(in) :: x !! real Kelvin argument
        real(dp) :: y

        complex(dp) :: be
        complex(dp) :: bep
        complex(dp) :: ke
        complex(dp) :: kep

        call kelvin(x, be, ke, bep, kep)
        y = aimag(bep)
    end function beip

    pure elemental function ker(x) result(y)
        real(dp), intent(in) :: x !! nonnegative real Kelvin argument
        real(dp) :: y

        complex(dp) :: be
        complex(dp) :: bep
        complex(dp) :: ke
        complex(dp) :: kep

        call kelvin(x, be, ke, bep, kep)
        y = real(ke, dp)
    end function ker

    pure elemental function kei(x) result(y)
        real(dp), intent(in) :: x !! nonnegative real Kelvin argument
        real(dp) :: y

        complex(dp) :: be
        complex(dp) :: bep
        complex(dp) :: ke
        complex(dp) :: kep

        call kelvin(x, be, ke, bep, kep)
        y = aimag(ke)
    end function kei

    pure elemental function kerp(x) result(y)
        real(dp), intent(in) :: x !! nonnegative real Kelvin argument
        real(dp) :: y

        complex(dp) :: be
        complex(dp) :: bep
        complex(dp) :: ke
        complex(dp) :: kep

        call kelvin(x, be, ke, bep, kep)
        y = real(kep, dp)
    end function kerp

    pure elemental function keip(x) result(y)
        real(dp), intent(in) :: x !! nonnegative real Kelvin argument
        real(dp) :: y

        complex(dp) :: be
        complex(dp) :: bep
        complex(dp) :: ke
        complex(dp) :: kep

        call kelvin(x, be, ke, bep, kep)
        y = aimag(kep)
    end function keip

    pure elemental subroutine kelvin(x, be, ke, bep, kep)
        real(dp), intent(in) :: x !! real Kelvin argument; K-family values require x >= 0
        complex(dp), intent(out) :: be !! ber(x) + i*bei(x)
        complex(dp), intent(out) :: ke !! ker(x) + i*kei(x)
        complex(dp), intent(out) :: bep !! derivative ber'(x) + i*bei'(x)
        complex(dp), intent(out) :: kep !! derivative ker'(x) + i*kei'(x)

        complex(dp), parameter :: ray_be = &
            cmplx(-inv_sqrt_two, inv_sqrt_two, kind=dp)
        complex(dp), parameter :: ray_ke = &
            cmplx(inv_sqrt_two, inv_sqrt_two, kind=dp)
        complex(dp) :: j0_value
        complex(dp) :: j1_value
        complex(dp) :: k0_value
        complex(dp) :: k1_value
        complex(dp) :: nan_complex
        real(dp) :: ax
        real(dp) :: sign_x

        nan_complex = cmplx(quiet_nan(0.0_dp), quiet_nan(0.0_dp), kind=dp)
        if (ieee_is_nan(x) .or. .not. ieee_is_finite(x)) then
            be = nan_complex
            ke = nan_complex
            bep = nan_complex
            kep = nan_complex
            return
        end if

        ax = abs(x)
        if (ax == 0.0_dp) then
            be = cmplx(1.0_dp, 0.0_dp, kind=dp)
            bep = cmplx(0.0_dp, 0.0_dp, kind=dp)
            ke = cmplx(positive_infinity(x), -0.25_dp * scifort_pi, kind=dp)
            kep = cmplx(negative_infinity(x), 0.0_dp, kind=dp)
            return
        end if

        if (ax < 24.0_dp) then
            call complex_j01_series(ax * ray_be, j0_value, j1_value)
            be = j0_value
            bep = -ray_be * j1_value
        else
            be = complex_j_asymptotic(0, ax * ray_be, cmplx(1.0_dp, 0.0_dp, kind=dp))
            bep = complex_j_asymptotic(1, ax * ray_be, -ray_be)
        end if

        if (x < 0.0_dp) then
            sign_x = -1.0_dp
            bep = sign_x * bep
            ke = nan_complex
            kep = nan_complex
            return
        end if

        if (ax < 2.0_dp / huge(1.0_dp)) then
            ke = cmplx(-log(ax) + log(2.0_dp) - euler_gamma, &
                -0.25_dp * scifort_pi, kind=dp)
            kep = cmplx(negative_infinity(x), 0.0_dp, kind=dp)
            return
        end if

        if (ax * inv_sqrt_two > 750.0_dp) then
            ke = cmplx(0.0_dp, 0.0_dp, kind=dp)
            kep = cmplx(0.0_dp, 0.0_dp, kind=dp)
            return
        end if

        if (ax < 10.0_dp) then
            call complex_k01_series(ax * ray_ke, k0_value, k1_value)
        else
            k0_value = complex_k_asymptotic(0, ax * ray_ke)
            k1_value = complex_k_asymptotic(1, ax * ray_ke)
        end if
        ke = k0_value
        kep = -ray_ke * k1_value
    end subroutine kelvin

    pure subroutine complex_j01_series(z, j0_value, j1_value)
        complex(dp), intent(in) :: z !! finite complex Bessel argument
        complex(dp), intent(out) :: j0_value !! J_0(z)
        complex(dp), intent(out) :: j1_value !! J_1(z)

        integer, parameter :: max_terms = 256
        integer :: k
        complex(dp) :: t
        complex(dp) :: term0
        complex(dp) :: term1

        t = -0.25_dp * z * z
        term0 = cmplx(1.0_dp, 0.0_dp, kind=dp)
        j0_value = term0
        do k = 1, max_terms
            term0 = term0 * t / real(k * k, dp)
            j0_value = j0_value + term0
            if (abs(term0) <= epsilon(1.0_dp) * max(1.0_dp, abs(j0_value))) exit
        end do

        term1 = 0.5_dp * z
        j1_value = term1
        do k = 1, max_terms
            term1 = term1 * t / real(k * (k + 1), dp)
            j1_value = j1_value + term1
            if (abs(term1) <= epsilon(1.0_dp) * max(1.0_dp, abs(j1_value))) exit
        end do
    end subroutine complex_j01_series

    pure function complex_j_asymptotic(nu, z, multiplier) result(value)
        integer, intent(in) :: nu !! cylindrical-Bessel order, zero or one
        complex(dp), intent(in) :: z !! large complex argument in the upper half-plane
        complex(dp), intent(in) :: multiplier !! finite factor applied before exponential scaling
        complex(dp) :: value

        integer, parameter :: max_terms = 80
        integer :: k
        real(dp) :: factor
        real(dp) :: log_scale
        real(dp) :: previous
        complex(dp) :: base
        complex(dp) :: omega
        complex(dp) :: series
        complex(dp) :: term

        series = cmplx(1.0_dp, 0.0_dp, kind=dp)
        term = cmplx(1.0_dp, 0.0_dp, kind=dp)
        previous = huge(1.0_dp)
        do k = 1, max_terms
            factor = (4.0_dp * real(nu * nu, dp) - real((2 * k - 1) ** 2, dp)) / &
                real(8 * k, dp)
            term = term * cmplx(0.0_dp, -factor, kind=dp) / z
            if (abs(term) > previous) exit
            previous = abs(term)
            series = series + term
        end do

        omega = z - 0.5_dp * real(nu, dp) * scifort_pi - 0.25_dp * scifort_pi
        base = multiplier * sqrt((0.5_dp / scifort_pi) / z) * &
            cmplx(cos(real(omega, dp)), -sin(real(omega, dp)), kind=dp) * series
        log_scale = aimag(omega)
        value = cmplx(scaled_real_component(real(base, dp), log_scale), &
            scaled_real_component(aimag(base), log_scale), kind=dp)
    end function complex_j_asymptotic

    pure function scaled_real_component(component, log_scale) result(value)
        real(dp), intent(in) :: component !! unscaled real or imaginary component
        real(dp), intent(in) :: log_scale !! logarithmic positive scale factor
        real(dp) :: value

        real(dp) :: log_abs
        real(dp) :: magnitude

        if (component == 0.0_dp) then
            value = component
            return
        end if
        log_abs = log(abs(component)) + log_scale
        if (log_abs > log(huge(1.0_dp))) then
            value = sign(positive_infinity(component), component)
        else
            magnitude = exp(log_abs)
            value = sign(magnitude, component)
        end if
    end function scaled_real_component

    pure subroutine complex_k01_series(z, k0_value, k1_value)
        complex(dp), intent(in) :: z !! finite nonzero complex argument on the principal branch
        complex(dp), intent(out) :: k0_value !! principal K_0(z)
        complex(dp), intent(out) :: k1_value !! principal K_1(z)

        integer, parameter :: max_terms = 256
        integer :: k
        real(dp) :: harmonic
        complex(dp) :: dseries
        complex(dp) :: i0_value
        complex(dp) :: i1_value
        complex(dp) :: log_term
        complex(dp) :: series
        complex(dp) :: t
        complex(dp) :: term0
        complex(dp) :: term1
        complex(dp) :: weighted

        t = 0.25_dp * z * z
        term0 = cmplx(1.0_dp, 0.0_dp, kind=dp)
        i0_value = term0
        series = cmplx(0.0_dp, 0.0_dp, kind=dp)
        dseries = cmplx(0.0_dp, 0.0_dp, kind=dp)
        harmonic = 0.0_dp
        do k = 1, max_terms
            term0 = term0 * t / real(k * k, dp)
            i0_value = i0_value + term0
            harmonic = harmonic + 1.0_dp / real(k, dp)
            weighted = harmonic * term0
            series = series + weighted
            dseries = dseries + 2.0_dp * real(k, dp) * weighted / z
            if (abs(term0) <= epsilon(1.0_dp) * max(1.0_dp, abs(i0_value)) .and. &
                    abs(weighted) <= epsilon(1.0_dp) * max(1.0_dp, abs(series))) exit
        end do

        term1 = 0.5_dp * z
        i1_value = term1
        do k = 1, max_terms
            term1 = term1 * t / real(k * (k + 1), dp)
            i1_value = i1_value + term1
            if (abs(term1) <= epsilon(1.0_dp) * max(1.0_dp, abs(i1_value))) exit
        end do

        log_term = log(0.5_dp * z) + euler_gamma
        k0_value = -log_term * i0_value + series
        k1_value = i0_value / z + log_term * i1_value - dseries
    end subroutine complex_k01_series

    pure function complex_k_asymptotic(nu, z) result(value)
        integer, intent(in) :: nu !! modified-Bessel order, zero or one
        complex(dp), intent(in) :: z !! large complex argument in the right half-plane
        complex(dp) :: value

        integer, parameter :: max_terms = 80
        integer :: k
        real(dp) :: factor
        real(dp) :: previous
        complex(dp) :: series
        complex(dp) :: term

        series = cmplx(1.0_dp, 0.0_dp, kind=dp)
        term = cmplx(1.0_dp, 0.0_dp, kind=dp)
        previous = huge(1.0_dp)
        do k = 1, max_terms
            factor = (4.0_dp * real(nu * nu, dp) - real((2 * k - 1) ** 2, dp)) / &
                real(8 * k, dp)
            term = term * factor / z
            if (abs(term) > previous) exit
            previous = abs(term)
            series = series + term
        end do
        value = sqrt((0.5_dp * scifort_pi) / z) * exp(-z) * series
    end function complex_k_asymptotic

    pure subroutine ber_zeros(nt, zeros, status)
        integer, intent(in) :: nt !! number of positive ber zeros requested, 1..kelvin_max_zeros
        real(dp), allocatable, intent(out) :: zeros(:) !! first nt positive zeros of ber
        integer, intent(out), optional :: status !! status: 0 success; nonzero otherwise

        call compute_kelvin_zeros(nt, root_ber, zeros, status)
    end subroutine ber_zeros

    pure subroutine bei_zeros(nt, zeros, status)
        integer, intent(in) :: nt !! number of positive bei zeros requested, 1..kelvin_max_zeros
        real(dp), allocatable, intent(out) :: zeros(:) !! first nt positive zeros of bei
        integer, intent(out), optional :: status !! status: 0 success; nonzero otherwise

        call compute_kelvin_zeros(nt, root_bei, zeros, status)
    end subroutine bei_zeros

    pure subroutine ker_zeros(nt, zeros, status)
        integer, intent(in) :: nt !! number of positive ker zeros requested, 1..kelvin_max_zeros
        real(dp), allocatable, intent(out) :: zeros(:) !! first nt positive zeros of ker
        integer, intent(out), optional :: status !! status: 0 success; nonzero otherwise

        call compute_kelvin_zeros(nt, root_ker, zeros, status)
    end subroutine ker_zeros

    pure subroutine kei_zeros(nt, zeros, status)
        integer, intent(in) :: nt !! number of positive kei zeros requested, 1..kelvin_max_zeros
        real(dp), allocatable, intent(out) :: zeros(:) !! first nt positive zeros of kei
        integer, intent(out), optional :: status !! status: 0 success; nonzero otherwise

        call compute_kelvin_zeros(nt, root_kei, zeros, status)
    end subroutine kei_zeros

    pure subroutine berp_zeros(nt, zeros, status)
        integer, intent(in) :: nt !! positive ber-prime zeros requested, 1..kelvin_max_zeros
        real(dp), allocatable, intent(out) :: zeros(:) !! first nt positive zeros of ber derivative
        integer, intent(out), optional :: status !! status: 0 success; nonzero otherwise

        call compute_kelvin_zeros(nt, root_berp, zeros, status)
    end subroutine berp_zeros

    pure subroutine beip_zeros(nt, zeros, status)
        integer, intent(in) :: nt !! positive bei-prime zeros requested, 1..kelvin_max_zeros
        real(dp), allocatable, intent(out) :: zeros(:) !! first nt positive zeros of bei derivative
        integer, intent(out), optional :: status !! status: 0 success; nonzero otherwise

        call compute_kelvin_zeros(nt, root_beip, zeros, status)
    end subroutine beip_zeros

    pure subroutine kerp_zeros(nt, zeros, status)
        integer, intent(in) :: nt !! positive ker-prime zeros requested, 1..kelvin_max_zeros
        real(dp), allocatable, intent(out) :: zeros(:) !! first nt positive zeros of ker derivative
        integer, intent(out), optional :: status !! status: 0 success; nonzero otherwise

        call compute_kelvin_zeros(nt, root_kerp, zeros, status)
    end subroutine kerp_zeros

    pure subroutine keip_zeros(nt, zeros, status)
        integer, intent(in) :: nt !! positive kei-prime zeros requested, 1..kelvin_max_zeros
        real(dp), allocatable, intent(out) :: zeros(:) !! first nt positive zeros of kei derivative
        integer, intent(out), optional :: status !! status: 0 success; nonzero otherwise

        call compute_kelvin_zeros(nt, root_keip, zeros, status)
    end subroutine keip_zeros

    pure subroutine kelvin_zeros(nt, berz, beiz, kerz, keiz, berpz, beipz, kerpz, keipz, status)
        integer, intent(in) :: nt !! zeros requested from each Kelvin family, 1..kelvin_max_zeros
        real(dp), allocatable, intent(out) :: berz(:) !! positive zeros of ber
        real(dp), allocatable, intent(out) :: beiz(:) !! positive zeros of bei
        real(dp), allocatable, intent(out) :: kerz(:) !! positive zeros of ker
        real(dp), allocatable, intent(out) :: keiz(:) !! positive zeros of kei
        real(dp), allocatable, intent(out) :: berpz(:) !! positive zeros of ber derivative
        real(dp), allocatable, intent(out) :: beipz(:) !! positive zeros of bei derivative
        real(dp), allocatable, intent(out) :: kerpz(:) !! positive zeros of ker derivative
        real(dp), allocatable, intent(out) :: keipz(:) !! positive zeros of kei derivative
        integer, intent(out), optional :: status !! status: 0 success; nonzero otherwise

        integer :: info

        call compute_kelvin_zeros(nt, root_ber, berz, info)
        if (info /= kelvin_status_success) then
            call allocate_empty_remaining(beiz, kerz, keiz, berpz, beipz, kerpz, keipz)
            if (present(status)) status = info
            return
        end if
        call compute_kelvin_zeros(nt, root_bei, beiz, info)
        if (info /= kelvin_status_success) then
            call allocate_empty_remaining(kerz, keiz, berpz, beipz, kerpz, keipz)
            if (present(status)) status = info
            return
        end if
        call compute_kelvin_zeros(nt, root_ker, kerz, info)
        if (info /= kelvin_status_success) then
            call allocate_empty_remaining(keiz, berpz, beipz, kerpz, keipz)
            if (present(status)) status = info
            return
        end if
        call compute_kelvin_zeros(nt, root_kei, keiz, info)
        if (info /= kelvin_status_success) then
            call allocate_empty_remaining(berpz, beipz, kerpz, keipz)
            if (present(status)) status = info
            return
        end if
        call compute_kelvin_zeros(nt, root_berp, berpz, info)
        if (info /= kelvin_status_success) then
            call allocate_empty_remaining(beipz, kerpz, keipz)
            if (present(status)) status = info
            return
        end if
        call compute_kelvin_zeros(nt, root_beip, beipz, info)
        if (info /= kelvin_status_success) then
            call allocate_empty_remaining(kerpz, keipz)
            if (present(status)) status = info
            return
        end if
        call compute_kelvin_zeros(nt, root_kerp, kerpz, info)
        if (info /= kelvin_status_success) then
            call allocate_empty_remaining(keipz)
            if (present(status)) status = info
            return
        end if
        call compute_kelvin_zeros(nt, root_keip, keipz, info)
        if (present(status)) status = info
    end subroutine kelvin_zeros

    pure subroutine allocate_empty_remaining(a1, a2, a3, a4, a5, a6, a7)
        real(dp), allocatable, intent(out) :: a1(:) !! first unused output array
        real(dp), allocatable, intent(out), optional :: a2(:) !! second unused output array
        real(dp), allocatable, intent(out), optional :: a3(:) !! third unused output array
        real(dp), allocatable, intent(out), optional :: a4(:) !! fourth unused output array
        real(dp), allocatable, intent(out), optional :: a5(:) !! fifth unused output array
        real(dp), allocatable, intent(out), optional :: a6(:) !! sixth unused output array
        real(dp), allocatable, intent(out), optional :: a7(:) !! seventh unused output array

        allocate(a1(0))
        if (present(a2)) allocate(a2(0))
        if (present(a3)) allocate(a3(0))
        if (present(a4)) allocate(a4(0))
        if (present(a5)) allocate(a5(0))
        if (present(a6)) allocate(a6(0))
        if (present(a7)) allocate(a7(0))
    end subroutine allocate_empty_remaining

    pure subroutine compute_kelvin_zeros(nt, root_kind, zeros, status)
        integer, intent(in) :: nt !! number of positive roots requested, 1..kelvin_max_zeros
        integer, intent(in) :: root_kind !! internal Kelvin root selector
        real(dp), allocatable, intent(out) :: zeros(:) !! first nt positive roots
        integer, intent(out), optional :: status !! status: 0 success; nonzero otherwise

        integer :: count
        integer :: iteration
        integer :: info
        logical :: root_found
        real(dp) :: a
        real(dp) :: b
        real(dp) :: fa
        real(dp) :: fb
        real(dp) :: fm
        real(dp) :: mid
        real(dp) :: x
        real(dp) :: x_next

        if (nt < 1 .or. nt > kelvin_max_zeros .or. root_kind < root_ber .or. &
                root_kind > root_keip) then
            allocate(zeros(0))
            if (present(status)) status = kelvin_status_invalid_argument
            return
        end if

        allocate(zeros(nt))
        zeros = quiet_nan(0.0_dp)
        count = 0
        x = 0.125_dp
        fa = kelvin_root_value(root_kind, x)
        if (.not. ieee_is_finite(fa)) then
            if (present(status)) status = kelvin_status_numerical_failure
            return
        end if

        do while (count < nt)
            x_next = x + 0.5_dp
            fb = kelvin_root_value(root_kind, x_next)
            if (.not. ieee_is_finite(fb)) then
                if (present(status)) status = kelvin_status_numerical_failure
                return
            end if
            if (fa == 0.0_dp .or. fb == 0.0_dp .or. &
                    (fa < 0.0_dp .and. fb > 0.0_dp) .or. &
                    (fa > 0.0_dp .and. fb < 0.0_dp)) then
                a = x
                b = x_next
                root_found = .false.
                do iteration = 1, 80
                    mid = 0.5_dp * (a + b)
                    fm = kelvin_root_value(root_kind, mid)
                    if (.not. ieee_is_finite(fm)) exit
                    if (fm == 0.0_dp .or. b - a <= &
                            16.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(mid))) then
                        count = count + 1
                        zeros(count) = mid
                        root_found = .true.
                        exit
                    end if
                    if ((fa < 0.0_dp .and. fm > 0.0_dp) .or. &
                            (fa > 0.0_dp .and. fm < 0.0_dp)) then
                        b = mid
                        fb = fm
                    else
                        a = mid
                        fa = fm
                    end if
                end do
                if (.not. root_found) then
                    if (present(status)) status = kelvin_status_numerical_failure
                    return
                end if
                x = x_next
                fa = fb
            else
                x = x_next
                fa = fb
            end if
        end do

        info = kelvin_status_success
        if (present(status)) status = info
    end subroutine compute_kelvin_zeros

    pure function kelvin_root_value(root_kind, x) result(value)
        integer, intent(in) :: root_kind !! internal Kelvin root selector
        real(dp), intent(in) :: x !! positive finite root-search argument
        real(dp) :: value

        select case (root_kind)
        case (root_ber)
            value = ber(x)
        case (root_bei)
            value = bei(x)
        case (root_ker)
            value = ker(x)
        case (root_kei)
            value = kei(x)
        case (root_berp)
            value = berp(x)
        case (root_beip)
            value = beip(x)
        case (root_kerp)
            value = kerp(x)
        case (root_keip)
            value = keip(x)
        case default
            value = quiet_nan(x)
        end select
    end function kelvin_root_value

end module scifort_kelvin
