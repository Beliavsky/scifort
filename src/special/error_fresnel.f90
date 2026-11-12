! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Real Dawson, imaginary-error, scaled-complementary-error, and Fresnel
! functions. The implementation is independent and follows the defining
! integrals and asymptotic expansions in NIST DLMF Chapter 7.

module scifort_error_fresnel
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_constants, only : scifort_log_two, scifort_pi
    use scifort_kinds, only : dp
    use scifort_math, only : positive_infinity
    implicit none
    private

    real(dp), parameter :: gl16_nodes(8) = [ &
        0.09501250983763744051_dp, &
        0.2816035507792589154_dp, &
        0.4580167776572273697_dp, &
        0.6178762444026437706_dp, &
        0.7554044083550029987_dp, &
        0.8656312023878317552_dp, &
        0.9445750230732326003_dp, &
        0.9894009349916499385_dp]
    real(dp), parameter :: gl16_weights(8) = [ &
        0.1894506104550686409_dp, &
        0.1826034150449236393_dp, &
        0.1691565193950026469_dp, &
        0.1495959888165767082_dp, &
        0.1246289712555340712_dp, &
        0.09515851168249260528_dp, &
        0.06225352393864745648_dp, &
        0.02715245941175417621_dp]

    public :: dawsn
    public :: erfcx
    public :: erfi
    public :: fresnel

contains

    pure elemental function dawsn(x) result(y)
        real(dp), intent(in) :: x !! real argument of Dawson's integral
        real(dp) :: y

        integer :: k
        integer :: n
        integer :: nseg
        real(dp) :: a
        real(dp) :: ax
        real(dp) :: b
        real(dp) :: half
        real(dp) :: mid
        real(dp) :: new_term
        real(dp) :: sumv
        real(dp) :: term
        real(dp) :: u1
        real(dp) :: u2
        real(dp) :: x2

        if (ieee_is_nan(x)) then
            y = x
            return
        end if
        if (.not. ieee_is_finite(x)) then
            y = sign(0.0_dp, x)
            return
        end if

        ax = abs(x)
        if (ax == 0.0_dp) then
            y = x
            return
        end if

        x2 = ax * ax
        if (ax <= 0.5_dp) then
            term = ax
            sumv = term
            do n = 1, 100
                term = term * (-2.0_dp * x2) / real(2 * n + 1, dp)
                sumv = sumv + term
                if (abs(term) <= epsilon(sumv) * abs(sumv)) exit
            end do
            y = sign(sumv, x)
            return
        end if

        if (ax <= 6.0_dp) then
            nseg = max(1, ceiling(x2 / 8.0_dp))
            sumv = 0.0_dp
            do n = 0, nseg - 1
                a = real(n, dp) / real(nseg, dp)
                b = real(n + 1, dp) / real(nseg, dp)
                mid = 0.5_dp * (a + b)
                half = 0.5_dp * (b - a)
                do k = 1, size(gl16_nodes)
                    u1 = mid - half * gl16_nodes(k)
                    u2 = mid + half * gl16_nodes(k)
                    sumv = sumv + half * gl16_weights(k) * ( &
                        exp(-x2 * (1.0_dp - u1 * u1)) + &
                        exp(-x2 * (1.0_dp - u2 * u2)))
                end do
            end do
            y = sign(ax * sumv, x)
            return
        end if

        term = 0.5_dp / ax
        sumv = term
        do n = 1, 100
            new_term = term * real(2 * n - 1, dp) / (2.0_dp * x2)
            if (new_term >= term) exit
            sumv = sumv + new_term
            if (new_term <= epsilon(sumv) * abs(sumv)) exit
            term = new_term
        end do
        y = sign(sumv, x)
    end function dawsn

    pure elemental function erfi(x) result(y)
        real(dp), intent(in) :: x !! real argument of the imaginary error function
        real(dp) :: y

        real(dp), parameter :: two_over_sqrt_pi = 2.0_dp / sqrt(scifort_pi)
        real(dp), parameter :: log_huge = log(huge(1.0_dp))
        real(dp) :: ax
        real(dp) :: log_abs

        if (ieee_is_nan(x)) then
            y = x
            return
        end if
        if (x == 0.0_dp) then
            y = x
            return
        end if
        if (.not. ieee_is_finite(x)) then
            y = sign(positive_infinity(x), x)
            return
        end if

        ax = abs(x)
        if (ax > 27.0_dp) then
            y = sign(positive_infinity(x), x)
            return
        end if

        log_abs = ax * ax + log(two_over_sqrt_pi * abs(dawsn(ax)))
        if (log_abs > log_huge) then
            y = sign(positive_infinity(x), x)
        else
            y = sign(exp(log_abs), x)
        end if
    end function erfi

    pure elemental function erfcx(x) result(y)
        real(dp), intent(in) :: x !! real argument of exp(x**2)*erfc(x)
        real(dp) :: y

        real(dp), parameter :: log_huge = log(huge(1.0_dp))
        real(dp) :: ax
        real(dp) :: log_twice_exp

        if (ieee_is_nan(x)) then
            y = x
            return
        end if
        if (.not. ieee_is_finite(x)) then
            if (x > 0.0_dp) then
                y = 0.0_dp
            else
                y = positive_infinity(x)
            end if
            return
        end if

        if (x >= 0.0_dp) then
            y = erfcx_positive(x)
            return
        end if

        ax = -x
        if (ax > 27.0_dp) then
            y = positive_infinity(x)
            return
        end if
        log_twice_exp = ax * ax + scifort_log_two
        if (log_twice_exp > log_huge) then
            y = positive_infinity(x)
        else
            y = exp(log_twice_exp) - erfcx_positive(ax)
        end if
    end function erfcx

    pure elemental subroutine fresnel(x, s, c)
        real(dp), intent(in) :: x !! real argument of the Fresnel integrals
        real(dp), intent(out) :: s !! Fresnel S integral from zero to x
        real(dp), intent(out) :: c !! Fresnel C integral from zero to x

        integer :: k
        integer :: n
        integer :: nseg
        real(dp) :: a
        real(dp) :: ax
        real(dp) :: b
        real(dp) :: f_aux
        real(dp) :: g_aux
        real(dp) :: half
        real(dp) :: mid
        real(dp) :: phase1
        real(dp) :: phase2
        real(dp) :: sin_phase
        real(dp) :: cos_phase
        real(dp) :: tc
        real(dp) :: ts
        real(dp) :: t1
        real(dp) :: t2
        real(dp) :: x4

        if (ieee_is_nan(x)) then
            s = x
            c = x
            return
        end if
        if (.not. ieee_is_finite(x)) then
            s = sign(0.5_dp, x)
            c = sign(0.5_dp, x)
            return
        end if

        ax = abs(x)
        if (ax == 0.0_dp) then
            s = 0.0_dp
            c = 0.0_dp
            return
        end if

        if (ax <= 1.5_dp) then
            x4 = ax**4
            tc = ax
            c = tc
            ts = 0.5_dp * scifort_pi * ax**3 / 3.0_dp
            s = ts
            do n = 0, 100
                tc = tc * (-(0.5_dp * scifort_pi)**2 * x4) * &
                    real(4 * n + 1, dp) / &
                    real((2 * n + 1) * (2 * n + 2) * (4 * n + 5), dp)
                c = c + tc
                ts = ts * (-(0.5_dp * scifort_pi)**2 * x4) * &
                    real(4 * n + 3, dp) / &
                    real((2 * n + 2) * (2 * n + 3) * (4 * n + 7), dp)
                s = s + ts
                if (max(abs(tc), abs(ts)) <= &
                        epsilon(s) * max(abs(s), abs(c))) exit
            end do
        else if (ax < 4.5_dp) then
            nseg = max(1, ceiling(ax))
            s = 0.0_dp
            c = 0.0_dp
            do n = 0, nseg - 1
                a = ax * real(n, dp) / real(nseg, dp)
                b = ax * real(n + 1, dp) / real(nseg, dp)
                mid = 0.5_dp * (a + b)
                half = 0.5_dp * (b - a)
                do k = 1, size(gl16_nodes)
                    t1 = mid - half * gl16_nodes(k)
                    t2 = mid + half * gl16_nodes(k)
                    phase1 = 0.5_dp * scifort_pi * t1 * t1
                    phase2 = 0.5_dp * scifort_pi * t2 * t2
                    s = s + half * gl16_weights(k) * &
                        (sin(phase1) + sin(phase2))
                    c = c + half * gl16_weights(k) * &
                        (cos(phase1) + cos(phase2))
                end do
            end do
        else
            call fresnel_auxiliary(ax, f_aux, g_aux)
            call fresnel_phase(ax, sin_phase, cos_phase)
            s = 0.5_dp - f_aux * cos_phase - g_aux * sin_phase
            c = 0.5_dp + f_aux * sin_phase - g_aux * cos_phase
        end if

        if (x < 0.0_dp) then
            s = -s
            c = -c
        end if
    end subroutine fresnel

    pure function erfcx_positive(x) result(y)
        real(dp), intent(in) :: x !! finite nonnegative erfcx argument
        real(dp) :: y

        integer :: n
        real(dp) :: new_term
        real(dp) :: sumv
        real(dp) :: term
        real(dp) :: x2

        if (x <= 6.0_dp) then
            y = exp(x * x) * erfc(x)
            return
        end if

        x2 = x * x
        term = 1.0_dp
        sumv = term
        do n = 1, 100
            new_term = -term * real(2 * n - 1, dp) / (2.0_dp * x2)
            if (abs(new_term) >= abs(term)) exit
            sumv = sumv + new_term
            if (abs(new_term) <= epsilon(sumv) * abs(sumv)) exit
            term = new_term
        end do
        y = sumv / (sqrt(scifort_pi) * x)
    end function erfcx_positive

    pure subroutine fresnel_auxiliary(x, f_aux, g_aux)
        real(dp), intent(in) :: x !! positive real argument, at least 4.5
        real(dp), intent(out) :: f_aux !! DLMF Fresnel auxiliary f(x)
        real(dp), intent(out) :: g_aux !! DLMF Fresnel auxiliary g(x)

        integer :: m
        real(dp) :: new_term
        real(dp) :: q
        real(dp) :: q2
        real(dp) :: sumv
        real(dp) :: term

        q = 0.5_dp * scifort_pi * x * x
        q2 = q * q

        term = 1.0_dp
        sumv = term
        do m = 0, 100
            new_term = -term * real((4 * m + 1) * (4 * m + 3), dp) / &
                (4.0_dp * q2)
            if (abs(new_term) >= abs(term)) exit
            sumv = sumv + new_term
            if (abs(new_term) <= epsilon(sumv) * abs(sumv)) exit
            term = new_term
        end do
        f_aux = sumv / (scifort_pi * x)

        term = 0.5_dp / q
        sumv = term
        do m = 0, 100
            new_term = -term * real((4 * m + 3) * (4 * m + 5), dp) / &
                (4.0_dp * q2)
            if (abs(new_term) >= abs(term)) exit
            sumv = sumv + new_term
            if (abs(new_term) <= epsilon(sumv) * abs(sumv)) exit
            term = new_term
        end do
        g_aux = sumv / (scifort_pi * x)
    end subroutine fresnel_auxiliary

    pure subroutine fresnel_phase(x, sin_phase, cos_phase)
        real(dp), intent(in) :: x !! positive finite Fresnel argument, at least 4.5
        real(dp), intent(out) :: sin_phase !! sin(pi*x**2/2) after phase reduction
        real(dp), intent(out) :: cos_phase !! cos(pi*x**2/2) after phase reduction

        real(dp), parameter :: splitter = 134217729.0_dp
        real(dp), parameter :: two52 = 4503599627370496.0_dp
        real(dp), parameter :: two53 = 9007199254740992.0_dp
        real(dp) :: err
        real(dp) :: hi
        real(dp) :: lo
        real(dp) :: p
        real(dp) :: phase
        real(dp) :: rem4
        real(dp) :: split

        if (x >= two53) then
            rem4 = 0.0_dp
        else if (x >= two52) then
            if (modulo(x, 2.0_dp) == 0.0_dp) then
                rem4 = 0.0_dp
            else
                rem4 = 1.0_dp
            end if
        else
            split = splitter * x
            hi = split - (split - x)
            lo = x - hi
            p = x * x
            err = ((hi * hi - p) + 2.0_dp * hi * lo) + lo * lo
            rem4 = modulo(modulo(p, 4.0_dp) + modulo(err, 4.0_dp), 4.0_dp)
        end if

        phase = 0.5_dp * scifort_pi * rem4
        sin_phase = sin(phase)
        cos_phase = cos(phase)
    end subroutine fresnel_phase

end module scifort_error_fresnel
