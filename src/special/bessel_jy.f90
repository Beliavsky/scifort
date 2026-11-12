! SPDX-License-Identifier: BSD-3-Clause
! Copyright (c) 1996 California Institute of Technology, Pasadena, California.
! ALL RIGHTS RESERVED. Based on Government Sponsored Research NAS7-03001.
!
! Adapted in 2026 for SciFort from the MATH77 DBESJN, DBESYN, and DBESPQ
! routines by E. W. Ng, W. V. Snyder, S. Singletary, C. L. Lawson,
! S. Chan, F. T. Krogh, and other MATH77 contributors.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
! Redistributions of source code must retain the above copyright notice, this
! list of conditions and the following disclaimer. Redistributions in binary
! form must reproduce the above copyright notice, this list of conditions and
! the following disclaimer in the documentation and/or other materials provided
! with the distribution. Neither the name of the California Institute of
! Technology (Caltech) nor the names of its contributors may be used to endorse
! or promote products derived from this software without specific prior written
! permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
! ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
! LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
! CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
! SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
! INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
! CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
! ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.
!
! Significant SciFort changes: free-form Fortran 2018 module conversion,
! scalar elemental APIs, intrinsic machine/gamma services, IEEE special-value
! handling, negative-order reflection, and removal of MATH77 diagnostics.

module scifort_bessel_jy
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_constants, only : scifort_pi
    use scifort_kinds, only : dp
    use scifort_math, only : negative_infinity, positive_infinity, quiet_nan
    implicit none
    private

    integer, parameter :: middle_workspace = 128

    public :: jn
    public :: jv
    public :: jve
    public :: yv
    public :: yve

contains

    pure elemental function jn(v, x) result(y)
        real(dp), intent(in) :: v !! real cylindrical-Bessel order
        real(dp), intent(in) :: x !! real argument; negative values require integer order
        real(dp) :: y

        y = jv(v, x)
    end function jn

    pure elemental function jv(v, x) result(y)
        real(dp), intent(in) :: v !! real cylindrical-Bessel order
        real(dp), intent(in) :: x !! real argument; negative values require integer order
        real(dp) :: y

        logical :: integer_v
        real(dp) :: jpos
        real(dp) :: nu
        real(dp) :: ypos

        if (ieee_is_nan(v) .or. ieee_is_nan(x)) then
            y = quiet_nan(x)
            return
        end if
        if (.not. ieee_is_finite(x)) then
            y = quiet_nan(x)
            return
        end if

        integer_v = is_integer_value(v)
        if (x < 0.0_dp .and. .not. integer_v) then
            y = quiet_nan(x)
            return
        end if
        if (.not. ieee_is_finite(v)) then
            if (x >= 0.0_dp) then
                y = 0.0_dp
            else
                y = quiet_nan(x)
            end if
            return
        end if

        if (x == 0.0_dp) then
            y = jv_at_zero(v, x)
            return
        end if

        if (x < 0.0_dp) then
            y = integer_negative_x_jv(v, x)
            return
        end if

        if (v >= 0.0_dp) then
            y = jv_nonnegative(v, x)
            return
        end if

        nu = -v
        if (integer_v) then
            jpos = jv_nonnegative(nu, x)
            if (odd_integer_value(nu)) jpos = -jpos
            y = jpos
            return
        end if

        jpos = jv_nonnegative(nu, x)
        ypos = yv_nonnegative(nu, x)
        y = reflected_jv(jpos, ypos, nu, x)
    end function jv

    pure elemental function jve(v, x) result(y)
        real(dp), intent(in) :: v !! real order of exp(-abs(aimag(z)))*J_v(z)
        real(dp), intent(in) :: x !! real argument; scaling is one on the real axis
        real(dp) :: y

        if (ieee_is_nan(v) .or. ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (.not. ieee_is_finite(v) .or. .not. ieee_is_finite(x)) then
            y = quiet_nan(x)
        else if (x == 0.0_dp .and. v < 0.0_dp .and. .not. is_integer_value(v)) then
            y = quiet_nan(x)
        else
            y = jv(v, x)
        end if
    end function jve

    pure elemental function yv(v, x) result(y)
        real(dp), intent(in) :: v !! real cylindrical-Bessel order
        real(dp), intent(in) :: x !! nonnegative real argument
        real(dp) :: y

        logical :: integer_v
        real(dp) :: jpos
        real(dp) :: nu
        real(dp) :: ypos

        if (ieee_is_nan(v) .or. ieee_is_nan(x)) then
            y = quiet_nan(x)
            return
        end if
        if (x < 0.0_dp .or. .not. ieee_is_finite(x)) then
            y = quiet_nan(x)
            return
        end if
        if (x == 0.0_dp) then
            if (is_integer_value(v) .and. v < 0.0_dp .and. odd_integer_value(v)) then
                y = positive_infinity(x)
            else
                y = negative_infinity(x)
            end if
            return
        end if
        if (.not. ieee_is_finite(v)) then
            y = quiet_nan(x)
            return
        end if

        integer_v = is_integer_value(v)
        if (v >= 0.0_dp) then
            y = yv_nonnegative(v, x)
            return
        end if

        nu = -v
        if (integer_v) then
            ypos = yv_nonnegative(nu, x)
            if (odd_integer_value(nu)) ypos = -ypos
            y = ypos
            return
        end if

        jpos = jv_nonnegative(nu, x)
        ypos = yv_nonnegative(nu, x)
        y = reflected_yv(jpos, ypos, nu, x)
    end function yv

    pure elemental function yve(v, x) result(y)
        real(dp), intent(in) :: v !! real order of exp(-abs(aimag(z)))*Y_v(z)
        real(dp), intent(in) :: x !! positive real argument; scaling is one on the real axis
        real(dp) :: y

        if (ieee_is_nan(v) .or. ieee_is_nan(x)) then
            y = quiet_nan(x)
        else if (.not. ieee_is_finite(v) .or. .not. ieee_is_finite(x)) then
            y = quiet_nan(x)
        else if (x <= 0.0_dp) then
            y = quiet_nan(x)
        else
            y = yv(v, x)
        end if
    end function yve

    pure elemental function jv_at_zero(v, kind_arg) result(y)
        real(dp), intent(in) :: v !! finite real order evaluated at zero
        real(dp), intent(in) :: kind_arg !! real-kind selector for IEEE special values
        real(dp) :: y

        if (v == 0.0_dp) then
            y = 1.0_dp
        else if (v > 0.0_dp) then
            y = 0.0_dp
        else if (is_integer_value(v)) then
            y = 0.0_dp
            if (odd_integer_value(v)) y = -y
        else if (sinpi_value(-v) > 0.0_dp) then
            y = positive_infinity(kind_arg)
        else
            y = negative_infinity(kind_arg)
        end if
    end function jv_at_zero

    pure elemental function integer_negative_x_jv(v, x) result(y)
        real(dp), intent(in) :: v !! finite integer-valued real order
        real(dp), intent(in) :: x !! negative finite argument
        real(dp) :: y

        real(dp) :: nu
        real(dp) :: parity

        nu = abs(v)
        y = jv_nonnegative(nu, -x)
        parity = 1.0_dp
        if (odd_integer_value(nu)) parity = -1.0_dp
        if (v < 0.0_dp .and. odd_integer_value(nu)) parity = -parity
        y = parity * y
    end function integer_negative_x_jv

    pure elemental function reflected_jv(jpos, ypos, nu, kind_arg) result(y)
        real(dp), intent(in) :: jpos !! J_nu(x) for positive nu
        real(dp), intent(in) :: ypos !! Y_nu(x) for positive nu
        real(dp), intent(in) :: nu !! positive reflected order
        real(dp), intent(in) :: kind_arg !! real-kind selector for IEEE special values
        real(dp) :: y

        real(dp) :: c
        real(dp) :: s

        if (ieee_is_nan(jpos) .or. ieee_is_nan(ypos)) then
            y = quiet_nan(kind_arg)
            return
        end if
        c = cospi_value(nu)
        s = sinpi_value(nu)
        if (c == 0.0_dp) then
            y = -s * ypos
        else if (s == 0.0_dp) then
            y = c * jpos
        else if (.not. ieee_is_finite(jpos) .or. .not. ieee_is_finite(ypos)) then
            if (.not. ieee_is_finite(ypos)) then
                y = sign(positive_infinity(kind_arg), -s * ypos)
            else
                y = sign(positive_infinity(kind_arg), c * jpos)
            end if
        else
            y = c * jpos - s * ypos
        end if
    end function reflected_jv

    pure elemental function reflected_yv(jpos, ypos, nu, kind_arg) result(y)
        real(dp), intent(in) :: jpos !! J_nu(x) for positive nu
        real(dp), intent(in) :: ypos !! Y_nu(x) for positive nu
        real(dp), intent(in) :: nu !! positive reflected order
        real(dp), intent(in) :: kind_arg !! real-kind selector for IEEE special values
        real(dp) :: y

        real(dp) :: c
        real(dp) :: s

        if (ieee_is_nan(jpos) .or. ieee_is_nan(ypos)) then
            y = quiet_nan(kind_arg)
            return
        end if
        c = cospi_value(nu)
        s = sinpi_value(nu)
        if (c == 0.0_dp) then
            y = s * jpos
        else if (s == 0.0_dp) then
            y = c * ypos
        else if (.not. ieee_is_finite(jpos) .or. .not. ieee_is_finite(ypos)) then
            if (.not. ieee_is_finite(ypos)) then
                y = sign(positive_infinity(kind_arg), c * ypos)
            else
                y = sign(positive_infinity(kind_arg), s * jpos)
            end if
        else
            y = c * ypos + s * jpos
        end if
    end function reflected_yv

    pure elemental function jv_nonnegative(nu, x) result(y)
        real(dp), intent(in) :: nu !! finite nonnegative real order
        real(dp), intent(in) :: x !! positive finite real argument
        real(dp) :: y

        integer :: nmin
        real(dp) :: max_integer

        intrinsic :: bessel_jn

        max_integer = real(huge(nmin), dp)
        if (is_integer_value(nu) .and. nu <= max_integer) then
            ! Some libm implementations of the Fortran BESSEL_JN intrinsic
            ! overflow internally for large n with x << n even though the
            ! mathematical result is tiny.  In that regime the defining
            ! power series is leading-term dominated and rapidly convergent;
            ! it also detects binary64 underflow from the logarithmic leading
            ! factor before forming any large recurrence intermediates.
            if (x <= sqrt(0.5_dp * (nu + 1.0_dp))) then
                y = jv_power_series(nu, x)
            else
                y = bessel_jn(int(nu), x)
            end if
            return
        end if
        if (nu > max_integer) then
            if (nu > x) then
                y = 0.0_dp
            else
                y = quiet_nan(x)
            end if
            return
        end if

        y = math77_jv_positive(nu, x)
    end function jv_nonnegative

    pure elemental function yv_nonnegative(nu, x) result(y)
        real(dp), intent(in) :: nu !! finite nonnegative real order
        real(dp), intent(in) :: x !! positive finite real argument
        real(dp) :: y

        integer :: nmin
        real(dp) :: max_integer

        intrinsic :: bessel_yn

        max_integer = real(huge(nmin), dp)
        if (is_integer_value(nu) .and. nu <= max_integer) then
            y = bessel_yn(int(nu), x)
            return
        end if
        if (nu > max_integer) then
            if (nu > x) then
                y = negative_infinity(x)
            else
                y = quiet_nan(x)
            end if
            return
        end if

        y = math77_yv_positive(nu, x)
    end function yv_nonnegative

    pure elemental function math77_jv_positive(alpha, x) result(y)
        real(dp), intent(in) :: alpha !! finite nonnegative noninteger order
        real(dp), intent(in) :: x !! positive finite real argument
        real(dp) :: y

        integer :: em_int
        integer :: mu
        integer :: n
        integer :: nmin
        real(dp) :: besv
        real(dp) :: besvm1
        real(dp) :: besvp1
        real(dp) :: big
        real(dp) :: chi
        real(dp) :: dr
        real(dp) :: em
        real(dp) :: emu
        real(dp) :: eps
        real(dp) :: eta
        real(dp) :: fac
        real(dp) :: fk
        real(dp) :: fkm1
        real(dp) :: fkp1
        real(dp) :: fv
        real(dp) :: fvm1
        real(dp) :: fvp1
        real(dp) :: halfpi
        real(dp) :: hicut
        real(dp) :: p
        real(dp) :: q
        real(dp) :: scale
        real(dp) :: small
        real(dp) :: sum_value
        real(dp) :: target
        real(dp) :: test
        real(dp) :: twodx
        real(dp) :: v
        real(dp) :: vpmu
        real(dp) :: xpq
        logical :: target_set

        eps = 0.5_dp * epsilon(1.0_dp)
        halfpi = 0.5_dp * scifort_pi
        hicut = 1.0_dp / (16.0_dp * eps)
        xpq = 1.1293_dp * (0.60206_dp - log10(eps)) - 0.59_dp
        small = 16.0_dp * tiny(1.0_dp) / eps
        big = 0.5_dp * huge(1.0_dp)
        nmin = int(alpha)
        v = alpha - real(nmin, dp)

        if (x <= 0.1_dp) then
            y = jv_power_series(alpha, x)
            return
        end if

        if (x < max(real(nmin + 1, dp), xpq)) then
            twodx = 2.0_dp / x
            mu = max(int(x) + 1, nmin)
            dr = twodx * (v + real(mu, dp))
            fkp1 = 1.0_dp
            fk = 0.0_dp
            do while (eps * abs(fkp1) <= 1.0_dp)
                mu = mu + 1
                dr = dr + twodx
                fkm1 = fk
                fk = fkp1
                fkp1 = dr * fk - fkm1
            end do
            if (mod(mu, 2) /= 0) mu = mu + 1

            fvm1 = small
            fv = 0.0_dp
            eta = 1.0_dp
            sum_value = fvm1
            em_int = mu / 2
            em = real(em_int, dp)
            emu = real(mu, dp)
            fac = (v + emu) * twodx
            test = big / max(1.0_dp, fac)
            target = 0.0_dp
            target_set = .false.

            do
                fvp1 = fv
                fv = fvm1
                if (abs(fv) > test) then
                    fv = fv / sum_value
                    fvp1 = fvp1 / sum_value
                    if (target_set) target = target / sum_value
                    sum_value = 1.0_dp
                end if
                fvm1 = fac * fv - fvp1
                mu = mu - 1
                emu = emu - 1.0_dp
                fac = (v + emu) * twodx
                if (mu == nmin) then
                    target = fvm1
                    target_set = .true.
                end if
                if (mod(mu, 2) == 0) then
                    if (v == 0.0_dp) then
                        sum_value = sum_value + fvm1
                        if (mu == 0) then
                            scale = 1.0_dp / sum_value
                            exit
                        end if
                        sum_value = sum_value + fvm1
                    else if (mu /= 0) then
                        vpmu = v + emu
                        eta = eta * (em / (v + em - 1.0_dp)) * &
                            (vpmu / (vpmu + 2.0_dp))
                        sum_value = sum_value + fvm1 * eta
                        em = em - 1.0_dp
                    else
                        eta = eta / (v + 2.0_dp)
                        sum_value = sum_value + fvm1 * eta
                        scale = gamma(v + 1.0_dp) / eta * sum_value * twodx**v
                        scale = 1.0_dp / scale
                        exit
                    end if
                end if
            end do
            y = target * scale
            return
        end if

        if (x >= hicut) then
            y = quiet_nan(x)
            return
        end if

        call math77_bespq(x, v, p, q)
        chi = x - (v + 0.5_dp) * halfpi
        besv = sqrt(1.0_dp / (halfpi * x)) * (p * cos(chi) - q * sin(chi))
        if (nmin == 0) then
            y = besv
            return
        end if

        call math77_bespq(x, v + 1.0_dp, p, q)
        chi = x - (v + 1.5_dp) * halfpi
        besvp1 = sqrt(1.0_dp / (halfpi * x)) * (p * cos(chi) - q * sin(chi))
        if (nmin == 1) then
            y = besvp1
            return
        end if

        twodx = 2.0_dp / x
        do n = 1, nmin - 1
            besvm1 = besv
            besv = besvp1
            besvp1 = (v + real(n, dp)) * twodx * besv - besvm1
        end do
        y = besvp1
    end function math77_jv_positive

    pure elemental function math77_yv_positive(alpha, x) result(y)
        real(dp), intent(in) :: alpha !! finite positive noninteger order
        real(dp), intent(in) :: x !! positive finite real argument
        real(dp) :: y

        real(dp), parameter :: d2ser(0:10) = [ &
            0.3674669051966159615185_dp, -1.782903980807269842231_dp, &
            0.9411644685122855908427_dp, -1.958865250248747878077_dp, &
            1.557306621108283294475_dp, -2.521051413546812096437_dp, &
            2.184502685635110914145_dp, -3.140067153452674402872_dp, &
            2.817401038921461361158_dp, -3.772198775799670081858_dp, &
            3.452780604492584575060_dp]

        integer :: i
        integer :: m
        integer :: mu
        integer :: n
        integer :: nmin
        real(dp) :: aj(middle_workspace)
        real(dp) :: arg
        real(dp) :: bgam
        real(dp) :: bgamsq
        real(dp) :: big
        real(dp) :: biglog
        real(dp) :: chi
        real(dp) :: d1
        real(dp) :: d2
        real(dp) :: d2val
        real(dp) :: d3
        real(dp) :: em1
        real(dp) :: en1
        real(dp) :: en2
        real(dp) :: en3
        real(dp) :: eps
        real(dp) :: g
        real(dp) :: g0
        real(dp) :: g1
        real(dp) :: gmax
        real(dp) :: halfpi
        real(dp) :: hicut
        real(dp) :: p
        real(dp) :: piv2
        real(dp) :: psi1
        real(dp) :: psiz
        real(dp) :: q
        real(dp) :: q2dxpv
        real(dp) :: test
        real(dp) :: thvdx
        real(dp) :: twodx
        real(dp) :: v
        real(dp) :: v2
        real(dp) :: xpq
        real(dp) :: ylog
        real(dp) :: yv0
        real(dp) :: yv1
        real(dp) :: yvm1
        real(dp) :: z
        logical :: flag
        logical :: j1smal
        logical :: ok

        eps = 0.5_dp * epsilon(1.0_dp)
        halfpi = 0.5_dp * scifort_pi
        hicut = 1.0_dp / (16.0_dp * eps)
        xpq = 1.1293_dp * (0.60206_dp - log10(eps)) - 0.59_dp
        big = 0.5_dp * huge(1.0_dp)
        biglog = log(big)
        nmin = int(alpha)
        v = alpha - real(nmin, dp)

        if (x < eps) then
            ylog = log_gamma(alpha) - alpha * (log(x) - log(2.0_dp)) - log(scifort_pi)
            if (ylog >= biglog) then
                y = negative_infinity(x)
            else
                y = -exp(ylog)
            end if
            return
        end if

        twodx = 2.0_dp / x
        if (x <= xpq) then
            call math77_fractional_j_sequence(x, v, aj, mu, bgam, q2dxpv, ok)
            if (.not. ok) then
                y = quiet_nan(x)
                return
            end if

            v2 = v
            bgamsq = bgam**2
            if (v > 0.5_dp) then
                v2 = v - 1.0_dp
                q2dxpv = q2dxpv / twodx
                bgamsq = bgamsq / v**2
            end if
            piv2 = scifort_pi * v2
            if (abs(v2) < 0.012_dp) then
                arg = 2.0_dp * v2 * log(twodx)
                if (abs(arg) < log(2.0_dp)) then
                    z = tanh(0.5_dp * arg)
                    em1 = 2.0_dp * z / (1.0_dp - z)
                else
                    em1 = exp(arg) - 1.0_dp
                end if
                d2val = d2ser(10)
                do i = 9, 0, -1
                    d2val = d2ser(i) + v2 * d2val
                end do
                g0 = d2val - bgamsq * (em1 / piv2)
                g1 = (q2dxpv**2 / halfpi) * bgamsq * (2.0_dp + v2) / (1.0_dp - v2)
            else
                g0 = 1.0_dp / tan(piv2) - q2dxpv**2 * bgamsq / piv2
                g1 = (q2dxpv**2 / halfpi) * bgamsq * (2.0_dp + v2) / (1.0_dp - v2)
            end if

            en3 = v2 + 1.0_dp
            en2 = v2 + en3
            en1 = v2 + 4.0_dp
            d1 = 2.0_dp
            d2 = d1 - v2
            d3 = d1 + v2
            flag = .false.
            j1smal = abs(aj(1)) < 0.12_dp
            if (j1smal .or. v2 < 0.0_dp) then
                flag = .true.
                thvdx = 3.0_dp * v2 / x
                psiz = -bgamsq * q2dxpv**2 / (halfpi * x)
                psi1 = g0 - 0.5_dp * g1
            end if

            if (v2 >= 0.0_dp) then
                m = 3
                yv0 = g0 * aj(1)
                if (j1smal) yv1 = psiz * aj(1) + psi1 * aj(2)
            else
                z = twodx * v * aj(1) - aj(2)
                yv0 = g0 * z
                m = 2
                yv1 = psiz * z + psi1 * aj(1)
            end if

            do i = m, mu, 2
                yv0 = g1 * aj(i) + yv0
                g = g1
                g1 = -g1 * (en1 / d1) * (en2 / d2) * (en3 / d3)
                en1 = en1 + 2.0_dp
                en2 = en2 + 1.0_dp
                en3 = en3 + 1.0_dp
                d1 = d1 + 1.0_dp
                d2 = d2 + 1.0_dp
                d3 = d3 + 2.0_dp
                if (flag) then
                    yv1 = yv1 + thvdx * g * aj(i) + 0.5_dp * (g - g1) * aj(i + 1)
                end if
            end do

            if (v2 < 0.0_dp) then
                z = yv1
                yv1 = v * z * twodx - yv0
                yv0 = z
            else if (.not. j1smal) then
                yv1 = (yv0 * aj(2) - 1.0_dp / (x * halfpi)) / aj(1)
            end if
        else if (x <= hicut) then
            call math77_bespq(x, v, p, q)
            chi = x - (v + 0.5_dp) * halfpi
            yv0 = sqrt(1.0_dp / (halfpi * x)) * (p * sin(chi) + q * cos(chi))
            call math77_bespq(x, v + 1.0_dp, p, q)
            chi = x - (v + 1.5_dp) * halfpi
            yv1 = sqrt(1.0_dp / (halfpi * x)) * (p * sin(chi) + q * cos(chi))
        else
            y = quiet_nan(x)
            return
        end if

        if (nmin == 0) then
            y = yv0
            return
        else if (nmin == 1) then
            y = yv1
            return
        end if

        g = v * twodx
        gmax = g + twodx * real(nmin - 1, dp)
        test = big / max(1.0_dp, gmax)
        do n = 1, nmin - 1
            yvm1 = yv0
            yv0 = yv1
            if (abs(yv0) > test) then
                y = sign(positive_infinity(x), yv0)
                return
            end if
            g = g + twodx
            yv1 = g * yv0 - yvm1
        end do
        y = yv1
    end function math77_yv_positive

    pure subroutine math77_fractional_j_sequence(x, v, aj, mu_out, bgam, q2dxpv, ok)
        real(dp), intent(in) :: x !! positive argument, no larger than the MATH77 middle cutoff
        real(dp), intent(in) :: v !! fractional order in the open interval (0,1)
        real(dp), intent(out) :: aj(:) !! workspace receiving J_v through J_(v+mu)
        integer, intent(out) :: mu_out !! highest fractional-order offset stored in aj
        real(dp), intent(out) :: bgam !! gamma(v+1) used by the companion Y algorithm
        real(dp), intent(out) :: q2dxpv !! (2/x)**v used by the companion Y algorithm
        logical, intent(out) :: ok !! true when the fixed workspace was sufficient

        integer :: em_int
        integer :: ii
        integer :: mu
        integer :: musave
        real(dp) :: big
        real(dp) :: dr
        real(dp) :: em
        real(dp) :: emu
        real(dp) :: eps
        real(dp) :: eta
        real(dp) :: fac
        real(dp) :: fk
        real(dp) :: fkm1
        real(dp) :: fkp1
        real(dp) :: fv
        real(dp) :: fvm1
        real(dp) :: fvp1
        real(dp) :: scale
        real(dp) :: small
        real(dp) :: sum_value
        real(dp) :: test
        real(dp) :: twodx
        real(dp) :: vpmu

        eps = 0.5_dp * epsilon(1.0_dp)
        small = 16.0_dp * tiny(1.0_dp) / eps
        big = 0.5_dp * huge(1.0_dp)
        twodx = 2.0_dp / x
        mu = int(x) + 1
        dr = twodx * (v + real(mu, dp))
        fkp1 = 1.0_dp
        fk = 0.0_dp
        do while (eps * abs(fkp1) <= 1.0_dp)
            mu = mu + 1
            dr = dr + twodx
            fkm1 = fk
            fk = fkp1
            fkp1 = dr * fk - fkm1
        end do
        if (mod(mu, 2) /= 0) mu = mu + 1
        musave = mu
        if (musave + 1 > size(aj)) then
            aj = 0.0_dp
            mu_out = 0
            bgam = 0.0_dp
            q2dxpv = 0.0_dp
            ok = .false.
            return
        end if

        aj = 0.0_dp
        fvm1 = small
        aj(mu + 1) = fvm1
        fv = 0.0_dp
        eta = 1.0_dp
        sum_value = fvm1
        em_int = mu / 2
        em = real(em_int, dp)
        emu = real(mu, dp)
        fac = (v + emu) * twodx
        test = big / max(1.0_dp, fac)

        do
            fvp1 = fv
            fv = fvm1
            if (abs(fv) > test) then
                fv = fv / sum_value
                fvp1 = fvp1 / sum_value
                do ii = mu + 1, musave
                    aj(ii) = aj(ii) / sum_value
                end do
                sum_value = 1.0_dp
            end if
            fvm1 = fac * fv - fvp1
            mu = mu - 1
            emu = emu - 1.0_dp
            fac = (v + emu) * twodx
            aj(mu + 1) = fvm1
            if (mod(mu, 2) == 0) then
                if (mu /= 0) then
                    vpmu = v + emu
                    eta = eta * (em / (v + em - 1.0_dp)) * &
                        (vpmu / (vpmu + 2.0_dp))
                    sum_value = sum_value + fvm1 * eta
                    em = em - 1.0_dp
                else
                    eta = eta / (v + 2.0_dp)
                    sum_value = sum_value + fvm1 * eta
                    bgam = gamma(v + 1.0_dp)
                    q2dxpv = twodx**v
                    scale = 1.0_dp / ((bgam / eta) * sum_value * q2dxpv)
                    exit
                end if
            end if
        end do

        aj(1:musave + 1) = aj(1:musave + 1) * scale
        mu_out = musave
        ok = .true.
    end subroutine math77_fractional_j_sequence

    pure elemental function jv_power_series(nu, x) result(y)
        real(dp), intent(in) :: nu !! finite nonnegative order
        real(dp), intent(in) :: x !! positive argument no larger than 0.1
        real(dp) :: y

        integer :: k
        real(dp) :: log_leading
        real(dp) :: sum_value
        real(dp) :: term
        real(dp) :: xx

        log_leading = nu * log(0.5_dp * x) - log_gamma(nu + 1.0_dp)
        if (log_leading < log(tiny(1.0_dp)) + log(epsilon(1.0_dp))) then
            y = 0.0_dp
            return
        end if

        xx = -0.25_dp * x * x
        term = 1.0_dp
        sum_value = 1.0_dp
        do k = 1, 1000
            term = term * xx / (real(k, dp) * (nu + real(k, dp)))
            sum_value = sum_value + term
            if (abs(term) <= epsilon(1.0_dp) * abs(sum_value)) exit
        end do
        y = exp(log_leading) * sum_value
    end function jv_power_series

    pure subroutine math77_bespq(x, v, p, q)
        real(dp), intent(in) :: x !! large positive argument
        real(dp), intent(in) :: v !! order in the closed interval [0,2]
        real(dp), intent(out) :: p !! asymptotic P auxiliary function
        real(dp), intent(out) :: q !! asymptotic Q auxiliary function

        real(dp) :: a0
        real(dp) :: a1
        real(dp) :: b
        real(dp) :: emu
        real(dp) :: psum
        real(dp) :: q1
        real(dp) :: qsum
        real(dp) :: sig
        real(dp) :: small
        real(dp) :: term
        real(dp) :: x8

        small = 0.25_dp * epsilon(1.0_dp)
        emu = 4.0_dp * v * v
        x8 = 8.0_dp * x
        a0 = 9.0_dp
        a1 = 8.0_dp
        b = 2.0_dp
        sig = -1.0_dp
        psum = 0.0_dp
        term = (emu - 1.0_dp) / x8
        q1 = term
        qsum = 0.0_dp

        do
            term = term * (emu - a0) / (b * x8)
            if (abs(term) <= small) exit
            psum = psum + sig * term
            a1 = a1 + 8.0_dp
            a0 = a0 + a1
            b = b + 1.0_dp

            term = term * (emu - a0) / (b * x8)
            if (abs(term) <= small) exit
            qsum = qsum + sig * term
            a1 = a1 + 8.0_dp
            a0 = a0 + a1
            b = b + 1.0_dp
            sig = -sig
        end do

        p = 1.0_dp + psum
        q = q1 + qsum
    end subroutine math77_bespq

    pure elemental logical function is_integer_value(v) result(is_integer)
        real(dp), intent(in) :: v !! finite real order to classify

        is_integer = v == anint(v)
    end function is_integer_value

    pure elemental logical function odd_integer_value(v) result(is_odd)
        real(dp), intent(in) :: v !! finite integer-valued real order

        is_odd = modulo(abs(v), 2.0_dp) == 1.0_dp
    end function odd_integer_value

    pure elemental function sinpi_value(v) result(y)
        real(dp), intent(in) :: v !! finite real multiplier of pi
        real(dp) :: y
        real(dp) :: r

        r = modulo(v, 2.0_dp)
        if (r == 0.0_dp .or. r == 1.0_dp) then
            y = 0.0_dp
        else if (r == 0.5_dp) then
            y = 1.0_dp
        else if (r == 1.5_dp) then
            y = -1.0_dp
        else if (r < 0.5_dp) then
            y = sin(scifort_pi * r)
        else if (r < 1.5_dp) then
            y = sin(scifort_pi * (1.0_dp - r))
        else
            y = sin(scifort_pi * (r - 2.0_dp))
        end if
    end function sinpi_value

    pure elemental function cospi_value(v) result(y)
        real(dp), intent(in) :: v !! finite real multiplier of pi
        real(dp) :: y
        real(dp) :: r

        r = modulo(v, 2.0_dp)
        if (r == 0.0_dp) then
            y = 1.0_dp
        else if (r == 0.5_dp .or. r == 1.5_dp) then
            y = 0.0_dp
        else if (r == 1.0_dp) then
            y = -1.0_dp
        else if (r < 1.0_dp) then
            y = cos(scifort_pi * r)
        else
            y = cos(scifort_pi * (2.0_dp - r))
        end if
    end function cospi_value

end module scifort_bessel_jy
