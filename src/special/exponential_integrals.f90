! SPDX-License-Identifier: BSD-3-Clause
! Copyright (c) 1996 California Institute of Technology, Pasadena, California.
! ALL RIGHTS RESERVED. Based on Government Sponsored Research NAS7-03001.
!
! The sici implementation adapts the MATH77 DSI, DCI, and DCPVAL routines by
! W. V. Snyder, C. L. Lawson, F. T. Krogh, and other MATH77 contributors.
! The exp1, expi, expn, and shichi implementations are original SciFort code.
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
! Significant SciFort changes to the MATH77 portion: free-form Fortran 2018,
! pure elemental tuple-style wrapper, IEEE return semantics, abs(x) continuation
! for the real finite cosine-integral branch, and local Chebyshev evaluation.

module scifort_exponential_integrals
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_constants, only : scifort_pi
    use scifort_kinds, only : dp
    use scifort_math, only : negative_infinity, positive_infinity, quiet_nan
    implicit none
    private

    real(dp), parameter :: euler_gamma = &
        0.577215664901532860606512090082402431_dp
    integer, parameter :: max_iterations = 10000

    real(dp), parameter :: si_coefficients(23) = [ &
        0.5_dp, 0.5_dp, &
        4.052926477680623_dp, -4.063980844911986_dp, &
        2.778756381742663_dp, -1.926565091150656_dp, &
        1.389308771171888_dp, -0.968322236987086_dp, &
        0.530148847916522_dp, -0.211263780976555_dp, &
        0.062033679432003_dp, -0.013867445589417_dp, &
        0.002436221404749_dp, -0.000345469155569_dp, &
        0.000040420271419_dp, -0.000003972908746_dp, &
        0.000000332988589_dp, -0.000000024100076_dp, &
        0.000000001522370_dp, -0.000000000084710_dp, &
        0.000000000004185_dp, -0.000000000000185_dp, &
        0.000000000000007_dp ]

    real(dp), parameter :: ci_coefficients(23) = [ &
        0.5_dp, 0.5_dp, &
        14.992589367813409_dp, -19.386124096607770_dp, &
        12.741870869758071_dp, -8.107903970562531_dp, &
        4.862022348500627_dp, -2.497505088539025_dp, &
        1.008660787358110_dp, -0.312080924825428_dp, &
        0.074678255294576_dp, -0.014110865253535_dp, &
        0.002152046752074_dp, -0.000270212331184_dp, &
        0.000028416945498_dp, -0.000002540125611_dp, &
        0.000000195437144_dp, -0.000000013084020_dp, &
        0.000000000769379_dp, -0.000000000040066_dp, &
        0.000000000001861_dp, -0.000000000000078_dp, &
        0.000000000000003_dp ]

    real(dp), parameter :: f_coefficients(13) = [ &
        0.5_dp, 0.5_dp, &
        0.062263729028927_dp, -0.000233756041393_dp, &
        0.000002453755677_dp, -0.000000058670317_dp, &
        0.000000002356196_dp, -0.000000000136096_dp, &
        0.000000000010308_dp, -0.000000000000964_dp, &
        0.000000000000107_dp, -0.000000000000014_dp, &
        0.000000000000002_dp ]

    real(dp), parameter :: g_coefficients(13) = [ &
        0.5_dp, 0.5_dp, &
        0.003862856096703_dp, -0.000042644182622_dp, &
        0.000000724995950_dp, -0.000000023468225_dp, &
        0.000000001169202_dp, -0.000000000079604_dp, &
        0.000000000006875_dp, -0.000000000000717_dp, &
        0.000000000000087_dp, -0.000000000000012_dp, &
        0.000000000000002_dp ]


    public :: exp1
    public :: expi
    public :: expn
    public :: shichi
    public :: sici

contains

    pure elemental function exp1(x) result(y)
        real(dp), intent(in) :: x !! real argument; supported domain x >= 0
        real(dp) :: y

        y = expn(1, x)
    end function exp1

    pure elemental function expi(x) result(y)
        real(dp), intent(in) :: x !! real argument of the principal real exponential integral
        real(dp) :: y

        integer :: k
        real(dp) :: delta
        real(dp) :: log_value
        real(dp) :: previous
        real(dp) :: sum_value
        real(dp) :: term

        if (ieee_is_nan(x)) then
            y = quiet_nan(x)
            return
        end if
        if (.not. ieee_is_finite(x)) then
            if (x > 0.0_dp) then
                y = positive_infinity(x)
            else
                y = -0.0_dp
            end if
            return
        end if
        if (x == 0.0_dp) then
            y = negative_infinity(x)
            return
        end if
        if (x < 0.0_dp) then
            y = -exp1(-x)
            return
        end if

        if (x <= 64.0_dp) then
            term = x
            sum_value = term
            do k = 2, max_iterations
                term = term * x / real(k, dp)
                delta = term / real(k, dp)
                sum_value = sum_value + delta
                if (abs(delta) <= 2.0_dp * epsilon(sum_value) * &
                        max(1.0_dp, abs(sum_value))) exit
            end do
            y = euler_gamma + log(x) + sum_value
            return
        end if

        term = 1.0_dp
        sum_value = 1.0_dp
        previous = term
        do k = 1, max_iterations
            term = term * real(k, dp) / x
            if (term > previous) exit
            sum_value = sum_value + term
            if (term <= 2.0_dp * epsilon(sum_value) * abs(sum_value)) exit
            previous = term
        end do

        log_value = x - log(x) + log(sum_value)
        if (log_value > log(huge(1.0_dp))) then
            y = positive_infinity(x)
        else
            y = exp(log_value)
        end if
    end function expi

    pure elemental function expn(n, x) result(y)
        integer, intent(in) :: n !! nonnegative integer order
        real(dp), intent(in) :: x !! real argument; supported domain x >= 0
        real(dp) :: y

        integer :: i
        integer :: k
        real(dp) :: a
        real(dp) :: b
        real(dp) :: c
        real(dp) :: d
        real(dp) :: delta
        real(dp) :: fact
        real(dp) :: fpmin
        real(dp) :: h
        real(dp) :: harmonic
        real(dp) :: psi

        if (n < 0 .or. ieee_is_nan(x)) then
            y = quiet_nan(x)
            return
        end if
        if (x < 0.0_dp) then
            y = quiet_nan(x)
            return
        end if
        if (.not. ieee_is_finite(x)) then
            y = 0.0_dp
            return
        end if
        if (x == 0.0_dp) then
            if (n <= 1) then
                y = positive_infinity(x)
            else
                y = 1.0_dp / real(n - 1, dp)
            end if
            return
        end if
        if (n == 0) then
            if (x < 1.0_dp / huge(1.0_dp)) then
                y = positive_infinity(x)
            else
                y = exp(-x) / x
            end if
            return
        end if

        if (x <= 1.0_dp .and. n < 50) then
            if (n == 1) then
                y = -log(x) - euler_gamma
            else
                y = 1.0_dp / real(n - 1, dp)
            end if

            fact = 1.0_dp
            harmonic = 0.0_dp
            if (n > 1) then
                do k = 1, n - 1
                    harmonic = harmonic + 1.0_dp / real(k, dp)
                end do
            end if
            psi = -euler_gamma + harmonic

            do i = 1, max_iterations
                fact = fact * (-x) / real(i, dp)
                if (i == n - 1) then
                    delta = fact * (-log(x) + psi)
                else
                    delta = -fact / real(i - (n - 1), dp)
                end if
                y = y + delta
                if (abs(delta) <= 2.0_dp * epsilon(y) * max(1.0_dp, abs(y))) exit
            end do
            return
        end if

        fpmin = tiny(1.0_dp) / epsilon(1.0_dp)
        b = x + real(n, dp)
        c = 1.0_dp / fpmin
        d = 1.0_dp / b
        h = d

        do i = 1, max_iterations
            a = -real(i, dp) * (real(n - 1, dp) + real(i, dp))
            b = b + 2.0_dp
            d = a * d + b
            if (abs(d) < fpmin) d = sign(fpmin, d)
            c = b + a / c
            if (abs(c) < fpmin) c = sign(fpmin, c)
            d = 1.0_dp / d
            delta = c * d
            h = h * delta
            if (abs(delta - 1.0_dp) <= 8.0_dp * epsilon(delta)) exit
        end do
        y = h * exp(-x)
    end function expn

    pure elemental subroutine shichi(x, shi, chi)
        real(dp), intent(in) :: x !! real argument
        real(dp), intent(out) :: shi !! hyperbolic sine integral Shi(x)
        real(dp), intent(out) :: chi !! real hyperbolic cosine integral Chi(x)

        real(dp) :: ax
        real(dp) :: e1_value
        real(dp) :: ei_value

        if (ieee_is_nan(x)) then
            shi = quiet_nan(x)
            chi = quiet_nan(x)
            return
        end if
        if (x == 0.0_dp) then
            shi = x
            chi = negative_infinity(x)
            return
        end if

        ax = abs(x)
        if (.not. ieee_is_finite(ax)) then
            shi = sign(positive_infinity(ax), x)
            chi = positive_infinity(ax)
            return
        end if

        ei_value = expi(ax)
        e1_value = exp1(ax)
        shi = sign(0.5_dp * (ei_value + e1_value), x)
        chi = 0.5_dp * (ei_value - e1_value)
    end subroutine shichi

    pure elemental subroutine sici(x, si, ci)
        real(dp), intent(in) :: x !! real argument
        real(dp), intent(out) :: si !! sine integral Si(x)
        real(dp), intent(out) :: ci !! real cosine-integral branch for finite x

        real(dp) :: ax
        real(dp) :: fz
        real(dp) :: gz
        real(dp) :: z
        real(dp) :: zw

        if (ieee_is_nan(x)) then
            si = quiet_nan(x)
            ci = quiet_nan(x)
            return
        end if
        if (.not. ieee_is_finite(x)) then
            si = sign(0.5_dp * scifort_pi, x)
            if (x > 0.0_dp) then
                ci = 0.0_dp
            else
                ci = quiet_nan(x)
            end if
            return
        end if
        if (x == 0.0_dp) then
            si = x
            ci = negative_infinity(x)
            return
        end if

        ax = abs(x)
        if (ax < 16.0_dp) then
            z = x / 16.0_dp
            zw = z * z
            si = z * chebyshev_value(si_coefficients, 20, zw)

            z = ax / 16.0_dp
            zw = z * z
            ci = log(ax) - zw * chebyshev_value(ci_coefficients, 20, zw) + &
                euler_gamma
            return
        end if

        z = 16.0_dp / ax
        zw = z * z
        fz = z * chebyshev_value(f_coefficients, 10, zw)
        gz = zw * chebyshev_value(g_coefficients, 10, zw)
        si = sign(0.5_dp * scifort_pi - fz * cos(ax) - gz * sin(ax), x)
        ci = fz * sin(ax) - gz * cos(ax)
    end subroutine sici

    pure function chebyshev_value(p, degree, x) result(value)
        real(dp), intent(in) :: p(:) !! transform pair followed by Chebyshev coefficients
        integer, intent(in) :: degree !! degree of the Chebyshev expansion
        real(dp), intent(in) :: x !! unscaled polynomial argument
        real(dp) :: value

        integer :: j
        real(dp) :: s
        real(dp) :: s2
        real(dp) :: w1
        real(dp) :: w2
        real(dp) :: w3

        s = (x - p(1)) / p(2)
        s2 = 2.0_dp * s
        w1 = 0.0_dp
        w2 = 0.0_dp
        do j = degree + 3, 4, -1
            w3 = w2
            w2 = w1
            w1 = s2 * w2 - w3 + p(j)
        end do
        value = s * w1 - w2 + p(3)
    end function chebyshev_value


end module scifort_exponential_integrals
