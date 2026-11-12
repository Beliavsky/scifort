! SPDX-License-Identifier: BSD-3-Clause
! Copyright (c) 2016-2022, Jacob Williams
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!
! * Redistributions of source code must retain the above copyright notice,
!   this list of conditions and the following disclaimer.
! * Redistributions in binary form must reproduce the above copyright notice,
!   this list of conditions and the following disclaimer in the documentation
!   and/or other materials provided with the distribution.
! * The names of its contributors may not be used to endorse or promote
!   products derived from this software without specific prior written
!   permission.
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
! Adapted in 2026 for SciFort from Jacob Williams' carlson-elliptic-integrals
! release 1.0.0, a modern Fortran refactoring of the public-domain SLATEC DRC,
! DRD, DRF, and DRJ routines. Significant SciFort changes: fixed binary64 kind, pure
! elemental SciPy-style APIs, IEEE special-value/domain handling, homogeneity
! scaling, Cauchy-principal-value continuations for real RC/RJ, and removal of
! diagnostics and error-stop behavior. Complex extensions are independent
! SciFort implementations from Carlson/DLMF definitions and identities.

module scifort_carlson_elliptic
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : positive_infinity, quiet_nan
    implicit none
    private

    integer, parameter :: max_iterations = 100
    real(dp), parameter :: unit_roundoff = &
        real(radix(1.0_dp), dp)**(-digits(1.0_dp))

    public :: elliprc
    public :: elliprd
    public :: elliprf
    public :: elliprg
    public :: elliprj

    interface elliprc
        module procedure elliprc_real
        module procedure elliprc_complex
    end interface elliprc

    interface elliprd
        module procedure elliprd_real
        module procedure elliprd_complex
    end interface elliprd

    interface elliprf
        module procedure elliprf_real
        module procedure elliprf_complex
    end interface elliprf

    interface elliprg
        module procedure elliprg_real
        module procedure elliprg_complex
    end interface elliprg

    interface elliprj
        module procedure elliprj_real
        module procedure elliprj_complex
    end interface elliprj

contains

    pure elemental function elliprc_real(x, y) result(value)
        real(dp), intent(in) :: x !! first RC argument; must be nonnegative on the real axis
        real(dp), intent(in) :: y !! nonzero second RC argument; negative values use principal value
        real(dp) :: value

        real(dp) :: scale_value
        real(dp) :: xp
        real(dp) :: yp

        if (ieee_is_nan(x) .or. ieee_is_nan(y)) then
            value = quiet_nan(x)
            return
        end if
        if (x < 0.0_dp .or. y == 0.0_dp) then
            value = quiet_nan(x)
            return
        end if
        if (.not. ieee_is_finite(x) .or. .not. ieee_is_finite(y)) then
            value = 0.0_dp
            return
        end if

        if (y < 0.0_dp) then
            yp = -y
            if (x == 0.0_dp) then
                value = 0.0_dp
                return
            end if
            scale_value = max(x, yp)
            xp = x / scale_value
            yp = yp / scale_value
            value = sqrt(xp / (xp + yp)) * elliprc_positive(xp + yp, yp) / &
                sqrt(scale_value)
        else
            value = elliprc_positive(x, y)
        end if

    end function elliprc_real

    pure elemental function elliprf_real(x, y, z) result(value)
        real(dp), intent(in) :: x !! first RF argument; nonnegative
        real(dp), intent(in) :: y !! second RF argument; nonnegative
        real(dp), intent(in) :: z !! third RF argument; nonnegative
        real(dp) :: value

        real(dp) :: scale_value
        real(dp) :: core

        if (ieee_is_nan(x) .or. ieee_is_nan(y) .or. ieee_is_nan(z)) then
            value = quiet_nan(x)
            return
        end if
        if (x < 0.0_dp .or. y < 0.0_dp .or. z < 0.0_dp) then
            value = quiet_nan(x)
            return
        end if
        if (.not. ieee_is_finite(x) .or. .not. ieee_is_finite(y) .or. &
                .not. ieee_is_finite(z)) then
            value = 0.0_dp
            return
        end if
        if ((x == 0.0_dp .and. y == 0.0_dp) .or. &
                (x == 0.0_dp .and. z == 0.0_dp) .or. &
                (y == 0.0_dp .and. z == 0.0_dp)) then
            value = positive_infinity(x)
            return
        end if

        scale_value = max(x, y, z)
        if (scale_value == 0.0_dp) then
            value = positive_infinity(x)
            return
        end if
        core = carlson_rf_core(x / scale_value, y / scale_value, z / scale_value)
        value = core / sqrt(scale_value)
    end function elliprf_real

    pure elemental function elliprd_real(x, y, z) result(value)
        real(dp), intent(in) :: x !! first RD argument; nonnegative
        real(dp), intent(in) :: y !! second RD argument; nonnegative
        real(dp), intent(in) :: z !! positive distinguished RD argument
        real(dp) :: value

        real(dp) :: core
        real(dp) :: scaled_once
        real(dp) :: scale_value

        if (ieee_is_nan(x) .or. ieee_is_nan(y) .or. ieee_is_nan(z)) then
            value = quiet_nan(x)
            return
        end if
        if (x < 0.0_dp .or. y < 0.0_dp .or. z < 0.0_dp) then
            value = quiet_nan(x)
            return
        end if
        if (z == 0.0_dp .or. (x == 0.0_dp .and. y == 0.0_dp)) then
            value = positive_infinity(x)
            return
        end if
        if (.not. ieee_is_finite(x) .or. .not. ieee_is_finite(y) .or. &
                .not. ieee_is_finite(z)) then
            value = 0.0_dp
            return
        end if

        scale_value = max(x, y, z)
        core = carlson_rd_core(x / scale_value, y / scale_value, z / scale_value)
        scaled_once = core / sqrt(scale_value)
        value = divide_by_positive_scale(scaled_once, scale_value, x)
    end function elliprd_real

    pure elemental function elliprj_real(x, y, z, p) result(value)
        real(dp), intent(in) :: x !! first RJ argument; nonnegative
        real(dp), intent(in) :: y !! second RJ argument; nonnegative
        real(dp), intent(in) :: z !! third RJ argument; nonnegative
        real(dp), intent(in) :: p !! nonzero fourth RJ argument; negative values use principal value
        real(dp) :: value

        real(dp) :: core
        real(dp) :: scaled_once
        real(dp) :: scale_value
        real(dp) :: xs
        real(dp) :: ys
        real(dp) :: zs

        if (ieee_is_nan(x) .or. ieee_is_nan(y) .or. ieee_is_nan(z) .or. ieee_is_nan(p)) then
            value = quiet_nan(x)
            return
        end if
        if (x < 0.0_dp .or. y < 0.0_dp .or. z < 0.0_dp .or. p == 0.0_dp) then
            value = quiet_nan(x)
            return
        end if
        if ((x == 0.0_dp .and. y == 0.0_dp) .or. &
                (x == 0.0_dp .and. z == 0.0_dp) .or. &
                (y == 0.0_dp .and. z == 0.0_dp)) then
            value = positive_infinity(x)
            return
        end if
        if (.not. ieee_is_finite(x) .or. .not. ieee_is_finite(y) .or. &
                .not. ieee_is_finite(z) .or. .not. ieee_is_finite(p)) then
            value = 0.0_dp
            return
        end if

        if (p < 0.0_dp) then
            scale_value = max(x, y, z, -p)
            call sort_three(x / scale_value, y / scale_value, z / scale_value, xs, ys, zs)
            core = elliprj_principal_value(xs, ys, zs, p / scale_value)
            scaled_once = core / sqrt(scale_value)
            value = divide_by_positive_scale(scaled_once, scale_value, x)
            return
        end if

        value = elliprj_positive(x, y, z, p)
    end function elliprj_real

    pure elemental function elliprg_real(x, y, z) result(value)
        real(dp), intent(in) :: x !! first RG argument; nonnegative
        real(dp), intent(in) :: y !! second RG argument; nonnegative
        real(dp), intent(in) :: z !! third RG argument; nonnegative
        real(dp) :: value

        real(dp) :: a
        real(dp) :: b
        real(dp) :: c
        real(dp) :: rd_value
        real(dp) :: rf_value
        real(dp) :: scale_value
        real(dp) :: term_root

        if (ieee_is_nan(x) .or. ieee_is_nan(y) .or. ieee_is_nan(z)) then
            value = quiet_nan(x)
            return
        end if
        if (x < 0.0_dp .or. y < 0.0_dp .or. z < 0.0_dp) then
            value = quiet_nan(x)
            return
        end if
        if (.not. ieee_is_finite(x) .or. .not. ieee_is_finite(y) .or. &
                .not. ieee_is_finite(z)) then
            value = positive_infinity(x)
            return
        end if
        if (x == 0.0_dp .and. y == 0.0_dp .and. z == 0.0_dp) then
            value = 0.0_dp
            return
        end if

        scale_value = max(x, y, z)
        call choose_rg_pivot(x / scale_value, y / scale_value, z / scale_value, a, b, c)
        if (a == 0.0_dp .and. b == 0.0_dp) then
            value = 0.5_dp * sqrt(c * scale_value)
            return
        end if

        rf_value = carlson_rf_core(a, b, c)
        rd_value = carlson_rd_core(a, b, c)
        if (a == 0.0_dp .or. b == 0.0_dp) then
            term_root = 0.0_dp
        else
            term_root = sqrt((a * b) / c)
        end if
        value = 0.5_dp * sqrt(scale_value) * &
            (c * rf_value - ((a - c) * (b - c) / 3.0_dp) * rd_value + term_root)
    end function elliprg_real

    pure elemental function elliprc_complex(x, y) result(value)
        complex(dp), intent(in) :: x !! first complex RC argument; off the negative-real cut
        complex(dp), intent(in) :: y !! nonzero second complex RC argument
        complex(dp) :: value

        real(dp) :: scale_value
        real(dp) :: y_abs

        if (complex_has_nan(x) .or. complex_has_nan(y)) then
            value = complex_nan()
            return
        end if
        if (on_negative_real_cut(x) .or. complex_is_zero(y)) then
            value = complex_nan()
            return
        end if
        if (aimag(x) == 0.0_dp .and. aimag(y) == 0.0_dp) then
            value = cmplx(elliprc_real(real(x, dp), real(y, dp)), 0.0_dp, kind=dp)
            return
        end if
        if (complex_has_ambiguous_infinity(x)) then
            value = complex_nan()
            return
        end if
        if (complex_has_infinite(x) .or. complex_has_infinite(y)) then
            value = cmplx(0.0_dp, 0.0_dp, kind=dp)
            return
        end if

        if (aimag(y) == 0.0_dp .and. real(y, dp) < 0.0_dp) then
            ! Carlson's Cauchy-principal-value continuation for a negative
            ! real second argument.  This branch is also needed when x is complex.
            y_abs = -real(y, dp)
            scale_value = max(abs(x), y_abs)
            value = sqrt((x / scale_value) / &
                (x / scale_value + y_abs / scale_value)) * &
                carlson_rc_core_complex(x / scale_value + y_abs / scale_value, &
                cmplx(y_abs / scale_value, 0.0_dp, kind=dp)) / sqrt(scale_value)
            return
        end if

        scale_value = max(abs(x), abs(y))
        value = carlson_rc_core_complex(x / scale_value, y / scale_value) / &
            sqrt(scale_value)
    end function elliprc_complex

    pure elemental function elliprf_complex(x, y, z) result(value)
        complex(dp), intent(in) :: x !! first complex RF argument; off the negative-real cut
        complex(dp), intent(in) :: y !! second complex RF argument; off the negative-real cut
        complex(dp), intent(in) :: z !! third complex RF argument; off the negative-real cut
        complex(dp) :: value

        integer :: zero_count
        real(dp) :: scale_value

        if (complex_has_nan(x) .or. complex_has_nan(y) .or. complex_has_nan(z)) then
            value = complex_nan()
            return
        end if
        if (on_negative_real_cut(x) .or. on_negative_real_cut(y) .or. &
                on_negative_real_cut(z)) then
            value = complex_nan()
            return
        end if
        if (aimag(x) == 0.0_dp .and. aimag(y) == 0.0_dp .and. aimag(z) == 0.0_dp) then
            value = cmplx(elliprf_real(real(x, dp), real(y, dp), real(z, dp)), &
                0.0_dp, kind=dp)
            return
        end if
        if (complex_has_ambiguous_infinity(x) .or. &
                complex_has_ambiguous_infinity(y) .or. &
                complex_has_ambiguous_infinity(z)) then
            value = complex_nan()
            return
        end if
        if (complex_has_infinite(x) .or. complex_has_infinite(y) .or. &
                complex_has_infinite(z)) then
            value = cmplx(0.0_dp, 0.0_dp, kind=dp)
            return
        end if

        zero_count = merge(1, 0, complex_is_zero(x)) + &
            merge(1, 0, complex_is_zero(y)) + merge(1, 0, complex_is_zero(z))
        if (zero_count >= 2) then
            value = complex_positive_infinity()
            return
        end if

        scale_value = max(abs(x), abs(y), abs(z))
        value = carlson_rf_core_complex(x / scale_value, y / scale_value, &
            z / scale_value) / sqrt(scale_value)
    end function elliprf_complex

    pure elemental function elliprd_complex(x, y, z) result(value)
        complex(dp), intent(in) :: x !! first complex RD argument; off the negative-real cut
        complex(dp), intent(in) :: y !! second complex RD argument; off the negative-real cut
        complex(dp), intent(in) :: z !! nonzero complex distinguished RD argument
        complex(dp) :: value

        real(dp) :: scale_value

        if (complex_has_nan(x) .or. complex_has_nan(y) .or. complex_has_nan(z)) then
            value = complex_nan()
            return
        end if
        if (on_negative_real_cut(x) .or. on_negative_real_cut(y) .or. &
                on_negative_real_cut(z)) then
            value = complex_nan()
            return
        end if
        if (aimag(x) == 0.0_dp .and. aimag(y) == 0.0_dp .and. aimag(z) == 0.0_dp) then
            value = cmplx(elliprd_real(real(x, dp), real(y, dp), real(z, dp)), &
                0.0_dp, kind=dp)
            return
        end if
        if (complex_has_ambiguous_infinity(x) .or. &
                complex_has_ambiguous_infinity(y) .or. &
                complex_has_ambiguous_infinity(z)) then
            value = complex_nan()
            return
        end if
        if (complex_has_infinite(x) .or. complex_has_infinite(y) .or. &
                complex_has_infinite(z)) then
            value = cmplx(0.0_dp, 0.0_dp, kind=dp)
            return
        end if
        if (complex_is_zero(z) .or. &
                (complex_is_zero(x) .and. complex_is_zero(y))) then
            value = complex_positive_infinity()
            return
        end if

        scale_value = max(abs(x), abs(y), abs(z))
        value = carlson_rd_core_complex(x / scale_value, y / scale_value, &
            z / scale_value) / sqrt(scale_value)
        value = value / scale_value
    end function elliprd_complex

    pure elemental function elliprg_complex(x, y, z) result(value)
        complex(dp), intent(in) :: x !! first complex RG argument; off the negative-real cut
        complex(dp), intent(in) :: y !! second complex RG argument; off the negative-real cut
        complex(dp), intent(in) :: z !! third complex RG argument; off the negative-real cut
        complex(dp) :: value

        complex(dp) :: a
        complex(dp) :: b
        complex(dp) :: c
        complex(dp) :: rd_value
        complex(dp) :: rf_value
        complex(dp) :: term_root
        real(dp) :: scale_value

        if (complex_has_nan(x) .or. complex_has_nan(y) .or. complex_has_nan(z)) then
            value = complex_nan()
            return
        end if
        if (on_negative_real_cut(x) .or. on_negative_real_cut(y) .or. &
                on_negative_real_cut(z)) then
            value = complex_nan()
            return
        end if
        if (aimag(x) == 0.0_dp .and. aimag(y) == 0.0_dp .and. aimag(z) == 0.0_dp) then
            value = cmplx(elliprg_real(real(x, dp), real(y, dp), real(z, dp)), &
                0.0_dp, kind=dp)
            return
        end if
        if (complex_has_ambiguous_infinity(x) .or. &
                complex_has_ambiguous_infinity(y) .or. &
                complex_has_ambiguous_infinity(z)) then
            value = complex_nan()
            return
        end if
        if (complex_has_infinite(x) .or. complex_has_infinite(y) .or. &
                complex_has_infinite(z)) then
            value = complex_positive_infinity()
            return
        end if
        if (complex_is_zero(x) .and. complex_is_zero(y) .and. complex_is_zero(z)) then
            value = cmplx(0.0_dp, 0.0_dp, kind=dp)
            return
        end if

        scale_value = max(abs(x), abs(y), abs(z))
        call choose_rg_pivot_complex(x / scale_value, y / scale_value, &
            z / scale_value, a, b, c)
        if (complex_is_zero(a) .and. complex_is_zero(b)) then
            value = 0.5_dp * sqrt(c * scale_value)
            return
        end if

        rf_value = carlson_rf_core_complex(a, b, c)
        rd_value = carlson_rd_core_complex(a, b, c)
        if (complex_is_zero(a) .or. complex_is_zero(b)) then
            term_root = cmplx(0.0_dp, 0.0_dp, kind=dp)
        else
            ! Keep the principal roots as separate factors.  Collapsing them to
            ! sqrt(a*b/c) can select the wrong branch for complex arguments.
            term_root = sqrt(a) * sqrt(b) / sqrt(c)
        end if
        value = 0.5_dp * sqrt(scale_value) * &
            (c * rf_value - ((a - c) * (b - c) / 3.0_dp) * rd_value + term_root)
    end function elliprg_complex

    pure elemental function elliprj_complex(x, y, z, p) result(value)
        complex(dp), intent(in) :: x !! first complex RJ symmetric argument
        complex(dp), intent(in) :: y !! second complex RJ symmetric argument
        complex(dp), intent(in) :: z !! third complex RJ symmetric argument
        complex(dp), intent(in) :: p !! nonzero complex RJ pole argument
        complex(dp) :: value

        complex(dp) :: core
        real(dp) :: scale_value

        if (complex_has_nan(x) .or. complex_has_nan(y) .or. &
                complex_has_nan(z) .or. complex_has_nan(p)) then
            value = complex_nan()
            return
        end if
        if (on_negative_real_cut(x) .or. on_negative_real_cut(y) .or. &
                on_negative_real_cut(z) .or. complex_is_zero(p)) then
            value = complex_nan()
            return
        end if
        if (aimag(x) == 0.0_dp .and. aimag(y) == 0.0_dp .and. &
                aimag(z) == 0.0_dp .and. aimag(p) == 0.0_dp) then
            value = cmplx(elliprj_real(real(x, dp), real(y, dp), real(z, dp), &
                real(p, dp)), 0.0_dp, kind=dp)
            return
        end if
        if (complex_has_ambiguous_infinity(x) .or. &
                complex_has_ambiguous_infinity(y) .or. &
                complex_has_ambiguous_infinity(z)) then
            value = complex_nan()
            return
        end if
        if (complex_has_infinite(x) .or. complex_has_infinite(y) .or. &
                complex_has_infinite(z) .or. complex_has_infinite(p)) then
            value = cmplx(0.0_dp, 0.0_dp, kind=dp)
            return
        end if
        if ((complex_is_zero(x) .and. complex_is_zero(y)) .or. &
                (complex_is_zero(x) .and. complex_is_zero(z)) .or. &
                (complex_is_zero(y) .and. complex_is_zero(z))) then
            value = complex_positive_infinity()
            return
        end if
        if (on_negative_real_cut(p) .or. .not. rj_sufficient_domain(x, y, z, p)) then
            value = complex_nan()
            return
        end if

        scale_value = max(abs(x), abs(y), abs(z), abs(p))
        core = carlson_rj_integral_complex(x / scale_value, y / scale_value, &
            z / scale_value, p / scale_value)
        value = core / sqrt(scale_value)
        value = value / scale_value
    end function elliprj_complex

    pure function elliprc_positive(x, y) result(value)
        real(dp), intent(in) :: x !! finite nonnegative first RC argument
        real(dp), intent(in) :: y !! finite positive second RC argument
        real(dp) :: value

        real(dp) :: core
        real(dp) :: scale_value

        if (x == y) then
            value = 1.0_dp / sqrt(y)
            return
        end if
        scale_value = max(x, y)
        if (y / scale_value == 0.0_dp) then
            value = (log(sqrt(x) + sqrt(x - y)) - 0.5_dp * log(y)) / sqrt(x - y)
            return
        end if
        core = carlson_rc_core(x / scale_value, y / scale_value)
        value = core / sqrt(scale_value)
    end function elliprc_positive


    pure function elliprj_positive(x, y, z, p) result(value)
        real(dp), intent(in) :: x !! finite nonnegative first RJ argument
        real(dp), intent(in) :: y !! finite nonnegative second RJ argument
        real(dp), intent(in) :: z !! finite nonnegative third RJ argument
        real(dp), intent(in) :: p !! finite positive fourth RJ argument
        real(dp) :: value

        real(dp) :: core
        real(dp) :: scaled_once
        real(dp) :: scale_value

        scale_value = max(x, y, z, p)
        core = carlson_rj_core(x / scale_value, y / scale_value, &
            z / scale_value, p / scale_value)
        scaled_once = core / sqrt(scale_value)
        value = divide_by_positive_scale(scaled_once, scale_value, x)
    end function elliprj_positive

    pure function elliprj_principal_value(x, y, z, p) result(value)
        real(dp), intent(in) :: x !! smallest nonnegative symmetric RJ argument
        real(dp), intent(in) :: y !! middle nonnegative symmetric RJ argument
        real(dp), intent(in) :: z !! largest nonnegative symmetric RJ argument
        real(dp), intent(in) :: p !! negative RJ pole argument
        real(dp) :: value

        real(dp) :: correction
        real(dp) :: p_abs
        real(dp) :: p_minus_y
        real(dp) :: p_positive
        real(dp) :: rc_arg1
        real(dp) :: rc_arg2
        real(dp) :: root_factor

        p_abs = -p
        p_minus_y = (z - y) * ((y - x) / (y + p_abs))
        p_positive = y + p_minus_y
        if (p_positive <= 0.0_dp) then
            value = quiet_nan(x)
            return
        end if

        rc_arg1 = x * z + p_positive * p_abs
        rc_arg2 = p_positive * p_abs
        if (rc_arg1 == 0.0_dp .or. rc_arg2 == 0.0_dp .or. x * y * z == 0.0_dp) then
            root_factor = 0.0_dp
        else
            root_factor = sqrt((x * y * z) / rc_arg1)
        end if

        correction = p_minus_y * elliprj_positive(x, y, z, p_positive) - &
            3.0_dp * elliprf(x, y, z) + &
            3.0_dp * root_factor * elliprc(rc_arg1, rc_arg2)
        value = correction / (y + p_abs)
    end function elliprj_principal_value

    pure function divide_by_positive_scale(numerator, scale_value, kind_arg) result(value)
        real(dp), intent(in) :: numerator !! finite signed numerator before division
        real(dp), intent(in) :: scale_value !! finite positive homogeneity scale
        real(dp), intent(in) :: kind_arg !! real-kind selector for positive infinity
        real(dp) :: value

        if (scale_value < 1.0_dp) then
            if (numerator > huge(1.0_dp) * scale_value) then
                value = positive_infinity(kind_arg)
                return
            end if
            if (numerator < -huge(1.0_dp) * scale_value) then
                value = -positive_infinity(kind_arg)
                return
            end if
        end if
        value = numerator / scale_value
    end function divide_by_positive_scale

    pure subroutine sort_three(x, y, z, a, b, c)
        real(dp), intent(in) :: x !! first value to sort
        real(dp), intent(in) :: y !! second value to sort
        real(dp), intent(in) :: z !! third value to sort
        real(dp), intent(out) :: a !! smallest sorted value
        real(dp), intent(out) :: b !! middle sorted value
        real(dp), intent(out) :: c !! largest sorted value

        real(dp) :: t1
        real(dp) :: t2
        real(dp) :: t3

        t1 = x
        t2 = y
        t3 = z
        if (t1 > t2) call swap_values(t1, t2)
        if (t2 > t3) call swap_values(t2, t3)
        if (t1 > t2) call swap_values(t1, t2)
        a = t1
        b = t2
        c = t3
    end subroutine sort_three

    pure subroutine choose_rg_pivot(x, y, z, a, b, c)
        real(dp), intent(in) :: x !! normalized first RG argument
        real(dp), intent(in) :: y !! normalized second RG argument
        real(dp), intent(in) :: z !! normalized third RG argument
        real(dp), intent(out) :: a !! first non-pivot argument
        real(dp), intent(out) :: b !! second non-pivot argument
        real(dp), intent(out) :: c !! positive pivot chosen as the largest argument

        if (z >= x .and. z >= y) then
            a = x
            b = y
            c = z
        else if (y >= x) then
            a = x
            b = z
            c = y
        else
            a = y
            b = z
            c = x
        end if
    end subroutine choose_rg_pivot

    pure subroutine swap_values(a, b)
        real(dp), intent(inout) :: a !! first scalar, replaced by second on return
        real(dp), intent(inout) :: b !! second scalar, replaced by first on return

        real(dp) :: temp

        temp = a
        a = b
        b = temp
    end subroutine swap_values

    pure function carlson_rc_core(x, y) result(value)
        real(dp), intent(in) :: x !! normalized nonnegative RC argument
        real(dp), intent(in) :: y !! normalized positive RC argument
        real(dp) :: value

        integer :: iter
        real(dp), parameter :: c1 = 1.0_dp / 7.0_dp
        real(dp), parameter :: c2 = 9.0_dp / 22.0_dp
        real(dp), parameter :: errtol = (unit_roundoff / 16.0_dp)**(1.0_dp / 6.0_dp)
        real(dp) :: lambda
        real(dp) :: mu
        real(dp) :: s
        real(dp) :: sn
        real(dp) :: xn
        real(dp) :: yn

        xn = x
        yn = y
        do iter = 1, max_iterations
            mu = (xn + 2.0_dp * yn) / 3.0_dp
            sn = (yn + mu) / mu - 2.0_dp
            if (abs(sn) < errtol) then
                s = sn * sn * (0.30_dp + sn * (c1 + sn * (0.375_dp + sn * c2)))
                value = (1.0_dp + s) / sqrt(mu)
                return
            end if
            lambda = 2.0_dp * sqrt(xn) * sqrt(yn) + yn
            xn = 0.25_dp * (xn + lambda)
            yn = 0.25_dp * (yn + lambda)
        end do
        value = quiet_nan(x)
    end function carlson_rc_core

    pure function carlson_rf_core(x, y, z) result(value)
        real(dp), intent(in) :: x !! normalized nonnegative RF argument
        real(dp), intent(in) :: y !! normalized nonnegative RF argument
        real(dp), intent(in) :: z !! normalized nonnegative RF argument
        real(dp) :: value

        integer :: iter
        real(dp), parameter :: c1 = 1.0_dp / 24.0_dp
        real(dp), parameter :: c2 = 3.0_dp / 44.0_dp
        real(dp), parameter :: c3 = 1.0_dp / 14.0_dp
        real(dp), parameter :: errtol = (4.0_dp * unit_roundoff)**(1.0_dp / 6.0_dp)
        real(dp) :: e2
        real(dp) :: e3
        real(dp) :: epsilon_value
        real(dp) :: lambda
        real(dp) :: mu
        real(dp) :: s
        real(dp) :: xn
        real(dp) :: xndev
        real(dp) :: yn
        real(dp) :: yndev
        real(dp) :: zn
        real(dp) :: zndev

        xn = x
        yn = y
        zn = z
        do iter = 1, max_iterations
            mu = (xn + yn + zn) / 3.0_dp
            xndev = (mu - xn) / mu
            yndev = (mu - yn) / mu
            zndev = (mu - zn) / mu
            epsilon_value = max(abs(xndev), abs(yndev), abs(zndev))
            if (epsilon_value < errtol) then
                e2 = xndev * yndev - zndev * zndev
                e3 = xndev * yndev * zndev
                s = 1.0_dp + (c1 * e2 - 0.10_dp - c2 * e3) * e2 + c3 * e3
                value = s / sqrt(mu)
                return
            end if
            lambda = sqrt(xn) * (sqrt(yn) + sqrt(zn)) + sqrt(yn) * sqrt(zn)
            xn = 0.25_dp * (xn + lambda)
            yn = 0.25_dp * (yn + lambda)
            zn = 0.25_dp * (zn + lambda)
        end do
        value = quiet_nan(x)
    end function carlson_rf_core

    pure function carlson_rd_core(x, y, z) result(value)
        real(dp), intent(in) :: x !! normalized nonnegative RD argument
        real(dp), intent(in) :: y !! normalized nonnegative RD argument
        real(dp), intent(in) :: z !! normalized positive RD argument
        real(dp) :: value

        integer :: iter
        real(dp), parameter :: c1 = 3.0_dp / 14.0_dp
        real(dp), parameter :: c2 = 1.0_dp / 6.0_dp
        real(dp), parameter :: c3 = 9.0_dp / 22.0_dp
        real(dp), parameter :: c4 = 3.0_dp / 26.0_dp
        real(dp), parameter :: errtol = (unit_roundoff / 3.0_dp)**(1.0_dp / 6.0_dp)
        real(dp) :: ea
        real(dp) :: eb
        real(dp) :: ec
        real(dp) :: ed
        real(dp) :: ef
        real(dp) :: epsilon_value
        real(dp) :: lambda
        real(dp) :: mu
        real(dp) :: power4
        real(dp) :: s1
        real(dp) :: s2
        real(dp) :: sigma
        real(dp) :: xn
        real(dp) :: xndev
        real(dp) :: yn
        real(dp) :: yndev
        real(dp) :: zn
        real(dp) :: zndev

        xn = x
        yn = y
        zn = z
        sigma = 0.0_dp
        power4 = 1.0_dp
        do iter = 1, max_iterations
            mu = 0.20_dp * (xn + yn + 3.0_dp * zn)
            xndev = (mu - xn) / mu
            yndev = (mu - yn) / mu
            zndev = (mu - zn) / mu
            epsilon_value = max(abs(xndev), abs(yndev), abs(zndev))
            if (epsilon_value < errtol) then
                ea = xndev * yndev
                eb = zndev * zndev
                ec = ea - eb
                ed = ea - 6.0_dp * eb
                ef = ed + 2.0_dp * ec
                s1 = ed * (-c1 + 0.25_dp * c3 * ed - 1.50_dp * c4 * zndev * ef)
                s2 = zndev * (c2 * ef + zndev * (-c3 * ec + zndev * c4 * ea))
                value = 3.0_dp * sigma + power4 * (1.0_dp + s1 + s2) / (mu * sqrt(mu))
                return
            end if
            lambda = sqrt(xn) * (sqrt(yn) + sqrt(zn)) + sqrt(yn) * sqrt(zn)
            sigma = sigma + power4 / (sqrt(zn) * (zn + lambda))
            power4 = 0.25_dp * power4
            xn = 0.25_dp * (xn + lambda)
            yn = 0.25_dp * (yn + lambda)
            zn = 0.25_dp * (zn + lambda)
        end do
        value = quiet_nan(x)
    end function carlson_rd_core

    pure function carlson_rj_core(x, y, z, p) result(value)
        real(dp), intent(in) :: x !! normalized nonnegative RJ argument
        real(dp), intent(in) :: y !! normalized nonnegative RJ argument
        real(dp), intent(in) :: z !! normalized nonnegative RJ argument
        real(dp), intent(in) :: p !! normalized positive RJ argument
        real(dp) :: value

        integer :: iter
        real(dp), parameter :: c1 = 3.0_dp / 14.0_dp
        real(dp), parameter :: c2 = 1.0_dp / 3.0_dp
        real(dp), parameter :: c3 = 3.0_dp / 22.0_dp
        real(dp), parameter :: c4 = 3.0_dp / 26.0_dp
        real(dp), parameter :: errtol = (unit_roundoff / 3.0_dp)**(1.0_dp / 6.0_dp)
        real(dp) :: alpha
        real(dp) :: beta
        real(dp) :: ea
        real(dp) :: eb
        real(dp) :: ec
        real(dp) :: e2
        real(dp) :: e3
        real(dp) :: epsilon_value
        real(dp) :: lambda
        real(dp) :: mu
        real(dp) :: pn
        real(dp) :: pndev
        real(dp) :: power4
        real(dp) :: s1
        real(dp) :: s2
        real(dp) :: s3
        real(dp) :: sigma
        real(dp) :: xn
        real(dp) :: xndev
        real(dp) :: yn
        real(dp) :: yndev
        real(dp) :: zn
        real(dp) :: zndev

        xn = x
        yn = y
        zn = z
        pn = p
        sigma = 0.0_dp
        power4 = 1.0_dp
        do iter = 1, max_iterations
            mu = 0.20_dp * (xn + yn + zn + 2.0_dp * pn)
            xndev = (mu - xn) / mu
            yndev = (mu - yn) / mu
            zndev = (mu - zn) / mu
            pndev = (mu - pn) / mu
            epsilon_value = max(abs(xndev), abs(yndev), abs(zndev), abs(pndev))
            if (epsilon_value < errtol) then
                ea = xndev * (yndev + zndev) + yndev * zndev
                eb = xndev * yndev * zndev
                ec = pndev * pndev
                e2 = ea - 3.0_dp * ec
                e3 = eb + 2.0_dp * pndev * (ea - ec)
                s1 = 1.0_dp + e2 * (-c1 + 0.75_dp * c3 * e2 - 1.50_dp * c4 * e3)
                s2 = eb * (0.50_dp * c2 + pndev * (-2.0_dp * c3 + pndev * c4))
                s3 = pndev * ea * (c2 - pndev * c3) - c2 * pndev * ec
                value = 3.0_dp * sigma + power4 * (s1 + s2 + s3) / (mu * sqrt(mu))
                return
            end if
            lambda = sqrt(xn) * (sqrt(yn) + sqrt(zn)) + sqrt(yn) * sqrt(zn)
            alpha = pn * (sqrt(xn) + sqrt(yn) + sqrt(zn)) + &
                sqrt(xn) * sqrt(yn) * sqrt(zn)
            alpha = alpha * alpha
            beta = pn * (pn + lambda) * (pn + lambda)
            sigma = sigma + power4 * carlson_rc_core(alpha, beta)
            power4 = 0.25_dp * power4
            xn = 0.25_dp * (xn + lambda)
            yn = 0.25_dp * (yn + lambda)
            zn = 0.25_dp * (zn + lambda)
            pn = 0.25_dp * (pn + lambda)
        end do
        value = quiet_nan(x)
    end function carlson_rj_core


    pure function carlson_rc_core_complex(x, y) result(value)
        complex(dp), intent(in) :: x !! normalized complex first RC argument
        complex(dp), intent(in) :: y !! normalized nonzero complex second RC argument
        complex(dp) :: value

        integer :: iter
        real(dp), parameter :: c1 = 1.0_dp / 7.0_dp
        real(dp), parameter :: c2 = 9.0_dp / 22.0_dp
        real(dp), parameter :: errtol = (unit_roundoff / 16.0_dp)**(1.0_dp / 6.0_dp)
        complex(dp) :: lambda
        complex(dp) :: mu
        complex(dp) :: s
        complex(dp) :: sn
        complex(dp) :: xn
        complex(dp) :: yn

        xn = x
        yn = y
        do iter = 1, max_iterations
            mu = (xn + 2.0_dp * yn) / 3.0_dp
            if (complex_is_zero(mu)) then
                value = complex_nan()
                return
            end if
            sn = (yn + mu) / mu - 2.0_dp
            if (abs(sn) < errtol) then
                s = sn * sn * (0.30_dp + sn * (c1 + sn * (0.375_dp + sn * c2)))
                value = (1.0_dp + s) / sqrt(mu)
                return
            end if
            lambda = 2.0_dp * sqrt(xn) * sqrt(yn) + yn
            xn = 0.25_dp * (xn + lambda)
            yn = 0.25_dp * (yn + lambda)
        end do
        value = complex_nan()
    end function carlson_rc_core_complex

    pure function carlson_rf_core_complex(x, y, z) result(value)
        complex(dp), intent(in) :: x !! normalized complex first RF argument
        complex(dp), intent(in) :: y !! normalized complex second RF argument
        complex(dp), intent(in) :: z !! normalized complex third RF argument
        complex(dp) :: value

        complex(dp) :: e2
        complex(dp) :: e3
        complex(dp) :: lambda
        complex(dp) :: mu
        complex(dp) :: s
        complex(dp) :: xn
        complex(dp) :: xndev
        complex(dp) :: yn
        complex(dp) :: yndev
        complex(dp) :: zn
        complex(dp) :: zndev
        integer :: iter
        real(dp), parameter :: c1 = 1.0_dp / 24.0_dp
        real(dp), parameter :: c2 = 3.0_dp / 44.0_dp
        real(dp), parameter :: c3 = 1.0_dp / 14.0_dp
        real(dp), parameter :: errtol = (4.0_dp * unit_roundoff)**(1.0_dp / 6.0_dp)
        real(dp) :: epsilon_value

        xn = x
        yn = y
        zn = z
        do iter = 1, max_iterations
            mu = (xn + yn + zn) / 3.0_dp
            if (complex_is_zero(mu)) then
                value = complex_nan()
                return
            end if
            xndev = (mu - xn) / mu
            yndev = (mu - yn) / mu
            zndev = (mu - zn) / mu
            epsilon_value = max(abs(xndev), abs(yndev), abs(zndev))
            if (epsilon_value < errtol) then
                e2 = xndev * yndev - zndev * zndev
                e3 = xndev * yndev * zndev
                s = 1.0_dp + (c1 * e2 - 0.10_dp - c2 * e3) * e2 + c3 * e3
                value = s / sqrt(mu)
                return
            end if
            lambda = sqrt(xn) * (sqrt(yn) + sqrt(zn)) + sqrt(yn) * sqrt(zn)
            xn = 0.25_dp * (xn + lambda)
            yn = 0.25_dp * (yn + lambda)
            zn = 0.25_dp * (zn + lambda)
        end do
        value = complex_nan()
    end function carlson_rf_core_complex

    pure function carlson_rd_core_complex(x, y, z) result(value)
        complex(dp), intent(in) :: x !! normalized complex first RD argument
        complex(dp), intent(in) :: y !! normalized complex second RD argument
        complex(dp), intent(in) :: z !! normalized nonzero complex distinguished RD argument
        complex(dp) :: value

        complex(dp) :: ea
        complex(dp) :: eb
        complex(dp) :: ec
        complex(dp) :: ed
        complex(dp) :: ef
        complex(dp) :: lambda
        complex(dp) :: mu
        complex(dp) :: s1
        complex(dp) :: s2
        complex(dp) :: sigma
        complex(dp) :: xn
        complex(dp) :: xndev
        complex(dp) :: yn
        complex(dp) :: yndev
        complex(dp) :: zn
        complex(dp) :: zndev
        integer :: iter
        real(dp), parameter :: c1 = 3.0_dp / 14.0_dp
        real(dp), parameter :: c2 = 1.0_dp / 6.0_dp
        real(dp), parameter :: c3 = 9.0_dp / 22.0_dp
        real(dp), parameter :: c4 = 3.0_dp / 26.0_dp
        real(dp), parameter :: errtol = (unit_roundoff / 3.0_dp)**(1.0_dp / 6.0_dp)
        real(dp) :: epsilon_value
        real(dp) :: power4

        xn = x
        yn = y
        zn = z
        sigma = cmplx(0.0_dp, 0.0_dp, kind=dp)
        power4 = 1.0_dp
        do iter = 1, max_iterations
            mu = 0.20_dp * (xn + yn + 3.0_dp * zn)
            if (complex_is_zero(mu)) then
                value = complex_nan()
                return
            end if
            xndev = (mu - xn) / mu
            yndev = (mu - yn) / mu
            zndev = (mu - zn) / mu
            epsilon_value = max(abs(xndev), abs(yndev), abs(zndev))
            if (epsilon_value < errtol) then
                ea = xndev * yndev
                eb = zndev * zndev
                ec = ea - eb
                ed = ea - 6.0_dp * eb
                ef = ed + 2.0_dp * ec
                s1 = ed * (-c1 + 0.25_dp * c3 * ed - 1.50_dp * c4 * zndev * ef)
                s2 = zndev * (c2 * ef + zndev * (-c3 * ec + zndev * c4 * ea))
                value = 3.0_dp * sigma + power4 * (1.0_dp + s1 + s2) / &
                    (mu * sqrt(mu))
                return
            end if
            lambda = sqrt(xn) * (sqrt(yn) + sqrt(zn)) + sqrt(yn) * sqrt(zn)
            sigma = sigma + power4 / (sqrt(zn) * (zn + lambda))
            power4 = 0.25_dp * power4
            xn = 0.25_dp * (xn + lambda)
            yn = 0.25_dp * (yn + lambda)
            zn = 0.25_dp * (zn + lambda)
        end do
        value = complex_nan()
    end function carlson_rd_core_complex

    pure function carlson_rj_integral_complex(x, y, z, p) result(value)
        complex(dp), intent(in) :: x !! normalized complex first RJ symmetric argument
        complex(dp), intent(in) :: y !! normalized complex second RJ symmetric argument
        complex(dp), intent(in) :: z !! normalized complex third RJ symmetric argument
        complex(dp), intent(in) :: p !! normalized nonzero complex pole argument
        complex(dp) :: value

        integer, parameter :: max_depth = 32
        integer, parameter :: max_subintervals = 2048
        real(dp), parameter :: absolute_tolerance = 5.0e-15_dp
        real(dp), parameter :: relative_tolerance = 5.0e-14_dp
        real(dp) :: a
        real(dp) :: b
        complex(dp) :: fine
        complex(dp) :: left_value
        complex(dp) :: right_value
        complex(dp) :: whole
        integer :: depth
        integer :: depth_stack(max_subintervals)
        real(dp) :: left_stack(max_subintervals)
        real(dp) :: midpoint
        real(dp) :: right_stack(max_subintervals)
        integer :: top

        value = cmplx(0.0_dp, 0.0_dp, kind=dp)
        top = 1
        left_stack(1) = 0.0_dp
        right_stack(1) = 1.0_dp
        depth_stack(1) = 0

        do while (top > 0)
            a = left_stack(top)
            b = right_stack(top)
            depth = depth_stack(top)
            top = top - 1
            midpoint = 0.5_dp * (a + b)

            whole = carlson_rj_gl16_complex(a, b, x, y, z, p)
            left_value = carlson_rj_gl16_complex(a, midpoint, x, y, z, p)
            right_value = carlson_rj_gl16_complex(midpoint, b, x, y, z, p)
            fine = left_value + right_value

            if (abs(fine - whole) <= absolute_tolerance * (b - a) + &
                    relative_tolerance * abs(fine)) then
                value = value + fine
            else
                if (depth >= max_depth .or. top + 2 > max_subintervals) then
                    value = complex_nan()
                    return
                end if
                top = top + 1
                left_stack(top) = midpoint
                right_stack(top) = b
                depth_stack(top) = depth + 1
                top = top + 1
                left_stack(top) = a
                right_stack(top) = midpoint
                depth_stack(top) = depth + 1
            end if
        end do
    end function carlson_rj_integral_complex

    pure function carlson_rj_gl16_complex(a, b, x, y, z, p) result(value)
        real(dp), intent(in) :: a !! lower transformed quadrature endpoint in [0,1]
        real(dp), intent(in) :: b !! upper transformed quadrature endpoint in [0,1]
        complex(dp), intent(in) :: x !! normalized complex first RJ argument
        complex(dp), intent(in) :: y !! normalized complex second RJ argument
        complex(dp), intent(in) :: z !! normalized complex third RJ argument
        complex(dp), intent(in) :: p !! normalized nonzero complex pole argument
        complex(dp) :: value

        real(dp), parameter :: nodes(8) = [ &
            0.0950125098376374401853193354250_dp, &
            0.281603550779258913230460501460_dp, &
            0.458016777657227386342419442984_dp, &
            0.617876244402643748446671764049_dp, &
            0.755404408355003033895101194847_dp, &
            0.865631202387831743880467897712_dp, &
            0.944575023073232576077988415535_dp, &
            0.989400934991649932596154173450_dp]
        real(dp), parameter :: weights(8) = [ &
            0.189450610455068496285396723208_dp, &
            0.182603415044923588866763667969_dp, &
            0.169156519395002538189312079030_dp, &
            0.149595988816576732081501730547_dp, &
            0.124628971255533872052476282192_dp, &
            0.0951585116824927848099251076022_dp, &
            0.0622535239386478928628438369944_dp, &
            0.0271524594117540948517805724560_dp]
        real(dp) :: half_width
        integer :: i
        real(dp) :: midpoint
        real(dp) :: u_left
        real(dp) :: u_right

        midpoint = 0.5_dp * (a + b)
        half_width = 0.5_dp * (b - a)
        value = cmplx(0.0_dp, 0.0_dp, kind=dp)
        do i = 1, size(nodes)
            u_left = midpoint - half_width * nodes(i)
            u_right = midpoint + half_width * nodes(i)
            value = value + weights(i) * &
                (carlson_rj_transformed_integrand(u_left, x, y, z, p) + &
                carlson_rj_transformed_integrand(u_right, x, y, z, p))
        end do
        value = half_width * value
    end function carlson_rj_gl16_complex

    pure function carlson_rj_transformed_integrand(u, x, y, z, p) result(value)
        real(dp), intent(in) :: u !! transformed integration coordinate strictly between zero and one
        complex(dp), intent(in) :: x !! normalized complex first RJ argument
        complex(dp), intent(in) :: y !! normalized complex second RJ argument
        complex(dp), intent(in) :: z !! normalized complex third RJ argument
        complex(dp), intent(in) :: p !! normalized nonzero complex pole argument
        complex(dp) :: value

        real(dp) :: one_minus_u
        real(dp) :: q
        real(dp) :: t

        one_minus_u = 1.0_dp - u
        q = u / one_minus_u
        t = q * q
        value = 3.0_dp * q / (one_minus_u * one_minus_u) / &
            (sqrt(cmplx(t, 0.0_dp, kind=dp) + x) * &
            sqrt(cmplx(t, 0.0_dp, kind=dp) + y) * &
            sqrt(cmplx(t, 0.0_dp, kind=dp) + z) * &
            (cmplx(t, 0.0_dp, kind=dp) + p))
    end function carlson_rj_transformed_integrand

    pure subroutine choose_rg_pivot_complex(x, y, z, a, b, c)
        complex(dp), intent(in) :: x !! normalized first complex RG argument
        complex(dp), intent(in) :: y !! normalized second complex RG argument
        complex(dp), intent(in) :: z !! normalized third complex RG argument
        complex(dp), intent(out) :: a !! first non-pivot argument
        complex(dp), intent(out) :: b !! second non-pivot argument
        complex(dp), intent(out) :: c !! nonzero pivot with largest magnitude

        if (abs(z) >= abs(x) .and. abs(z) >= abs(y)) then
            a = x
            b = y
            c = z
        else if (abs(y) >= abs(x)) then
            a = x
            b = z
            c = y
        else
            a = y
            b = z
            c = x
        end if
    end subroutine choose_rg_pivot_complex

    pure elemental logical function rj_sufficient_domain(x, y, z, p) result(is_valid)
        complex(dp), intent(in) :: x !! finite first RJ symmetric argument
        complex(dp), intent(in) :: y !! finite second RJ symmetric argument
        complex(dp), intent(in) :: z !! finite third RJ symmetric argument
        complex(dp), intent(in) :: p !! finite nonzero RJ pole argument

        ! Carlson's sufficient complex-domain conditions as exposed by SciPy:
        ! real nonnegative x/y/z and the conjugate-pair case permit a general
        ! off-pole p, while a genuinely complex generic triple requires Re(p)>0.
        if (aimag(x) == 0.0_dp .and. aimag(y) == 0.0_dp .and. &
                aimag(z) == 0.0_dp .and. real(x, dp) >= 0.0_dp .and. &
                real(y, dp) >= 0.0_dp .and. real(z, dp) >= 0.0_dp) then
            is_valid = .true.
            return
        end if

        is_valid = (y == conjg(x) .and. aimag(z) == 0.0_dp .and. real(z, dp) >= 0.0_dp) .or. &
            (z == conjg(x) .and. aimag(y) == 0.0_dp .and. real(y, dp) >= 0.0_dp) .or. &
            (z == conjg(y) .and. aimag(x) == 0.0_dp .and. real(x, dp) >= 0.0_dp)
        if (is_valid) return

        is_valid = real(x, dp) >= 0.0_dp .and. real(y, dp) >= 0.0_dp .and. &
            real(z, dp) >= 0.0_dp .and. real(p, dp) > 0.0_dp
    end function rj_sufficient_domain

    pure elemental logical function on_negative_real_cut(z) result(on_cut)
        complex(dp), intent(in) :: z !! finite complex value to test against the principal cut

        on_cut = aimag(z) == 0.0_dp .and. real(z, dp) < 0.0_dp
    end function on_negative_real_cut

    pure elemental logical function complex_is_zero(z) result(is_zero)
        complex(dp), intent(in) :: z !! finite complex value to test for exact zero

        is_zero = real(z, dp) == 0.0_dp .and. aimag(z) == 0.0_dp
    end function complex_is_zero

    pure elemental logical function complex_has_nan(z) result(has_nan)
        complex(dp), intent(in) :: z !! complex value whose components are checked for NaN

        has_nan = ieee_is_nan(real(z, dp)) .or. ieee_is_nan(aimag(z))
    end function complex_has_nan

    pure elemental logical function complex_has_infinite(z) result(has_infinite)
        complex(dp), intent(in) :: z !! complex value whose components are checked for infinity

        has_infinite = (.not. ieee_is_finite(real(z, dp)) .and. &
            .not. ieee_is_nan(real(z, dp))) .or. &
            (.not. ieee_is_finite(aimag(z)) .and. .not. ieee_is_nan(aimag(z)))
    end function complex_has_infinite

    pure elemental logical function complex_has_ambiguous_infinity(z) result(is_ambiguous)
        complex(dp), intent(in) :: z !! complex value tested for two infinite components

        is_ambiguous = .not. ieee_is_finite(real(z, dp)) .and. &
            .not. ieee_is_nan(real(z, dp)) .and. .not. ieee_is_finite(aimag(z)) .and. &
            .not. ieee_is_nan(aimag(z))
    end function complex_has_ambiguous_infinity

    pure function complex_nan() result(value)
        complex(dp) :: value

        value = cmplx(quiet_nan(0.0_dp), quiet_nan(0.0_dp), kind=dp)
    end function complex_nan

    pure function complex_positive_infinity() result(value)
        complex(dp) :: value

        value = cmplx(positive_infinity(0.0_dp), 0.0_dp, kind=dp)
    end function complex_positive_infinity

end module scifort_carlson_elliptic
