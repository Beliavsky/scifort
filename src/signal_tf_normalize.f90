! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Transfer-function coefficient normalization compatible with signal.normalize.
module scifort_signal_tf_normalize
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    implicit none
    private
    integer, parameter, public :: normalize_status_success = 0
    integer, parameter, public :: normalize_status_invalid_input = 1

    type, public :: normalized_tf_result
        complex(dp), allocatable :: b(:,:) !! numerator(s), one polynomial per row
        complex(dp), allocatable :: a(:) !! monic denominator, descending powers
        integer :: trimmed = 0 !! count of discarded leading numerator columns
        integer :: status = normalize_status_invalid_input !! zero on success
    end type normalized_tf_result

    interface normalize
        module procedure normalize_rr, normalize_cc, normalize_rc, normalize_cr
        module procedure normalize_rr2, normalize_cc2, normalize_rc2, normalize_cr2
    end interface normalize
    public :: normalize
contains

    pure logical function is_finite(z) result(ok)
        complex(dp), intent(in) :: z(:) !! complex vector
        ok = all(ieee_is_finite(real(z,dp))) .and. all(ieee_is_finite(aimag(z)))
    end function is_finite

    function normalize_rr(b,a) result(out)
        real(dp), intent(in) :: b(:) !! real numerator coefficients
        real(dp), intent(in) :: a(:) !! real denominator coefficients
        type(normalized_tf_result) :: out
        complex(dp) :: rows(1,size(b))
        rows(1,:) = cmplx(b,0.0_dp,dp)
        out = normalize_internal(rows,cmplx(a,0.0_dp,dp))
    end function normalize_rr

    function normalize_cc(b,a) result(out)
        complex(dp), intent(in) :: b(:) !! numerator coefficients
        complex(dp), intent(in) :: a(:) !! denominator coefficients
        type(normalized_tf_result) :: out
        complex(dp) :: rows(1,size(b))
        rows(1,:) = b
        out = normalize_internal(rows,a)
    end function normalize_cc

    function normalize_rr2(b,a) result(out)
        real(dp), intent(in) :: b(:,:) !! numerator coefficients, per output row
        real(dp), intent(in) :: a(:) !! denominator coefficients
        type(normalized_tf_result) :: out
        out = normalize_internal(cmplx(b,0.0_dp,dp),cmplx(a,0.0_dp,dp))
    end function normalize_rr2

    function normalize_cc2(b,a) result(out)
        complex(dp), intent(in) :: b(:,:) !! numerator coefficients, per output row
        complex(dp), intent(in) :: a(:) !! denominator coefficients
        type(normalized_tf_result) :: out
        out = normalize_internal(b,a)
    end function normalize_cc2


    function normalize_rc(b,a) result(out)
        real(dp), intent(in) :: b(:) !! real numerator coefficients
        complex(dp), intent(in) :: a(:) !! complex denominator coefficients
        type(normalized_tf_result) :: out
        complex(dp) :: rows(1,size(b))
        rows(1,:) = cmplx(b,0.0_dp,dp)
        out = normalize_internal(rows,a)
    end function normalize_rc

    function normalize_cr(b,a) result(out)
        complex(dp), intent(in) :: b(:) !! complex numerator coefficients
        real(dp), intent(in) :: a(:) !! real denominator coefficients
        type(normalized_tf_result) :: out
        complex(dp) :: rows(1,size(b))
        rows(1,:) = b
        out = normalize_internal(rows,cmplx(a,0.0_dp,dp))
    end function normalize_cr

    function normalize_rc2(b,a) result(out)
        real(dp), intent(in) :: b(:,:) !! real numerator polynomials per row
        complex(dp), intent(in) :: a(:) !! complex denominator
        type(normalized_tf_result) :: out
        out = normalize_internal(cmplx(b,0.0_dp,dp),a)
    end function normalize_rc2

    function normalize_cr2(b,a) result(out)
        complex(dp), intent(in) :: b(:,:) !! complex numerator polynomials per row
        real(dp), intent(in) :: a(:) !! real denominator
        type(normalized_tf_result) :: out
        out = normalize_internal(b,cmplx(a,0.0_dp,dp))
    end function normalize_cr2

    function normalize_internal(b,a) result(out)
        complex(dp), intent(in) :: b(:,:) !! numerator rows
        complex(dp), intent(in) :: a(:) !! denominator
        type(normalized_tf_result) :: out
        integer :: i, start, trim_count
        complex(dp) :: lead

        allocate(out%b(0,0),out%a(0))
        if (size(a)<1 .or. size(b,1)<1 .or. size(b,2)<1) return
        if (.not. is_finite(a) .or. .not. is_finite(reshape(b,[size(b)]))) return
        start = 1
        do while (start < size(a))
            if (abs(a(start)) > 0.0_dp) exit
            start = start+1
        end do
        if (abs(a(start)) <= 0.0_dp) return
        lead = a(start)
        out%a = a(start:)/lead
        trim_count = 0
        do i = 1,size(b,2)-1
            if (any(abs(b(:,i)/lead)>1.0e-14_dp)) exit
            trim_count = trim_count+1
        end do
        out%b = b(:,trim_count+1:)/lead
        out%trimmed = trim_count
        out%status = normalize_status_success
    end function normalize_internal
end module scifort_signal_tf_normalize
