! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! State-space matrix shape normalization (scipy.signal.abcd_normalize).
module scifort_signal_abcd_normalize
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_signal_state_space, only : state_space_result, state_space_success, &
        state_space_invalid_input
    implicit none
    private
    public :: abcd_normalize
contains

    function abcd_normalize(a,b,c,d) result(out)
        class(*), optional, intent(in) :: a(:,:) !! state matrix (n,n), optional
        class(*), optional, intent(in) :: b(:,:) !! input matrix (n,p), optional
        class(*), optional, intent(in) :: c(:,:) !! output matrix (q,n), optional
        class(*), optional, intent(in) :: d(:,:) !! feedthrough matrix (q,p), optional
        type(state_space_result) :: out
        complex(dp), allocatable :: aa(:,:), bb(:,:), cc(:,:), dd(:,:)
        logical :: ok
        integer :: n, p, q

        if (.not. present(a) .and. .not. present(b) .and. .not. present(c)) return
        if (.not. present(b) .and. .not. present(d)) return
        if (.not. present(c) .and. .not. present(d)) return

        call convert_matrix(a, aa, ok)
        if (.not. ok) return
        call convert_matrix(b, bb, ok)
        if (.not. ok) return
        call convert_matrix(c, cc, ok)
        if (.not. ok) return
        call convert_matrix(d, dd, ok)
        if (.not. ok) return

        n = size(aa,1)
        if (n == 0) n = size(bb,1)
        if (n == 0) n = size(cc,2)
        p = size(bb,2)
        if (p == 0) p = size(dd,2)
        q = size(cc,1)
        if (q == 0) q = size(dd,1)

        if (size(aa) > 0) then
            if (size(aa,1)/=n .or. size(aa,2)/=n) return
        end if
        if (size(bb) > 0) then
            if (size(bb,1)/=n .or. size(bb,2)/=p) return
        end if
        if (size(cc) > 0) then
            if (size(cc,1)/=q .or. size(cc,2)/=n) return
        end if
        if (size(dd) > 0) then
            if (size(dd,1)/=q .or. size(dd,2)/=p) return
        end if
        allocate(out%a(n,n),out%b(n,p),out%c(q,n),out%d(q,p))
        out%a = (0.0_dp,0.0_dp)
        out%b = (0.0_dp,0.0_dp)
        out%c = (0.0_dp,0.0_dp)
        out%d = (0.0_dp,0.0_dp)
        if (size(aa)>0) out%a = aa
        if (size(bb)>0) out%b = bb
        if (size(cc)>0) out%c = cc
        if (size(dd)>0) out%d = dd
        out%status = state_space_success
    end function abcd_normalize

    subroutine convert_matrix(mat, out, ok)
        class(*), optional, intent(in) :: mat(:,:) !! matrix of supported intrinsic numeric type
        complex(dp), allocatable, intent(out) :: out(:,:) !! converted values, (0,0) if absent
        logical, intent(out) :: ok !! true for a numeric finite matrix

        ok = .true.
        if (.not. present(mat)) then
            allocate(out(0,0))
            return
        end if
        select type(mat)
        type is (real(dp))
            out = cmplx(mat,0.0_dp,dp)
        type is (complex(dp))
            out = mat
        type is (real)
            out = cmplx(mat,0.0_dp,dp)
        type is (complex)
            out = cmplx(mat,kind=dp)
        type is (integer)
            out = cmplx(real(mat,dp),0.0_dp,dp)
        class default
            ok = .false.
            allocate(out(0,0))
            return
        end select
        ok = all(ieee_is_finite(real(out,dp))) .and. all(ieee_is_finite(aimag(out)))
    end subroutine convert_matrix
end module scifort_signal_abcd_normalize
