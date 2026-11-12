! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Conversion of an eigendecomposition of a real matrix from complex
! eigenvectors to a real invariant basis. The change of basis follows
! directly from the real and imaginary parts of A*v = lambda*v.
module scifort_linalg_cdf2rdf
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    implicit none
    private
    public :: cdf2rdf

    interface cdf2rdf
        module procedure cdf2rdf_complex
        module procedure cdf2rdf_complex_batch
        module procedure cdf2rdf_real
        module procedure cdf2rdf_real_batch
        module procedure cdf2rdf_real_values
        module procedure cdf2rdf_real_values_batch
    end interface cdf2rdf

contains

    pure subroutine cdf2rdf_complex(w, v, wr, vr, status)
        complex(dp), intent(in) :: w(:) !! eigenvalues of a real n-by-n matrix
        complex(dp), intent(in) :: v(:, :) !! corresponding right eigenvectors, columns (n,n)
        real(dp), intent(out) :: wr(:, :) !! real n-by-n block eigenvalue matrix
        real(dp), intent(out) :: vr(:, :) !! real n-by-n invariant-basis matrix
        integer, intent(out) :: status !! success, invalid shape/pairs, or nonfinite input

        integer :: n, i, first
        real(dp) :: global_scale
        logical, allocatable :: nonreal(:)
        real(dp), parameter :: real_noise_tolerance = 64.0_dp * epsilon(1.0_dp)

        wr = 0.0_dp
        vr = 0.0_dp
        status = linalg_status_invalid_shape
        n = size(w)
        if (size(v,1) /= n .or. size(v,2) /= n) return
        if (size(wr,1) /= n .or. size(wr,2) /= n) return
        if (size(vr,1) /= n .or. size(vr,2) /= n) return
        status = linalg_status_no_convergence
        if (.not. all(ieee_is_finite(real(w, dp)))) return
        if (.not. all(ieee_is_finite(aimag(w)))) return
        if (.not. all(ieee_is_finite(real(v, dp)))) return
        if (.not. all(ieee_is_finite(aimag(v)))) return

        ! Match nonreal eigenvalues in their occurrence order, permitting real
        ! eigenvalues between pair members. Interleaved conjugate pairs are
        ! intentionally rejected instead of returning a false decomposition.
        ! First preserve every exact nonreal pair, including tiny true pairs.
        allocate(nonreal(n))
        nonreal = aimag(w) /= 0.0_dp
        if (.not. valid_conjugate_pairs(w,nonreal)) then
            ! The self-contained SciFort eigensolver can sometimes return an
            ! effectively real eigenvalue with an imaginary rounding error.
            ! Only try the roundoff fallback when exact pairing fails.
            global_scale = tiny(1.0_dp)
            if (n > 0) global_scale = max(global_scale,maxval(abs(w)))
            do i = 1, n
                if (abs(aimag(w(i))) <= real_noise_tolerance * &
                    max(abs(w(i)),global_scale)) nonreal(i) = .false.
            end do
            status = linalg_status_invalid_shape
            if (.not. valid_conjugate_pairs(w,nonreal)) return
        end if

        do i = 1, n
            wr(i,i) = real(w(i),dp)
            if (.not. nonreal(i)) vr(:,i) = real(v(:,i),dp)
        end do
        first = 0
        do i = 1, n
            if (.not. nonreal(i)) cycle
            if (first == 0) then
                first = i
            else
                ! For v=a+ib, A*(-b) = alpha*(-b) + beta*a and
                ! A*a = -beta*(-b) + alpha*a. Only the first complex
                ! eigenvector is required: the conjugate eigenvector may
                ! carry an unrelated unit complex phase.
                vr(:,first) = -aimag(v(:,first))
                vr(:,i) = real(v(:,first),dp)
                wr(first,i) = aimag(w(first))
                wr(i,first) = aimag(w(i))
                first = 0
            end if
        end do
        status = linalg_status_success
    end subroutine cdf2rdf_complex

    pure logical function valid_conjugate_pairs(w, nonreal) result(ok)
        complex(dp), intent(in) :: w(:) !! finite eigenvalues
        logical, intent(in) :: nonreal(:) !! candidate members of conjugate pairs
        real(dp), parameter :: pair_tolerance = 128.0_dp * epsilon(1.0_dp)
        real(dp) :: scale
        integer :: i, first

        ok = .false.
        first = 0
        do i = 1, size(w)
            if (.not. nonreal(i)) cycle
            if (first == 0) then
                first = i
            else
                scale = max(abs(w(first)),abs(w(i)),tiny(1.0_dp))
                if (abs(w(first)/scale-conjg(w(i))/scale) > pair_tolerance) return
                first = 0
            end if
        end do
        ok = first == 0
    end function valid_conjugate_pairs

    pure subroutine cdf2rdf_complex_batch(w, v, wr, vr, status)
        complex(dp), intent(in) :: w(:, :) !! eigenvalues (batch,n), batch >=1
        complex(dp), intent(in) :: v(:, :, :) !! right eigenvectors (batch,n,n)
        real(dp), intent(out) :: wr(:, :, :) !! real block matrices (batch,n,n)
        real(dp), intent(out) :: vr(:, :, :) !! real basis matrices (batch,n,n)
        integer, intent(out) :: status !! shared success or first batch failure

        integer :: nb, n, b

        wr = 0.0_dp
        vr = 0.0_dp
        status = linalg_status_invalid_shape
        nb = size(w,1)
        n = size(w,2)
        if (nb == 0) return
        if (size(v,1) /= nb .or. size(v,2) /= n .or. size(v,3) /= n) return
        if (size(wr,1) /= nb .or. size(wr,2) /= n .or. size(wr,3) /= n) return
        if (size(vr,1) /= nb .or. size(vr,2) /= n .or. size(vr,3) /= n) return
        do b = 1, nb
            call cdf2rdf_complex(w(b,:), v(b,:,:), wr(b,:,:), vr(b,:,:), status)
            if (status /= linalg_status_success) then
                wr = 0.0_dp
                vr = 0.0_dp
                return
            end if
        end do
        status = linalg_status_success
    end subroutine cdf2rdf_complex_batch

    pure subroutine cdf2rdf_real(w, v, wr, vr, status)
        real(dp), intent(in) :: w(:) !! real eigenvalues (n)
        real(dp), intent(in) :: v(:, :) !! real right eigenvectors (n,n)
        real(dp), intent(out) :: wr(:, :) !! diagonal eigenvalue matrix (n,n)
        real(dp), intent(out) :: vr(:, :) !! unchanged real eigenvectors (n,n)
        integer, intent(out) :: status !! status of shape and finite checks

        integer :: n, i

        wr = 0.0_dp
        vr = 0.0_dp
        status = linalg_status_invalid_shape
        n = size(w)
        if (size(v,1) /= n .or. size(v,2) /= n) return
        if (size(wr,1) /= n .or. size(wr,2) /= n) return
        if (size(vr,1) /= n .or. size(vr,2) /= n) return
        status = linalg_status_no_convergence
        if (.not. all(ieee_is_finite(w))) return
        if (.not. all(ieee_is_finite(v))) return
        do i = 1, n
            wr(i,i) = w(i)
        end do
        vr = v
        status = linalg_status_success
    end subroutine cdf2rdf_real

    pure subroutine cdf2rdf_real_batch(w, v, wr, vr, status)
        real(dp), intent(in) :: w(:, :) !! eigenvalues (batch,n), batch >= 1
        real(dp), intent(in) :: v(:, :, :) !! eigenvectors (batch,n,n)
        real(dp), intent(out) :: wr(:, :, :) !! diagonal matrices (batch,n,n)
        real(dp), intent(out) :: vr(:, :, :) !! unchanged eigenvectors (batch,n,n)
        integer, intent(out) :: status !! success, shape error or nonfinite input

        integer :: nb, n, b

        wr = 0.0_dp
        vr = 0.0_dp
        status = linalg_status_invalid_shape
        nb = size(w,1)
        n = size(w,2)
        if (nb == 0) return
        if (size(v,1) /= nb .or. size(v,2) /= n .or. size(v,3) /= n) return
        if (size(wr,1) /= nb .or. size(wr,2) /= n .or. size(wr,3) /= n) return
        if (size(vr,1) /= nb .or. size(vr,2) /= n .or. size(vr,3) /= n) return
        do b = 1, nb
            call cdf2rdf_real(w(b,:), v(b,:,:), wr(b,:,:), vr(b,:,:), status)
            if (status /= linalg_status_success) then
                wr = 0.0_dp
                vr = 0.0_dp
                return
            end if
        end do
        status = linalg_status_success
    end subroutine cdf2rdf_real_batch

    pure subroutine cdf2rdf_real_values(w, v, wr, vr, status)
        real(dp), intent(in) :: w(:) !! real eigenvalues (n)
        complex(dp), intent(in) :: v(:, :) !! corresponding complex right eigenvectors (n,n)
        real(dp), intent(out) :: wr(:, :) !! diagonal eigenvalue matrix (n,n)
        real(dp), intent(out) :: vr(:, :) !! real parts of the eigenvectors (n,n)
        integer, intent(out) :: status !! success, shape error or nonfinite input

        call cdf2rdf_complex(cmplx(w,0.0_dp,dp), v, wr, vr, status)
    end subroutine cdf2rdf_real_values

    pure subroutine cdf2rdf_real_values_batch(w, v, wr, vr, status)
        real(dp), intent(in) :: w(:, :) !! real eigenvalues (batch,n)
        complex(dp), intent(in) :: v(:, :, :) !! right eigenvectors (batch,n,n)
        real(dp), intent(out) :: wr(:, :, :) !! diagonal block matrices (batch,n,n)
        real(dp), intent(out) :: vr(:, :, :) !! real eigenvectors (batch,n,n)
        integer, intent(out) :: status !! success, shape error or nonfinite input

        call cdf2rdf_complex_batch(cmplx(w,0.0_dp,dp), v, wr, vr, status)
    end subroutine cdf2rdf_real_values_batch

end module scifort_linalg_cdf2rdf
