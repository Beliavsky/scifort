! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Upper Hessenberg reduction by two-sided Householder similarities.
! The implementation is independent and based on standard Householder theory.
module scifort_linalg_hessenberg
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    implicit none
    private

    public :: hessenberg

    interface hessenberg
        module procedure hessenberg_real
        module procedure hessenberg_complex
    end interface hessenberg

contains

    subroutine hessenberg_real(a, h, status, q)
        real(dp), intent(in) :: a(:, :) !! finite square input matrix (n,n), including n=0
        real(dp), intent(out) :: h(:, :) !! upper Hessenberg matrix (n,n), A = Q H Q^T
        integer, intent(out) :: status !! success, invalid shape, or numerical failure
        real(dp), intent(out), optional :: q(:, :) !! optional orthogonal similarity matrix (n,n)

        integer :: i, j, k, n, m
        real(dp) :: amax, norm_x, norm_v, phase, xscale, projection
        real(dp), allocatable :: v(:)

        status = linalg_status_invalid_shape
        h = 0.0_dp
        if (present(q)) q = 0.0_dp
        n = size(a, 1)
        if (size(a, 2) /= n .or. size(h, 1) /= n .or. size(h, 2) /= n) return
        if (present(q)) then
            if (size(q, 1) /= n .or. size(q, 2) /= n) return
        end if
        status = linalg_status_no_convergence
        if (.not. all(ieee_is_finite(a))) return
        if (present(q)) then
            do i = 1, n
                q(i, i) = 1.0_dp
            end do
        end if
        if (n == 0) then
            status = linalg_status_success
            return
        end if

        ! Global normalization prevents overflow in the two-sided updates.
        amax = maxval(abs(a))
        if (amax == 0.0_dp) amax = 1.0_dp
        h = a / amax
        do k = 1, n - 2
            m = n - k
            allocate(v(m))
            xscale = maxval(abs(h(k + 1:n, k)))
            if (xscale /= 0.0_dp) then
                v = h(k + 1:n, k) / xscale
                norm_x = sqrt(sum(v * v))
                phase = 1.0_dp
                if (v(1) < 0.0_dp) phase = -1.0_dp
                v(1) = v(1) + phase * norm_x
                norm_v = sqrt(sum(v * v))
                if (norm_v == 0.0_dp) return
                v = v / norm_v

                ! H = (I - 2 vv^T) H (I - 2 vv^T).
                do j = k, n
                    projection = 2.0_dp * dot_product(v, h(k + 1:n, j))
                    h(k + 1:n, j) = h(k + 1:n, j) - projection * v
                end do
                do i = 1, n
                    projection = 2.0_dp * dot_product(h(i, k + 1:n), v)
                    h(i, k + 1:n) = h(i, k + 1:n) - projection * v
                end do
                if (present(q)) then
                    do i = 1, n
                        projection = 2.0_dp * dot_product(q(i, k + 1:n), v)
                        q(i, k + 1:n) = q(i, k + 1:n) - projection * v
                    end do
                end if
            end if
            h(k + 2:n, k) = 0.0_dp
            deallocate(v)
        end do
        h = h * amax
        if (.not. all(ieee_is_finite(h))) return
        if (present(q)) then
            if (.not. all(ieee_is_finite(q))) return
        end if
        status = linalg_status_success
    end subroutine hessenberg_real

    subroutine hessenberg_complex(a, h, status, q)
        complex(dp), intent(in) :: a(:, :) !! finite square input matrix (n,n), including n=0
        complex(dp), intent(out) :: h(:, :) !! upper Hessenberg matrix (n,n), A = Q H Q^H
        integer, intent(out) :: status !! success, invalid shape, or numerical failure
        complex(dp), intent(out), optional :: q(:, :) !! optional unitary similarity matrix (n,n)

        integer :: i, j, k, n, m
        real(dp) :: amax, norm_x, norm_v, xscale
        complex(dp) :: phase, projection
        complex(dp), allocatable :: v(:)

        status = linalg_status_invalid_shape
        h = cmplx(0.0_dp, 0.0_dp, dp)
        if (present(q)) q = cmplx(0.0_dp, 0.0_dp, dp)
        n = size(a, 1)
        if (size(a, 2) /= n .or. size(h, 1) /= n .or. size(h, 2) /= n) return
        if (present(q)) then
            if (size(q, 1) /= n .or. size(q, 2) /= n) return
        end if
        status = linalg_status_no_convergence
        if (.not. all(ieee_is_finite(real(a, dp))) .or. &
            .not. all(ieee_is_finite(aimag(a)))) return
        if (present(q)) then
            do i = 1, n
                q(i, i) = cmplx(1.0_dp, 0.0_dp, dp)
            end do
        end if
        if (n == 0) then
            status = linalg_status_success
            return
        end if

        ! Componentwise maximum avoids an overflowing complex magnitude.
        amax = max(maxval(abs(real(a, dp))), maxval(abs(aimag(a))))
        if (amax == 0.0_dp) amax = 1.0_dp
        h = a / amax
        do k = 1, n - 2
            m = n - k
            allocate(v(m))
            xscale = max(maxval(abs(real(h(k + 1:n, k), dp))), &
                maxval(abs(aimag(h(k + 1:n, k)))))
            if (xscale /= 0.0_dp) then
                v = h(k + 1:n, k) / xscale
                norm_x = sqrt(sum(real(v, dp)**2 + aimag(v)**2))
                phase = cmplx(1.0_dp, 0.0_dp, dp)
                if (abs(v(1)) /= 0.0_dp) phase = v(1) / abs(v(1))
                v(1) = v(1) + phase * norm_x
                norm_v = sqrt(sum(real(v, dp)**2 + aimag(v)**2))
                if (norm_v == 0.0_dp) return
                v = v / norm_v

                ! Householder reflector P = I - 2 vv^H is unitary and Hermitian.
                do j = k, n
                    projection = 2.0_dp * dot_product(v, h(k + 1:n, j))
                    h(k + 1:n, j) = h(k + 1:n, j) - projection * v
                end do
                do i = 1, n
                    projection = 2.0_dp * sum(h(i, k + 1:n) * v)
                    h(i, k + 1:n) = h(i, k + 1:n) - projection * conjg(v)
                end do
                if (present(q)) then
                    do i = 1, n
                        projection = 2.0_dp * sum(q(i, k + 1:n) * v)
                        q(i, k + 1:n) = q(i, k + 1:n) - projection * conjg(v)
                    end do
                end if
            end if
            h(k + 2:n, k) = cmplx(0.0_dp, 0.0_dp, dp)
            deallocate(v)
        end do
        h = h * amax
        if (.not. all(ieee_is_finite(real(h, dp))) .or. &
            .not. all(ieee_is_finite(aimag(h)))) return
        if (present(q)) then
            if (.not. all(ieee_is_finite(real(q, dp))) .or. &
                .not. all(ieee_is_finite(aimag(q)))) return
        end if
        status = linalg_status_success
    end subroutine hessenberg_complex

end module scifort_linalg_hessenberg
