! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Orthonormal Helmert contrasts, corresponding to scipy.linalg.helmert.
module scifort_linalg_helmert
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape
    implicit none
    private

    public :: helmert

contains

    pure subroutine helmert(n, h, status, full)
        integer, intent(in) :: n !! order of the Helmert matrix; must be >= 1
        real(dp), intent(out) :: h(:,:) !! shape (n,n) if full, otherwise (n-1,n)
        integer, intent(out) :: status !! success or invalid_shape
        logical, intent(in), optional :: full !! include the normalized all-ones first row; default false
        logical :: include_mean
        integer :: k, shift
        real(dp) :: normalization

        h = 0.0_dp
        status = linalg_status_invalid_shape
        include_mean = .false.
        if (present(full)) include_mean = full
        if (n < 1) return
        shift = 1
        if (include_mean) shift = 0
        if (size(h,1) /= n-shift .or. size(h,2) /= n) return

        if (include_mean) h(1,:) = 1.0_dp/sqrt(real(n,dp))
        do k=2,n
            ! Row k has (k-1) equal positive coefficients and one negative
            ! coefficient; using real arithmetic avoids integer overflow.
            normalization = sqrt(real(k-1,dp)*real(k,dp))
            h(k-shift,1:k-1) = 1.0_dp/normalization
            h(k-shift,k) = -real(k-1,dp)/normalization
        end do
        status = linalg_status_success
    end subroutine helmert

end module scifort_linalg_helmert
