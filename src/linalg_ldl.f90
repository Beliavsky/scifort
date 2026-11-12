! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Symmetric/Hermitian indefinite LDL factorization with Bunch-Kaufman pivoting.
! Independently implemented from the diagonal-pivoting mathematical algorithm
! described in the LAPACK DSYTF2 documentation (not translated LAPACK code).
! https://www.netlib.org/lapack/explore-html/db/d8d/
! group__hetf2_ga91a9fdb3b98a38e1a91689b65c366569.html
module scifort_linalg_ldl
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_singular, linalg_status_no_convergence
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    implicit none
    private

    public :: ldl

    interface ldl
        module procedure ldl_real
        module procedure ldl_complex
    end interface ldl

contains

    subroutine ldl_real(a, lu, d, perm, status, lower)
        real(dp), intent(in) :: a(:, :) !! symmetric n-by-n matrix; selected triangle used
        real(dp), intent(out) :: lu(:, :) !! n-by-n outer factor; lu(perm,:) is triangular
        real(dp), intent(out) :: d(:, :) !! block-diagonal n-by-n factor, blocks of order 1 or 2
        integer, intent(out) :: perm(:) !! one-based row indices for triangular lu(perm,:)
        integer, intent(out) :: status !! linalg_status_* result
        logical, intent(in), optional :: lower !! use lower triangle (default true), else upper
        complex(dp), allocatable :: luc(:, :), dc(:, :)
        integer :: n

        lu = 0.0_dp
        d = 0.0_dp
        perm = 0
        status = linalg_status_invalid_shape
        n = size(a, 1)
        if (size(a, 2) /= n) return
        if (size(lu, 1) /= n .or. size(lu, 2) /= n) return
        if (size(d, 1) /= n .or. size(d, 2) /= n .or. size(perm) /= n) return

        allocate(luc(n,n), dc(n,n))
        call ldl_complex(cmplx(a,kind=dp), luc, dc, perm, status, lower=lower)
        if (status /= linalg_status_success) then
            perm = 0
            return
        end if
        lu = real(luc,dp)
        d = real(dc,dp)
    end subroutine ldl_real

    subroutine ldl_complex(a, lu, d, perm, status, lower, hermitian)
        complex(dp), intent(in) :: a(:, :) !! symmetric or Hermitian n-by-n matrix
        complex(dp), intent(out) :: lu(:, :) !! outer factor; lu(perm,:) is unit triangular
        complex(dp), intent(out) :: d(:, :) !! block-diagonal factor, at most 2-by-2 blocks
        integer, intent(out) :: perm(:) !! one-based permutation that triangularizes lu
        integer, intent(out) :: status !! linalg_status_* result
        logical, intent(in), optional :: lower !! use input lower triangle (default true)
        logical, intent(in), optional :: hermitian !! use conjugate transpose (default true)

        ! Work in lower form. The reverse-identity transform maps an upper
        ! factorization to a lower one without reading the unused input triangle.
        real(dp), parameter :: alpha = (1.0_dp + sqrt(17.0_dp))/8.0_dp
        complex(dp), allocatable :: w(:, :), l(:, :), blocks(:, :)
        complex(dp) :: diag, det, l1, l2, p11, p12, p21, p22, v1, v2
        complex(dp) :: entry, term, tmp
        real(dp) :: absakk, colmax, rowmax, scale, tolerance
        integer, allocatable :: p(:)
        integer :: i, imax, j, k, kp, n, pivot_size, pivot_row, item
        logical :: use_lower, use_hermitian

        lu = cmplx(0.0_dp, 0.0_dp, dp)
        d = cmplx(0.0_dp, 0.0_dp, dp)
        perm = 0
        status = linalg_status_invalid_shape
        n = size(a, 1)
        if (size(a, 2) /= n .or. size(lu, 1) /= n .or. size(lu, 2) /= n) return
        if (size(d, 1) /= n .or. size(d, 2) /= n .or. size(perm) /= n) return
        if (n == 0) then
            status = linalg_status_success
            return
        end if
        use_lower = .true.
        if (present(lower)) use_lower = lower
        use_hermitian = .true.
        if (present(hermitian)) use_hermitian = hermitian

        allocate(w(n,n),l(n,n),blocks(n,n),p(n))
        w = cmplx(0.0_dp, 0.0_dp, dp)
        l = cmplx(0.0_dp, 0.0_dp, dp)
        blocks = cmplx(0.0_dp, 0.0_dp, dp)
        do i = 1, n
            p(i) = i
            l(i,i) = cmplx(1.0_dp, 0.0_dp, dp)
            do j = 1, i
                if (use_lower) then
                    entry = a(i,j)
                else
                    entry = a(n+1-i,n+1-j)
                end if
                if (.not. finite_complex(entry)) then
                    status = linalg_status_no_convergence
                    return
                end if
                if (use_hermitian .and. i == j) then
                    ! A Hermitian matrix has a real diagonal. Do not silently
                    ! discard an imaginary diagonal supplied by the caller.
                    if (abs(aimag(entry)) > 0.0_dp) return
                    entry = cmplx(real(entry,dp),0.0_dp,dp)
                end if
                w(i,j) = entry
                if (use_hermitian) then
                    w(j,i) = conjg(entry)
                else
                    w(j,i) = entry
                end if
            end do
        end do

        status = linalg_status_no_convergence
        k = 1
        do while (k <= n)
            pivot_size = 1
            kp = k
            absakk = abs(w(k,k))
            colmax = 0.0_dp
            imax = k
            do i = k+1, n
                if (abs(w(i,k)) > colmax) then
                    colmax = abs(w(i,k))
                    imax = i
                end if
            end do

            if (absakk < alpha*colmax) then
                rowmax = 0.0_dp
                do j = k, n
                    if (j /= imax) rowmax = max(rowmax,abs(w(imax,j)))
                end do
                if (abs(w(imax,imax)) >= alpha*rowmax) then
                    kp = imax
                else if (rowmax > 0.0_dp) then
                    if (absakk/colmax >= alpha*(colmax/rowmax)) then
                        kp = k
                    else
                        kp = imax
                        pivot_size = 2
                    end if
                else
                    ! A nonzero colmax implies rowmax > 0, so this is
                    ! reachable only through underflow in the comparisons.
                    status = linalg_status_no_convergence
                    return
                end if
            end if

            pivot_row = k
            if (pivot_size == 2) pivot_row = k+1
            if (kp /= pivot_row) then
                ! Symmetric row/column interchange in the active Schur
                ! complement, and the corresponding rows of earlier L.
                do j = 1, n
                    tmp = w(pivot_row,j)
                    w(pivot_row,j) = w(kp,j)
                    w(kp,j) = tmp
                end do
                do i = 1, n
                    tmp = w(i,pivot_row)
                    w(i,pivot_row) = w(i,kp)
                    w(i,kp) = tmp
                end do
                if (k > 1) then
                    do j = 1, k-1
                        tmp = l(pivot_row,j)
                        l(pivot_row,j) = l(kp,j)
                        l(kp,j) = tmp
                    end do
                end if
                item = p(pivot_row)
                p(pivot_row) = p(kp)
                p(kp) = item
            end if

            if (pivot_size == 1) then
                diag = w(k,k)
                blocks(k,k) = diag
                if (abs(diag) > 0.0_dp) then
                    do i = k+1, n
                        l(i,k) = w(i,k)/diag
                    end do
                    do j = k+1, n
                        do i = j, n
                            term = w(j,k)
                            if (use_hermitian) term = conjg(term)
                            w(i,j) = w(i,j) - l(i,k)*term
                            call set_mirror(w,i,j,use_hermitian)
                        end do
                    end do
                else if (colmax > 0.0_dp) then
                    status = linalg_status_singular
                    return
                end if
                k = k+1
            else
                ! The pivot is D = [p11 p12; p21 p22]. For each trailing
                ! row v, solve L_row * D = v without assuming D is diagonal.
                p11 = w(k,k)
                p12 = w(k,k+1)
                p21 = w(k+1,k)
                p22 = w(k+1,k+1)
                blocks(k,k) = p11
                blocks(k,k+1) = p12
                blocks(k+1,k) = p21
                blocks(k+1,k+1) = p22
                scale = max(abs(p11),abs(p12),abs(p21),abs(p22))
                if (scale <= 0.0_dp) then
                    status = linalg_status_singular
                    return
                end if
                p11 = p11/scale
                p12 = p12/scale
                p21 = p21/scale
                p22 = p22/scale
                det = p11*p22 - p12*p21
                tolerance = 8.0_dp*epsilon(1.0_dp)* &
                    max(abs(p11*p22),abs(p12*p21))
                if (abs(det) <= tolerance) then
                    status = linalg_status_singular
                    return
                end if
                do i = k+2, n
                    v1 = w(i,k)/scale
                    v2 = w(i,k+1)/scale
                    l1 = (v1*p22-v2*p21)/det
                    l2 = (v2*p11-v1*p12)/det
                    l(i,k) = l1
                    l(i,k+1) = l2
                end do
                do j = k+2, n
                    do i = j, n
                        term = w(j,k)
                        tmp = w(j,k+1)
                        if (use_hermitian) then
                            term = conjg(term)
                            tmp = conjg(tmp)
                        end if
                        w(i,j) = w(i,j)-l(i,k)*term-l(i,k+1)*tmp
                        call set_mirror(w,i,j,use_hermitian)
                    end do
                end do
                k = k+2
            end if
            if (.not. finite_matrix(w)) then
                status = linalg_status_no_convergence
                return
            end if
        end do

        ! This factor is intentionally nontriangular until the caller
        ! applies perm: A = lu * D * lu**H (or lu * D * lu**T).
        do i = 1, n
            lu(p(i),:) = l(i,:)
        end do
        d = blocks
        perm = p
        if (.not. use_lower) then
            lu = lu(n:1:-1,n:1:-1)
            d = d(n:1:-1,n:1:-1)
            perm = n+1-perm(n:1:-1)
        end if
        status = linalg_status_success
    end subroutine ldl_complex

    pure subroutine set_mirror(w, i, j, hermitian)
        complex(dp), intent(inout) :: w(:, :) !! updated symmetric/Hermitian matrix
        integer, intent(in) :: i !! active row
        integer, intent(in) :: j !! active column
        logical, intent(in) :: hermitian !! conjugate when filling upper triangle

        if (hermitian) then
            if (i == j) w(i,j) = cmplx(real(w(i,j),dp),0.0_dp,dp)
            w(j,i) = conjg(w(i,j))
        else
            w(j,i) = w(i,j)
        end if
    end subroutine set_mirror

    pure logical function finite_complex(z) result(ok)
        complex(dp), intent(in) :: z !! complex scalar to check
        ok = ieee_is_finite(real(z,dp)) .and. ieee_is_finite(aimag(z))
    end function finite_complex

    pure logical function finite_matrix(z) result(ok)
        complex(dp), intent(in) :: z(:, :) !! complex matrix to check
        ok = all(ieee_is_finite(real(z,dp))) .and. all(ieee_is_finite(aimag(z)))
    end function finite_matrix

end module scifort_linalg_ldl
