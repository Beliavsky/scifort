! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Givens bulge chasing for symmetric/Hermitian band-to-tridiagonal reduction.
! Each arithmetic kind holds its lower band plus one transient bulge diagonal.
! There is no dense matrix allocation unless the transformation Q is requested.
module scifort_linalg_banded_tridiagonal
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    implicit none
    private
    public :: band_to_tridiagonal

    interface band_to_tridiagonal
        module procedure band_to_tridiagonal_real
        module procedure band_to_tridiagonal_complex
    end interface band_to_tridiagonal

contains

    subroutine band_to_tridiagonal_real(ab, d, e, status, q, lower)
        real(dp), intent(in) :: ab(:, :) !! SciPy upper/lower packed symmetric matrix, (kd+1,n)
        real(dp), intent(out) :: d(:) !! tridiagonal diagonal, shape (n)
        real(dp), intent(out) :: e(:) !! nonnegative tridiagonal offdiagonal, shape (n-1)
        integer, intent(out) :: status !! SciFort linalg status
        real(dp), intent(out), optional :: q(:, :) !! orthogonal similarity matrix, shape (n,n)
        logical, intent(in), optional :: lower !! lower-packed storage, default false

        real(dp), allocatable :: work(:, :), phases(:)
        real(dp) :: entry, a, b, c, s, tmp, divisor, r, sigma, off
        integer :: n, ld, kd, j, k, i, col, row, p, high, st
        logical :: low

        status = linalg_status_invalid_shape
        d = 0.0_dp
        e = 0.0_dp
        if (present(q)) q = 0.0_dp
        n = size(ab,2)
        ld = size(ab,1)-1
        if (ld < 0) return
        if (size(d) /= n .or. size(e) /= max(0,n-1)) return
        if (present(q)) then
            if (size(q,1) /= n .or. size(q,2) /= n) return
        end if
        low = .false.
        if (present(lower)) low = lower
        status = linalg_status_no_convergence
        kd = min(ld,max(0,n-1))
        allocate(work(0:kd+1,n),phases(n))
        work = 0.0_dp
        sigma = 0.0_dp
        do j = 1,n
            high = min(kd,n-j)
            do k = 0,high
                if (low) then
                    entry = ab(k+1,j)
                else
                    entry = ab(ld+1-k,j+k)
                end if
                if (.not. ieee_is_finite(entry)) return
                sigma = max(sigma,abs(entry))
                work(k,j) = entry
            end do
        end do
        if (present(q)) then
            do j = 1,n
                q(j,j) = 1.0_dp
            end do
        end if
        if (n == 0 .or. sigma == 0.0_dp) then
            status = linalg_status_success
            return
        end if
        do j = 1,n
            do k = 0,min(kd,n-j)
                work(k,j) = work(k,j)/sigma
            end do
        end do

        ! Eliminate each unwanted offdiagonal by an orthogonal Givens
        ! similarity, chasing the resulting bulge through the extra band.
        if (kd >= 2) then
            do i = 1,n-2
                do j = min(n,i+kd),i+2,-1
                    row = i
                    col = j
                    do while (col <= n)
                        p = col-1
                        a = get_band_real(work,row,p)
                        b = get_band_real(work,row,col)
                        if (b /= 0.0_dp) then
                            divisor = max(abs(a),abs(b))
                            if (divisor == 0.0_dp) exit
                            r = hypot(a/divisor,b/divisor)
                            c = (a/divisor)/r
                            s = (b/divisor)/r
                            call apply_rotation_real(work,p,c,s,st)
                            if (st /= linalg_status_success) return
                            call put_band_real(work,row,col,0.0_dp,st)
                            if (st /= linalg_status_success) return
                            if (present(q)) then
                                do k = 1,n
                                    tmp = q(k,p)
                                    q(k,p) = c*tmp+s*q(k,col)
                                    q(k,col) = -s*tmp+c*q(k,col)
                                end do
                            end if
                        end if
                        row = col-1
                        col = col+kd
                    end do
                end do
            end do
        end if
        do i = 1,n
            d(i) = work(0,i)*sigma
            if (.not. ieee_is_finite(d(i))) return
        end do
        if (n > 0) phases(1) = 1.0_dp
        do i = 1,n-1
            off = work(1,i)
            e(i) = abs(off)*sigma
            if (.not. ieee_is_finite(e(i))) return
            phases(i+1) = phases(i)
            if (off < 0.0_dp) phases(i+1) = -phases(i)
        end do
        if (present(q)) then
            do j = 1,n
                q(:,j) = q(:,j)*phases(j)
            end do
        end if
        status = linalg_status_success
    end subroutine band_to_tridiagonal_real

    subroutine band_to_tridiagonal_complex(ab, d, e, status, q, lower)
        complex(dp), intent(in) :: ab(:, :) !! SciPy-packed Hermitian bands, shape (kd+1,n)
        real(dp), intent(out) :: d(:) !! real diagonal of reduced tridiagonal, shape (n)
        real(dp), intent(out) :: e(:) !! nonnegative offdiagonal, shape (n-1)
        integer, intent(out) :: status !! SciFort linalg status
        complex(dp), intent(out), optional :: q(:, :) !! unitary similarity factor, shape (n,n)
        logical, intent(in), optional :: lower !! use lower packed storage, default false

        complex(dp), allocatable :: work(:, :), phases(:)
        complex(dp) :: entry, a, b, c, s, tmp, off
        real(dp) :: sigma, mag, divisor, r
        integer :: n, ld, kd, j, k, i, col, row, p, high, st
        logical :: low

        status = linalg_status_invalid_shape
        d = 0.0_dp
        e = 0.0_dp
        if (present(q)) q = cmplx(0.0_dp,0.0_dp,dp)
        n = size(ab,2)
        ld = size(ab,1)-1
        if (ld < 0) return
        if (size(d) /= n .or. size(e) /= max(0,n-1)) return
        if (present(q)) then
            if (size(q,1) /= n .or. size(q,2) /= n) return
        end if
        low = .false.
        if (present(lower)) low = lower
        status = linalg_status_no_convergence
        kd = min(ld,max(0,n-1))
        allocate(work(0:kd+1,n),phases(n))
        work = cmplx(0.0_dp,0.0_dp,dp)
        sigma = 0.0_dp
        do j = 1,n
            high = min(kd,n-j)
            do k = 0,high
                if (low) then
                    entry = ab(k+1,j)
                else
                    entry = conjg(ab(ld+1-k,j+k))
                end if
                if (k == 0) entry = cmplx(real(entry,dp),0.0_dp,dp)
                if (.not. finite_complex(entry)) return
                sigma = max(sigma,abs(real(entry,dp)),abs(aimag(entry)))
                work(k,j) = entry
            end do
        end do
        if (present(q)) then
            do j = 1,n
                q(j,j) = cmplx(1.0_dp,0.0_dp,dp)
            end do
        end if
        if (n == 0 .or. sigma == 0.0_dp) then
            status = linalg_status_success
            return
        end if
        do j = 1,n
            do k = 0,min(kd,n-j)
                entry = work(k,j)
                work(k,j) = cmplx(real(entry,dp)/sigma,aimag(entry)/sigma,dp)
            end do
        end do

        ! Annihilate row i beyond its first superdiagonal, farthest first.
        ! A Givens similarity (p,p+1) sends one bulge to distance kd+1;
        ! each successive rotation chases it kd columns to the right.
        ! The active bandwidth is therefore kd+1 throughout.
        if (kd >= 2) then
            do i = 1,n-2
                do j = min(n,i+kd),i+2,-1
                    row = i
                    col = j
                    do while (col <= n)
                        p = col-1
                        a = get_band(work,row,p)
                        b = get_band(work,row,col)
                        if (b /= cmplx(0.0_dp,0.0_dp,dp)) then
                            ! Form complex SU(2) Givens G with [a,b]*G=[r,0].
                            ! Compute ratios first to avoid a*a+b*b overflow.
                            divisor = max(abs(a),abs(b))
                            if (divisor == 0.0_dp) exit
                            r = hypot(abs(a/divisor),abs(b/divisor))
                            c = conjg(a/divisor)/r
                            s = conjg(b/divisor)/r
                            call apply_rotation(work,p,c,s,st)
                            if (st /= linalg_status_success) return
                            call put_band(work,row,col,cmplx(0.0_dp,0.0_dp,dp),st)
                            if (st /= linalg_status_success) return
                            if (present(q)) then
                                do k = 1,n
                                    tmp = q(k,p)
                                    q(k,p) = c*tmp+s*q(k,col)
                                    q(k,col) = -conjg(s)*tmp+conjg(c)*q(k,col)
                                end do
                            end if
                        end if
                        row = col-1
                        col = col+kd
                    end do
                end do
            end do
        end if
        do i = 1,n
            d(i) = real(work(0,i),dp)*sigma
            if (.not. ieee_is_finite(d(i))) return
        end do
        if (n > 0) phases(1) = cmplx(1.0_dp,0.0_dp,dp)
        do i = 1,n-1
            off = work(1,i)
            mag = abs(off)
            e(i) = mag*sigma
            if (.not. ieee_is_finite(e(i))) return
            if (mag == 0.0_dp) then
                phases(i+1) = phases(i)
            else
                phases(i+1) = phases(i)*(off/mag)
            end if
        end do
        if (present(q)) then
            do j = 1,n
                q(:,j) = q(:,j)*phases(j)
            end do
        end if
        status = linalg_status_success
    end subroutine band_to_tridiagonal_complex

    pure function get_band_real(work,i,j) result(value)
        real(dp), intent(in) :: work(0:,:) !! lower symmetric band including bulge
        integer, intent(in) :: i !! row index
        integer, intent(in) :: j !! column index
        real(dp) :: value
        integer :: distance

        value = 0.0_dp
        distance = abs(i-j)
        if (distance > ubound(work,1)) return
        if (i >= j) then
            value = work(distance,j)
        else
            value = work(distance,i)
        end if
    end function get_band_real

    subroutine put_band_real(work,i,j,value,status)
        real(dp), intent(inout) :: work(0:,:) !! lower symmetric band including bulge
        integer, intent(in) :: i !! row index
        integer, intent(in) :: j !! column index
        real(dp), intent(in) :: value !! symmetric matrix entry (i,j)
        integer, intent(out) :: status !! success or excessive discarded fill
        integer :: distance

        status = linalg_status_success
        distance = abs(i-j)
        if (distance > ubound(work,1)) then
            if (abs(value) > 2.0e-12_dp) status = linalg_status_no_convergence
            return
        end if
        if (i >= j) then
            work(distance,j) = value
        else
            work(distance,i) = value
        end if
    end subroutine put_band_real

    subroutine apply_rotation_real(work,p,c,s,status)
        real(dp), intent(inout) :: work(0:,:) !! lower symmetric band plus transient bulge
        integer, intent(in) :: p !! first rotated column, followed by p+1
        real(dp), intent(in) :: c !! cosine parameter of the Givens rotation
        real(dp), intent(in) :: s !! sine parameter of the Givens rotation
        integer, intent(out) :: status !! success or unexpected out-of-band fill

        real(dp) :: oldp, oldq, apq, pp, qq, qp, tmp1, tmp2, app, aqq
        integer :: k, q, n, width, st

        status = linalg_status_success
        q = p+1
        n = size(work,2)
        width = ubound(work,1)
        app = get_band_real(work,p,p)
        aqq = get_band_real(work,q,q)
        apq = get_band_real(work,p,q)
        do k = max(1,p-width),min(n,q+width)
            if (k == p .or. k == q) cycle
            oldp = get_band_real(work,k,p)
            oldq = get_band_real(work,k,q)
            call put_band_real(work,k,p,c*oldp+s*oldq,st)
            if (st /= linalg_status_success) then
                status = st
                return
            end if
            call put_band_real(work,k,q,-s*oldp+c*oldq,st)
            if (st /= linalg_status_success) then
                status = st
                return
            end if
        end do
        tmp1 = app*c+apq*s
        tmp2 = apq*c+aqq*s
        pp = c*tmp1+s*tmp2
        qp = -s*tmp1+c*tmp2
        tmp1 = -app*s+apq*c
        tmp2 = -apq*s+aqq*c
        qq = -s*tmp1+c*tmp2
        work(0,p) = pp
        work(0,q) = qq
        work(1,p) = qp
    end subroutine apply_rotation_real

    pure function get_band(work,i,j) result(value)
        complex(dp), intent(in) :: work(0:,:) !! lower Hermitian band including bulge row
        integer, intent(in) :: i !! row index
        integer, intent(in) :: j !! column index
        complex(dp) :: value
        integer :: distance

        value = cmplx(0.0_dp,0.0_dp,dp)
        distance = abs(i-j)
        if (distance > ubound(work,1)) return
        if (i >= j) then
            value = work(distance,j)
        else
            value = conjg(work(distance,i))
        end if
    end function get_band

    subroutine put_band(work,i,j,value,status)
        complex(dp), intent(inout) :: work(0:,:) !! lower Hermitian band including bulge row
        integer, intent(in) :: i !! row index
        integer, intent(in) :: j !! column index
        complex(dp), intent(in) :: value !! matrix entry A(i,j)
        integer, intent(out) :: status !! success if entry fits active band
        integer :: distance

        status = linalg_status_success
        distance = abs(i-j)
        if (distance > ubound(work,1)) then
            ! The discarded entries must be zero to rounding error; a larger
            ! value signals an unexpected fill-in and must not be hidden.
            if (abs(value) > 2.0e-12_dp) status = linalg_status_no_convergence
            return
        end if
        if (i >= j) then
            work(distance,j) = value
        else
            work(distance,i) = conjg(value)
        end if
    end subroutine put_band

    subroutine apply_rotation(work,p,c,s,status)
        complex(dp), intent(inout) :: work(0:,:) !! lower Hermitian band plus transient bulge
        integer, intent(in) :: p !! first rotated column, second is p+1
        complex(dp), intent(in) :: c !! G(1,1)=c, G(2,2)=conjg(c)
        complex(dp), intent(in) :: s !! G(2,1)=s, G(1,2)=-conjg(s)
        integer, intent(out) :: status !! success or unexpected out-of-band fill

        complex(dp) :: oldp, oldq, apq, pp, qq, qp, tmp1, tmp2
        real(dp) :: app, aqq
        integer :: k, q, n, width, st

        status = linalg_status_success
        q = p+1
        n = size(work,2)
        width = ubound(work,1)
        app = real(get_band(work,p,p),dp)
        aqq = real(get_band(work,q,q),dp)
        apq = get_band(work,p,q)
        do k = max(1,p-width),min(n,q+width)
            if (k == p .or. k == q) cycle
            oldp = get_band(work,k,p)
            oldq = get_band(work,k,q)
            call put_band(work,k,p,c*oldp+s*oldq,st)
            if (st /= linalg_status_success) then
                status = st
                return
            end if
            call put_band(work,k,q,-conjg(s)*oldp+conjg(c)*oldq,st)
            if (st /= linalg_status_success) then
                status = st
                return
            end if
        end do
        ! Two-sided unitary update to the 2x2 Hermitian pivot block.
        tmp1 = app*c+apq*s
        tmp2 = conjg(apq)*c+aqq*s
        pp = conjg(c)*tmp1+conjg(s)*tmp2
        qp = -s*tmp1+c*tmp2
        tmp1 = -app*conjg(s)+apq*conjg(c)
        tmp2 = -conjg(apq)*conjg(s)+aqq*conjg(c)
        qq = -s*tmp1+c*tmp2
        work(0,p) = cmplx(real(pp,dp),0.0_dp,dp)
        work(0,q) = cmplx(real(qq,dp),0.0_dp,dp)
        work(1,p) = qp
    end subroutine apply_rotation

    pure logical function finite_complex(z) result(ok)
        complex(dp), intent(in) :: z !! complex scalar checked for finiteness
        ok = ieee_is_finite(real(z,dp)) .and. ieee_is_finite(aimag(z))
    end function finite_complex

end module scifort_linalg_banded_tridiagonal
