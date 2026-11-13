! SPDX-License-Identifier: MIT AND BSD-3-Clause
! Copyright (c) 1992-2023 The University of Tennessee and The University
! of Tennessee Research Foundation. All rights reserved.
! Copyright (c) 2000-2023 The University of California Berkeley. All rights reserved.
! Copyright (c) 2006-2023 The University of Colorado Denver. All rights reserved.
! Copyright (c) 2026 SciFort contributors
!
! QZ and supporting kernels adapted from Reference-LAPACK v3.12.1.
! Source attribution includes NAG Ltd.; full terms: LICENSES/lapack-BSD-3-Clause.txt.
! ZLARTG credits Weslley Pereira (2021), using Anderson (2017) safe scaling.
! This internal backend is restricted to normalized AAA arrowhead pencils.
module scifort_interpolate_aaa_spectrum
    use scifort_kinds, only: dp
    use scifort_linalg_hessenberg, only: hessenberg
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: aaa_arrowhead_eigenvalues, aaa_spectral_max_terms
    integer, parameter :: aaa_spectral_max_terms = 512
    real(dp), parameter :: safmin = real(radix(0.0_dp), dp)** &
                           max(minexponent(0.0_dp) - 1, 1 - maxexponent(0.0_dp))
    real(dp), parameter :: safmax = 1.0_dp/safmin
    real(dp), parameter :: zero = 0.0_dp, one = 1.0_dp, two = 2.0_dp
    complex(dp), parameter :: czero = cmplx(0.0_dp, 0.0_dp, dp)
contains
    subroutine aaa_arrowhead_eigenvalues(s, u, alpha, beta, status)
        implicit none
        complex(dp), intent(in) :: s(:) !! Distinct normalized support coordinates, modulus <=2.
        complex(dp), intent(in) :: u(:) !! Normalized secular weights, modulus <=2; length size(s).
        complex(dp), intent(out) :: alpha(:) !! Homogeneous numerators, length size(s)+1.
        complex(dp), intent(out) :: beta(:) !! Homogeneous denominators, length size(s)+1.
        integer, intent(out) :: status !! 0 success; 1 invalid/size, 2 allocation, 3 QZ failure.
        complex(dp), allocatable :: a(:, :), h(:, :), b(:, :), work(:)
        real(dp), allocatable :: rwork(:)
        complex(dp) :: q(1, 1), z(1, 1)
        integer :: m, n, i, info
        alpha = cmplx(0.0_dp, 0.0_dp, dp)
        beta = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        m = size(s)
        n = m + 1
        if (m < 1 .or. m > aaa_spectral_max_terms) return
        if (size(u) /= m .or. size(alpha) /= n .or. size(beta) /= n) return
        if (.not. all(ieee_is_finite(real(s, dp))) .or. &
            .not. all(ieee_is_finite(aimag(s)))) return
        if (.not. all(ieee_is_finite(real(u, dp))) .or. &
            .not. all(ieee_is_finite(aimag(u)))) return
        if (maxval(abs(s)) > 2.0_dp .or. maxval(abs(u)) > 2.0_dp) return
        allocate (a(n, n), h(n, n), b(n, n), work(n), rwork(n), stat=info)
        status = 2
        if (info /= 0) return
        a = cmplx(0.0_dp, 0.0_dp, dp)
        b = a
        a(1, 2:) = u
        a(2:, 1) = cmplx(1.0_dp, 0.0_dp, dp)
        do i = 1, m
            a(i + 1, i + 1) = s(i)
            b(i + 1, i + 1) = cmplx(1.0_dp, 0.0_dp, dp)
        end do
! Hessenberg similarities fix coordinate one. They commute with
! B=diag(0,I), so the second pencil factor remains exactly unchanged.
        call hessenberg(a, h, info)
        status = 3
        if (info /= 0) return
      call zhgeqz('e', 'n', 'n', n, 1, n, h, n, b, n, alpha, beta, q, 1, z, 1, work, n, rwork, info)
        if (info /= 0) return
        if (.not. all(ieee_is_finite(real(alpha, dp))) .or. &
            .not. all(ieee_is_finite(aimag(alpha)))) return
        if (.not. all(ieee_is_finite(real(beta, dp))) .or. &
            .not. all(ieee_is_finite(aimag(beta)))) return
        if (any(abs(alpha) + abs(beta) == 0.0_dp)) return
        status = 0
    end subroutine aaa_arrowhead_eigenvalues

    pure function abs1(x) result(value)
        implicit none
        complex(dp), intent(in) :: x !! Complex scalar for the componentwise one-norm.
        real(dp) :: value
        value = abs(real(x, dp)) + abs(aimag(x))
    end function abs1

    pure function abssq(x) result(value)
        implicit none
        complex(dp), intent(in) :: x !! Safely scaled complex scalar to square.
        real(dp) :: value
        value = real(x, dp)**2 + aimag(x)**2
    end function abssq

    pure function lsame(a, b) result(same)
        implicit none
        character(len=*), intent(in) :: a !! Internal lowercase option.
        character(len=*), intent(in) :: b !! Internal lowercase option to compare.
        logical :: same
        same = a == b
    end function lsame

    pure function dlamch(option) result(value)
        implicit none
        character(len=*), intent(in) :: option !! Internal lowercase machine-constant selection.
        real(dp) :: value
        select case (option(1:1))
        case ('e')
            value = epsilon(1.0_dp)*0.5_dp
        case ('s')
            value = tiny(1.0_dp)
        case ('o')
            value = huge(1.0_dp)
        case ('b')
            value = real(radix(1.0_dp), dp)
        case default
            value = 0.0_dp
        end select
    end function dlamch

    subroutine zlaset(option, m, n, offdiagonal, diagonal, a, lda)
        implicit none
        character(len=*), intent(in) :: option !! Internal full-matrix initialization option.
        integer, intent(in) :: m !! Number of rows to initialize.
        integer, intent(in) :: n !! Number of columns to initialize.
        complex(dp), intent(in) :: offdiagonal !! Value outside the diagonal.
        complex(dp), intent(in) :: diagonal !! Value on the diagonal.
        integer, intent(in) :: lda !! Leading dimension, at least m.
        complex(dp), intent(out) :: &
            a(lda, *) !! Matrix initialized for the unused Schur-vector path.
        integer :: i
        a(:m, :n) = offdiagonal
        do i = 1, min(m, n)
            a(i, i) = diagonal
        end do
    end subroutine zlaset

    subroutine zhgeqz(job, compq, compz, n, ilo, ihi, h, ldh, t, ldt, alpha, beta, q, &
                      ldq, z, ldz, work, lwork, rwork, info)
        implicit none
!
! -- lapack computational routine --
! -- lapack is a software package provided by univ. of tennessee,    --
! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
!
! .. scalar arguments ..
        character(len=1), intent(in) :: compq !! Internal left-vector selection; wrapper uses n.
        character(len=1), intent(in) :: compz !! Internal right-vector selection; wrapper uses n.
        character(len=1), intent(in) :: job !! Internal eigenvalue/Schur selection; wrapper uses e.
        integer, intent(in) :: ihi !! Last active Hessenberg row.
        integer, intent(in) :: ilo !! First active Hessenberg row, one in this wrapper.
        integer, intent(out) :: info !! Zero on convergence; positive if QZ iteration fails.
        integer, intent(in) :: ldh !! Leading dimension of the Hessenberg factor.
        integer, intent(in) :: ldq !! Leading dimension of the unused left basis.
        integer, intent(in) :: ldt !! Leading dimension of the triangular factor.
        integer, intent(in) :: ldz !! Leading dimension of the unused right basis.
        integer, intent(in) :: lwork !! Available complex workspace length.
        integer, intent(in) :: n !! Pencil order, validated by the wrapper.
! ..
! .. array arguments ..
        real(dp), intent(inout) :: &
            rwork(*) !! Unused real workspace retained from the upstream interface.
        complex(dp), intent(out) :: alpha(*) !! Homogeneous eigenvalue numerators.
        complex(dp), intent(out) :: beta(*) !! Homogeneous eigenvalue denominators.
        complex(dp), intent(inout) :: h(ldh, *) !! Hessenberg factor, overwritten by QZ iteration.
        complex(dp), intent(inout) :: &
            q(ldq, *) !! Unused left Schur basis for the eigenvalue-only path.
        complex(dp), intent(inout) :: t(ldt, *) !! Triangular factor, overwritten by QZ iteration.
        complex(dp), intent(inout) :: work(*) !! Complex QZ workspace, at least n entries.
        complex(dp), intent(inout) :: &
            z(ldz, *) !! Unused right Schur basis for the eigenvalue-only path.
! ..
!
! =====================================================================
!
! .. parameters ..
        complex(dp) :: czero
        complex(dp) :: cone
        parameter(czero=(0.0e+0_dp, 0.0e+0_dp), cone=(1.0e+0_dp, 0.0e+0_dp))
        real(dp) :: zero
        real(dp) :: one
        parameter(zero=0.0e+0_dp, one=1.0e+0_dp)
        real(dp) :: half
        parameter(half=0.5e+0_dp)
! ..
! .. local scalars ..
        logical :: ilazr2
        logical :: ilazro
        logical :: ilq
        logical :: ilschr
        logical :: ilz
        logical :: lquery
        integer :: icompq
        integer :: icompz
        integer :: ifirst
        integer :: ifrstm
        integer :: iiter
        integer :: ilast
        integer :: ilastm
        integer :: in
        integer :: ischur
        integer :: istart
        integer :: j
        integer :: jc
        integer :: jch
        integer :: jiter
        integer :: jr
        integer :: maxit
        real(dp) :: absb
        real(dp) :: anorm
        real(dp) :: ascale
        real(dp) :: atol
        real(dp) :: bnorm
        real(dp) :: bscale
        real(dp) :: btol
        real(dp) :: c
        real(dp) :: safmin
        real(dp) :: temp
        real(dp) :: temp2
        real(dp) :: tempr
        real(dp) :: ulp
        complex(dp) :: abi22
        complex(dp) :: ad11
        complex(dp) :: ad12
        complex(dp) :: ad21
        complex(dp) :: ad22
        complex(dp) :: ctemp
        complex(dp) :: ctemp2
        complex(dp) :: ctemp3
        complex(dp) :: eshift
        complex(dp) :: s
        complex(dp) :: shift
        complex(dp) :: signbc
        complex(dp) :: u12
        complex(dp) :: x
        complex(dp) :: abi12
        complex(dp) :: y
! ..
! .. external functions ..
! ..
! .. external subroutines ..
! ..
! .. intrinsic functions ..
! ..
! .. statement functions ..
! ..
! .. statement function definitions ..
! ..
! .. executable statements ..
!
! decode job, compq, compz
!
        if (lsame(job, 'e')) then
            ilschr = .false.
            ischur = 1
        else if (lsame(job, 's')) then
            ilschr = .true.
            ischur = 2
        else
            ilschr = .true.
            ischur = 0
        end if
!
        if (lsame(compq, 'n')) then
            ilq = .false.
            icompq = 1
        else if (lsame(compq, 'v')) then
            ilq = .true.
            icompq = 2
        else if (lsame(compq, 'i')) then
            ilq = .true.
            icompq = 3
        else
            ilq = .true.
            icompq = 0
        end if
!
        if (lsame(compz, 'n')) then
            ilz = .false.
            icompz = 1
        else if (lsame(compz, 'v')) then
            ilz = .true.
            icompz = 2
        else if (lsame(compz, 'i')) then
            ilz = .true.
            icompz = 3
        else
            ilz = .true.
            icompz = 0
        end if
!
! check argument values
!
        info = 0
        work(1) = max(1, n)
        lquery = (lwork .eq. -1)
        if (ischur .eq. 0) then
            info = -1
        else if (icompq .eq. 0) then
            info = -2
        else if (icompz .eq. 0) then
            info = -3
        else if (n .lt. 0) then
            info = -4
        else if (ilo .lt. 1) then
            info = -5
        else if (ihi .gt. n .or. ihi .lt. ilo - 1) then
            info = -6
        else if (ldh .lt. n) then
            info = -8
        else if (ldt .lt. n) then
            info = -10
        else if (ldq .lt. 1 .or. (ilq .and. ldq .lt. n)) then
            info = -14
        else if (ldz .lt. 1 .or. (ilz .and. ldz .lt. n)) then
            info = -16
        else if (lwork .lt. max(1, n) .and. .not. lquery) then
            info = -18
        end if
        if (info .ne. 0) then
            return
        else if (lquery) then
            return
        end if
!
! quick return if possible
!
! work( 1 ) = cmplx( 1 )
        if (n .le. 0) then
            work(1) = cmplx(1, kind=dp)
            return
        end if
!
! initialize q and z
!
        if (icompq .eq. 3) call zlaset('full', n, n, czero, cone, q, ldq)
        if (icompz .eq. 3) call zlaset('full', n, n, czero, cone, z, ldz)
!
! machine constants
!
        in = ihi + 1 - ilo
        safmin = dlamch('s')
        ulp = dlamch('e')*dlamch('b')
        anorm = sqrt(sum(abs(h(ilo:ihi, ilo:ihi))**2))
        bnorm = sqrt(sum(abs(t(ilo:ihi, ilo:ihi))**2))
        atol = max(safmin, ulp*anorm)
        btol = max(safmin, ulp*bnorm)
        ascale = one/max(safmin, anorm)
        bscale = one/max(safmin, bnorm)
!
!
! set eigenvalues ihi+1:n
!
        do j = ihi + 1, n
            absb = abs(t(j, j))
            if (absb .gt. safmin) then
                signbc = conjg(t(j, j)/absb)
                t(j, j) = absb
                if (ilschr) then
                    call zscal(j - 1, signbc, t(1, j), 1)
                    call zscal(j, signbc, h(1, j), 1)
                else
                    call zscal(1, signbc, h(j, j), 1)
                end if
                if (ilz) call zscal(n, signbc, z(1, j), 1)
            else
                t(j, j) = czero
            end if
            alpha(j) = h(j, j)
            beta(j) = t(j, j)
10      end do
!
! if ihi < ilo, skip qz steps
!
        if (ihi .lt. ilo) go to 190
!
! main qz iteration loop
!
! initialize dynamic indices
!
! eigenvalues ilast+1:n have been found.
! column operations modify rows ifrstm:whatever
! row operations modify columns whatever:ilastm
!
! if only eigenvalues are being computed, then
! ifrstm is the row of the last splitting row above row ilast;
! this is always at least ilo.
! iiter counts iterations since the last eigenvalue was found,
! to tell when to use an extraordinary shift.
! maxit is the maximum number of qz sweeps allowed.
!
        ilast = ihi
        if (ilschr) then
            ifrstm = 1
            ilastm = n
        else
            ifrstm = ilo
            ilastm = ihi
        end if
        iiter = 0
        eshift = czero
        maxit = 30*(ihi - ilo + 1)
!
        do jiter = 1, maxit
!
! check for too many iterations.
!
            if (jiter .gt. maxit) go to 180
!
! split the matrix if possible.
!
! two tests:
! 1: h(j,j-1)=0  or  j=ilo
! 2: t(j,j)=0
!
! special case: j=ilast
!
            if (ilast .eq. ilo) then
                go to 60
            else
                if (abs1(h(ilast, ilast - 1)) .le. max(safmin, ulp*(abs1(h(ilast, &
                                                     ilast)) + abs1(h(ilast - 1, ilast - 1))))) then
                    h(ilast, ilast - 1) = czero
                    go to 60
                end if
            end if
!
            if (abs(t(ilast, ilast)) .le. btol) then
                t(ilast, ilast) = czero
                go to 50
            end if
!
! general case: j<ilast
!
            do j = ilast - 1, ilo, -1
!
! test 1: for h(j,j-1)=0 or j=ilo
!
                if (j .eq. ilo) then
                    ilazro = .true.
                else
                    if (abs1(h(j, j - 1)) .le. max(safmin, ulp*(abs1(h(j, j)) + abs1(h(j - 1, &
                                                                                     j - 1))))) then
                        h(j, j - 1) = czero
                        ilazro = .true.
                    else
                        ilazro = .false.
                    end if
                end if
!
! test 2: for t(j,j)=0
!
                if (abs(t(j, j)) .lt. btol) then
                    t(j, j) = czero
!
! test 1a: check for 2 consecutive small subdiagonals in a
!
                    ilazr2 = .false.
                    if (.not. ilazro) then
                        if (abs1(h(j, j - 1))*(ascale*abs1(h(j + 1, j))) .le. abs1(h(j, &
                                                                  j))*(ascale*atol)) ilazr2 = .true.
                    end if
!
! if both tests pass (1 & 2), i.e., the leading diagonal
! element of b in the block is zero, split a 1x1 block off
! at the top. (i.e., at the j-th row/column) the leading
! diagonal element of the remainder can also be zero, so
! this may have to be done repeatedly.
!
                    if (ilazro .or. ilazr2) then
                    do jch = j, ilast - 1
                        ctemp = h(jch, jch)
                        call zlartg(ctemp, h(jch + 1, jch), c, s, h(jch, jch))
                        h(jch + 1, jch) = czero
                       call zrot(ilastm - jch, h(jch, jch + 1), ldh, h(jch + 1, jch + 1), ldh, c, s)
                       call zrot(ilastm - jch, t(jch, jch + 1), ldt, t(jch + 1, jch + 1), ldt, c, s)
                        if (ilq) call zrot(n, q(1, jch), 1, q(1, jch + 1), 1, c, conjg(s))
                        if (ilazr2) h(jch, jch - 1) = h(jch, jch - 1)*c
                        ilazr2 = .false.
                        if (abs1(t(jch + 1, jch + 1)) .ge. btol) then
                        if (jch + 1 .ge. ilast) then
                            go to 60
                        else
                            ifirst = jch + 1
                            go to 70
                        end if
                        end if
                        t(jch + 1, jch + 1) = czero
20                  end do
                    go to 50
                    else
!
! only test 2 passed -- chase the zero to t(ilast,ilast)
! then process as in the case t(ilast,ilast)=0
!
                    do jch = j, ilast - 1
                        ctemp = t(jch, jch + 1)
                        call zlartg(ctemp, t(jch + 1, jch + 1), c, s, t(jch, jch + 1))
                        t(jch + 1, jch + 1) = czero
             if (jch .lt. ilastm - 1) call zrot(ilastm - jch - 1, t(jch, jch + 2), ldt, t(jch + 1, &
                                                                                jch + 2), ldt, c, s)
                   call zrot(ilastm - jch + 2, h(jch, jch - 1), ldh, h(jch + 1, jch - 1), ldh, c, s)
                        if (ilq) call zrot(n, q(1, jch), 1, q(1, jch + 1), 1, c, conjg(s))
                        ctemp = h(jch + 1, jch)
                        call zlartg(ctemp, h(jch + 1, jch - 1), c, s, h(jch + 1, jch))
                        h(jch + 1, jch - 1) = czero
                        call zrot(jch + 1 - ifrstm, h(ifrstm, jch), 1, h(ifrstm, jch - 1), 1, c, s)
                        call zrot(jch - ifrstm, t(ifrstm, jch), 1, t(ifrstm, jch - 1), 1, c, s)
                        if (ilz) call zrot(n, z(1, jch), 1, z(1, jch - 1), 1, c, s)
30                  end do
                    go to 50
                    end if
                else if (ilazro) then
!
! only test 1 passed -- work on j:ilast
!
                    ifirst = j
                    go to 70
                end if
!
! neither test passed -- try next j
!
40          end do
!
! (drop-through is "impossible")
!
            info = 2*n + 1
            go to 210
!
! t(ilast,ilast)=0 -- clear h(ilast,ilast-1) to split off a
! 1x1 block.
!
50          continue
            ctemp = h(ilast, ilast)
            call zlartg(ctemp, h(ilast, ilast - 1), c, s, h(ilast, ilast))
            h(ilast, ilast - 1) = czero
            call zrot(ilast - ifrstm, h(ifrstm, ilast), 1, h(ifrstm, ilast - 1), 1, c, s)
            call zrot(ilast - ifrstm, t(ifrstm, ilast), 1, t(ifrstm, ilast - 1), 1, c, s)
            if (ilz) call zrot(n, z(1, ilast), 1, z(1, ilast - 1), 1, c, s)
!
! h(ilast,ilast-1)=0 -- standardize b, set alpha and beta
!
60          continue
            absb = abs(t(ilast, ilast))
            if (absb .gt. safmin) then
                signbc = conjg(t(ilast, ilast)/absb)
                t(ilast, ilast) = absb
                if (ilschr) then
                    call zscal(ilast - ifrstm, signbc, t(ifrstm, ilast), 1)
                    call zscal(ilast + 1 - ifrstm, signbc, h(ifrstm, ilast), 1)
                else
                    call zscal(1, signbc, h(ilast, ilast), 1)
                end if
                if (ilz) call zscal(n, signbc, z(1, ilast), 1)
            else
                t(ilast, ilast) = czero
            end if
            alpha(ilast) = h(ilast, ilast)
            beta(ilast) = t(ilast, ilast)
!
! go to next block -- exit if finished.
!
            ilast = ilast - 1
            if (ilast .lt. ilo) go to 190
!
! reset counters
!
            iiter = 0
            eshift = czero
            if (.not. ilschr) then
                ilastm = ilast
                if (ifrstm .gt. ilast) ifrstm = ilo
            end if
            go to 160
!
! qz step
!
! this iteration only involves rows/columns ifirst:ilast.  we
! assume ifirst < ilast, and that the diagonal of b is non-zero.
!
70          continue
            iiter = iiter + 1
            if (.not. ilschr) then
                ifrstm = ifirst
            end if
!
! compute the shift.
!
! at this point, ifirst < ilast, and the diagonal elements of
! t(ifirst:ilast,ifirst,ilast) are larger than btol (in
! magnitude)
!
            if ((iiter/10)*10 .ne. iiter) then
!
! the wilkinson shift (aep p.512), i.e., the eigenvalue of
! the bottom-right 2x2 block of a inv(b) which is nearest to
! the bottom-right element.
!
! we factor b as u*d, where u has unit diagonals, and
! compute (a*inv(d))*inv(u).
!
                u12 = (bscale*t(ilast - 1, ilast))/(bscale*t(ilast, ilast))
                ad11 = (ascale*h(ilast - 1, ilast - 1))/(bscale*t(ilast - 1, ilast - 1))
                ad21 = (ascale*h(ilast, ilast - 1))/(bscale*t(ilast - 1, ilast - 1))
                ad12 = (ascale*h(ilast - 1, ilast))/(bscale*t(ilast, ilast))
                ad22 = (ascale*h(ilast, ilast))/(bscale*t(ilast, ilast))
                abi22 = ad22 - u12*ad21
                abi12 = ad12 - u12*ad11
!
                shift = abi22
                ctemp = sqrt(abi12)*sqrt(ad21)
                temp = abs1(ctemp)
                if (ctemp .ne. zero) then
                    x = half*(ad11 - shift)
                    temp2 = abs1(x)
                    temp = max(temp, abs1(x))
                    y = temp*sqrt((x/temp)**2 + (ctemp/temp)**2)
                    if (temp2 .gt. zero) then
                       if (real(x/temp2, dp)*real(y, dp) + aimag(x/temp2)*aimag(y) .lt. zero) y = -y
                    end if
                    shift = shift - ctemp*zladiv(ctemp, (x + y))
                end if
            else
!
! exceptional shift.  chosen for no particularly good reason.
!
                if ((iiter/20)*20 .eq. iiter .and. bscale*abs1(t(ilast, ilast)) .gt. safmin) then
                    eshift = eshift + (ascale*h(ilast, ilast))/(bscale*t(ilast, ilast))
                else
                    eshift = eshift + (ascale*h(ilast, ilast - 1))/(bscale*t(ilast - 1, ilast - 1))
                end if
                shift = eshift
            end if
!
! now check for two consecutive small subdiagonals.
!
            do j = ilast - 1, ifirst + 1, -1
                istart = j
                ctemp = ascale*h(j, j) - shift*(bscale*t(j, j))
                temp = abs1(ctemp)
                temp2 = ascale*abs1(h(j + 1, j))
                tempr = max(temp, temp2)
                if (tempr .lt. one .and. tempr .ne. zero) then
                    temp = temp/tempr
                    temp2 = temp2/tempr
                end if
                if (abs1(h(j, j - 1))*temp2 .le. temp*atol) go to 90
80          end do
!
            istart = ifirst
            ctemp = ascale*h(ifirst, ifirst) - shift*(bscale*t(ifirst, ifirst))
90          continue
!
! do an implicit-shift qz sweep.
!
! initial q
!
            ctemp2 = ascale*h(istart + 1, istart)
            call zlartg(ctemp, ctemp2, c, s, ctemp3)
!
! sweep
!
            do j = istart, ilast - 1
            if (j .gt. istart) then
                ctemp = h(j, j - 1)
                call zlartg(ctemp, h(j + 1, j - 1), c, s, h(j, j - 1))
                h(j + 1, j - 1) = czero
            end if
!
            do jc = j, ilastm
                ctemp = c*h(j, jc) + s*h(j + 1, jc)
                h(j + 1, jc) = -conjg(s)*h(j, jc) + c*h(j + 1, jc)
                h(j, jc) = ctemp
                ctemp2 = c*t(j, jc) + s*t(j + 1, jc)
                t(j + 1, jc) = -conjg(s)*t(j, jc) + c*t(j + 1, jc)
                t(j, jc) = ctemp2
100         end do
            if (ilq) then
            do jr = 1, n
                ctemp = c*q(jr, j) + conjg(s)*q(jr, j + 1)
                q(jr, j + 1) = -s*q(jr, j) + c*q(jr, j + 1)
                q(jr, j) = ctemp
110         end do
            end if
!
            ctemp = t(j + 1, j + 1)
            call zlartg(ctemp, t(j + 1, j), c, s, t(j + 1, j + 1))
            t(j + 1, j) = czero
!
            do jr = ifrstm, min(j + 2, ilast)
                ctemp = c*h(jr, j + 1) + s*h(jr, j)
                h(jr, j) = -conjg(s)*h(jr, j + 1) + c*h(jr, j)
                h(jr, j + 1) = ctemp
120         end do
            do jr = ifrstm, j
                ctemp = c*t(jr, j + 1) + s*t(jr, j)
                t(jr, j) = -conjg(s)*t(jr, j + 1) + c*t(jr, j)
                t(jr, j + 1) = ctemp
130         end do
            if (ilz) then
            do jr = 1, n
                ctemp = c*z(jr, j + 1) + s*z(jr, j)
                z(jr, j) = -conjg(s)*z(jr, j + 1) + c*z(jr, j)
                z(jr, j + 1) = ctemp
140         end do
            end if
150         end do
!
160         continue
!
170     end do
!
! drop-through = non-convergence
!
180     continue
        info = ilast
        go to 210
!
! successful completion of all qz steps
!
190     continue
!
! set eigenvalues 1:ilo-1
!
        do j = 1, ilo - 1
            absb = abs(t(j, j))
            if (absb .gt. safmin) then
                signbc = conjg(t(j, j)/absb)
                t(j, j) = absb
                if (ilschr) then
                    call zscal(j - 1, signbc, t(1, j), 1)
                    call zscal(j, signbc, h(1, j), 1)
                else
                    call zscal(1, signbc, h(j, j), 1)
                end if
                if (ilz) call zscal(n, signbc, z(1, j), 1)
            else
                t(j, j) = czero
            end if
            alpha(j) = h(j, j)
            beta(j) = t(j, j)
200     end do
!
! normal termination
!
        info = 0
!
! exit (other than argument error) -- return optimal workspace size
!
210     continue
        work(1) = cmplx(n, kind=dp)
        return
!
! end of zhgeqz
!
    end subroutine zhgeqz
    subroutine zrot(n, cx, incx, cy, incy, c, s)
        implicit none
!
! -- lapack auxiliary routine --
! -- lapack is a software package provided by univ. of tennessee,    --
! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
!
! .. scalar arguments ..
        integer, intent(in) :: incx !! Nonzero stride through the first vector.
        integer, intent(in) :: incy !! Nonzero stride through the second vector.
        integer, intent(in) :: n !! Number of rotated complex pairs.
        real(dp), intent(in) :: c !! Real cosine of the unitary rotation.
        complex(dp), intent(in) :: s !! Complex sine of the unitary rotation.
! ..
! .. array arguments ..
        complex(dp), intent(inout) :: cx(*) !! First strided vector, overwritten by its rotation.
        complex(dp), intent(inout) :: cy(*) !! Second strided vector, overwritten by its rotation.
! ..
!
! =====================================================================
!
! .. local scalars ..
        integer :: i
        integer :: ix
        integer :: iy
        complex(dp) :: stemp
! ..
! .. intrinsic functions ..
! ..
! .. executable statements ..
!
        if (n .le. 0) return
        if (incx .eq. 1 .and. incy .eq. 1) go to 20
!
! code for unequal increments or equal increments not equal to 1
!
        ix = 1
        iy = 1
        if (incx .lt. 0) ix = (-n + 1)*incx + 1
        if (incy .lt. 0) iy = (-n + 1)*incy + 1
        do i = 1, n
            stemp = c*cx(ix) + s*cy(iy)
            cy(iy) = c*cy(iy) - conjg(s)*cx(ix)
            cx(ix) = stemp
            ix = ix + incx
            iy = iy + incy
10      end do
        return
!
! code for both increments equal to 1
!
20      continue
        do i = 1, n
            stemp = c*cx(i) + s*cy(i)
            cy(i) = c*cy(i) - conjg(s)*cx(i)
            cx(i) = stemp
30      end do
        return
    end subroutine zrot
    subroutine zscal(n, za, zx, incx)
        implicit none
!
! -- reference blas level1 routine --
! -- reference blas is a software package provided by univ. of tennessee,    --
! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
!
! .. scalar arguments ..
        complex(dp), intent(in) :: za !! Complex scaling factor.
        integer, intent(in) :: incx !! Positive stride through the complex vector.
        integer, intent(in) :: n !! Number of scaled vector entries.
! ..
! .. array arguments ..
        complex(dp), intent(inout) :: zx(*) !! Strided complex vector, scaled in place.
! ..
!
! =====================================================================
!
! .. local scalars ..
        integer :: i
        integer :: nincx
! ..
! .. parameters ..
        complex(dp) :: one
        parameter(one=(1.0e+0_dp, 0.0e+0_dp))
! ..
        if (n .le. 0 .or. incx .le. 0 .or. za .eq. one) return
        if (incx .eq. 1) then
!
! code for increment equal to 1
!
            do i = 1, n
                zx(i) = za*zx(i)
            end do
        else
!
! code for increment not equal to 1
!
            nincx = n*incx
            do i = 1, nincx, incx
                zx(i) = za*zx(i)
            end do
        end if
        return
!
! end of zscal
!
    end subroutine zscal
    complex(dp) function zladiv(x, y) result(value)
        implicit none
!
! -- lapack auxiliary routine --
! -- lapack is a software package provided by univ. of tennessee,    --
! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
!
! .. scalar arguments ..
        complex(dp), intent(in) :: x !! Complex numerator.
        complex(dp), intent(in) :: y !! Nonzero complex denominator.
! ..
!
! =====================================================================
!
! .. local scalars ..
        real(dp) :: zi
        real(dp) :: zr
! ..
! .. external subroutines ..
! ..
! .. intrinsic functions ..
! ..
! .. executable statements ..
!
        call dladiv(real(x, dp), aimag(x), real(y, dp), aimag(y), zr, zi)
        value = cmplx(zr, zi, kind=dp)
!
        return
!
! end of zladiv
!
    end function zladiv
    subroutine dladiv(a, b, c, d, p, q)
        implicit none
!
! -- lapack auxiliary routine --
! -- lapack is a software package provided by univ. of tennessee,    --
! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
!
! .. scalar arguments ..
        real(dp), intent(in) :: a !! Real part of the numerator.
        real(dp), intent(in) :: b !! Imaginary part of the numerator.
        real(dp), intent(in) :: c !! Real part of the nonzero denominator.
        real(dp), intent(in) :: d !! Imaginary part of the nonzero denominator.
        real(dp), intent(out) :: p !! Real part of the scaled quotient.
        real(dp), intent(out) :: q !! Imaginary part of the scaled quotient.
! ..
!
! =====================================================================
!
! .. parameters ..
        real(dp) :: bs
        parameter(bs=2.0e0_dp)
        real(dp) :: half
        parameter(half=0.5e0_dp)
        real(dp) :: two
        parameter(two=2.0e0_dp)
!
! .. local scalars ..
        real(dp) :: aa
        real(dp) :: bb
        real(dp) :: cc
        real(dp) :: dd
        real(dp) :: ab
        real(dp) :: cd
        real(dp) :: s
        real(dp) :: ov
        real(dp) :: un
        real(dp) :: be
        real(dp) :: eps
! ..
! .. external functions ..
! ..
! .. external subroutines ..
! ..
! .. intrinsic functions ..
! ..
! .. executable statements ..
!
        aa = a
        bb = b
        cc = c
        dd = d
        ab = max(abs(a), abs(b))
        cd = max(abs(c), abs(d))
        s = 1.0e0_dp
        ov = dlamch('overflow threshold')
        un = dlamch('safe minimum')
        eps = dlamch('epsilon')
        be = bs/(eps*eps)
        if (ab >= half*ov) then
            aa = half*aa
            bb = half*bb
            s = two*s
        end if
        if (cd >= half*ov) then
            cc = half*cc
            dd = half*dd
            s = half*s
        end if
        if (ab <= un*bs/eps) then
            aa = aa*be
            bb = bb*be
            s = s/be
        end if
        if (cd <= un*bs/eps) then
            cc = cc*be
            dd = dd*be
            s = s*be
        end if
        if (abs(d) .le. abs(c)) then
            call dladiv1(aa, bb, cc, dd, p, q)
        else
            call dladiv1(bb, aa, dd, cc, p, q)
            q = -q
        end if
        p = p*s
        q = q*s
!
        return
!
! end of dladiv
!
    end subroutine dladiv
! > \ingroup ladiv
    subroutine dladiv1(a, b, c, d, p, q)
        implicit none
!
! -- lapack auxiliary routine --
! -- lapack is a software package provided by univ. of tennessee,    --
! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
!
! .. scalar arguments ..
        real(dp), intent(inout) :: &
            a !! Numerator component, temporarily negated for the second quotient.
        real(dp), intent(in) :: b !! Second numerator component.
        real(dp), intent(in) :: c !! Dominant nonzero denominator component.
        real(dp), intent(in) :: d !! Other denominator component.
        real(dp), intent(out) :: p !! First quotient component.
        real(dp), intent(out) :: q !! Second quotient component.
! ..
!
! =====================================================================
!
! .. parameters ..
        real(dp) :: one
        parameter(one=1.0e0_dp)
!
! .. local scalars ..
        real(dp) :: r
        real(dp) :: t
! ..
! .. external functions ..
! ..
! .. executable statements ..
!
        r = d/c
        t = one/(c + d*r)
        p = dladiv2(a, b, c, d, r, t)
        a = -a
        q = dladiv2(b, a, c, d, r, t)
!
        return
!
! end of dladiv1
!
    end subroutine dladiv1
! > \ingroup ladiv
    real(dp) function dladiv2(a, b, c, d, r, t) result(value)
        implicit none
!
! -- lapack auxiliary routine --
! -- lapack is a software package provided by univ. of tennessee,    --
! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
!
! .. scalar arguments ..
        real(dp), intent(in) :: a !! First numerator component.
        real(dp), intent(in) :: b !! Second numerator component.
        real(dp), intent(in) :: c !! Dominant nonzero denominator component.
        real(dp), intent(in) :: d !! Other denominator component.
        real(dp), intent(in) :: r !! Ratio of denominator components.
        real(dp), intent(in) :: t !! Reciprocal scaled denominator.
! ..
!
! =====================================================================
!
! .. parameters ..
        real(dp) :: zero
        parameter(zero=0.0e0_dp)
!
! .. local scalars ..
        real(dp) :: br
! ..
! .. executable statements ..
!
        if (r .ne. zero) then
            br = b*r
            if (br .ne. zero) then
                value = (a + br)*t
            else
                value = a*t + (b*t)*r
            end if
        else
            value = (a + d*(b/c))*t
        end if
!
        return
!
! end of dladiv2
!
    end function dladiv2
    subroutine zlartg(f, g, c, s, r)
        implicit none
!
!  -- lapack auxiliary routine --
!  -- lapack is a software package provided by univ. of tennessee,    --
!  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
!     february 2021
!
!  .. scalar arguments ..
        real(dp), intent(out) :: c !! Real rotation cosine.
        complex(dp), intent(in) :: f !! First complex component to rotate.
        complex(dp), intent(in) :: g !! Second complex component to eliminate.
        complex(dp), intent(out) :: r !! Rotated first component.
        complex(dp), intent(out) :: s !! Complex rotation sine.
!  ..
!  .. local scalars ..
        real(dp) :: d
        real(dp) :: f1
        real(dp) :: f2
        real(dp) :: g1
        real(dp) :: g2
        real(dp) :: h2
        real(dp) :: u
        real(dp) :: v
        real(dp) :: w
        real(dp) :: rtmin
        real(dp) :: rtmax
        complex(dp) :: fs
        complex(dp) :: gs
        complex(dp) :: t
!  ..
!  .. intrinsic functions ..
!  ..
!  .. statement functions ..
!  ..
!  .. statement function definitions ..
!  ..
!  .. constants ..
        rtmin = sqrt(safmin)
!  ..
!  .. executable statements ..
!
        if (g == czero) then
            c = one
            s = czero
            r = f
        else if (f == czero) then
            c = zero
            if (real(g) == zero) then
                r = abs(aimag(g))
                s = conjg(g)/r
            elseif (aimag(g) == zero) then
                r = abs(real(g))
                s = conjg(g)/r
            else
                g1 = max(abs(real(g)), abs(aimag(g)))
                rtmax = sqrt(safmax/2)
                if (g1 > rtmin .and. g1 < rtmax) then
!
!        use unscaled algorithm
!
!           the following two lines can be replaced by `d = abs( g )`.
!           this algorithm do not use the intrinsic complex abs.
                    g2 = abssq(g)
                    d = sqrt(g2)
                    s = conjg(g)/d
                    r = d
                else
!
!        use scaled algorithm
!
                    u = min(safmax, max(safmin, g1))
                    gs = g/u
!           the following two lines can be replaced by `d = abs( gs )`.
!           this algorithm do not use the intrinsic complex abs.
                    g2 = abssq(gs)
                    d = sqrt(g2)
                    s = conjg(gs)/d
                    r = d*u
                end if
            end if
        else
            f1 = max(abs(real(f)), abs(aimag(f)))
            g1 = max(abs(real(g)), abs(aimag(g)))
            rtmax = sqrt(safmax/4)
            if (f1 > rtmin .and. f1 < rtmax .and. &
                g1 > rtmin .and. g1 < rtmax) then
!
!        use unscaled algorithm
!
                f2 = abssq(f)
                g2 = abssq(g)
                h2 = f2 + g2
! safmin <= f2 <= h2 <= safmax
                if (f2 >= h2*safmin) then
! safmin <= f2/h2 <= 1, and h2/f2 is finite
                    c = sqrt(f2/h2)
                    r = f/c
                    rtmax = rtmax*2
                    if (f2 > rtmin .and. h2 < rtmax) then
! safmin <= sqrt( f2*h2 ) <= safmax
                        s = conjg(g)*(f/sqrt(f2*h2))
                    else
                        s = conjg(g)*(r/h2)
                    end if
                else
! f2/h2 <= safmin may be subnormal, and h2/f2 may overflow.
! moreover,
!  safmin <= f2*f2 * safmax < f2 * h2 < h2*h2 * safmin <= safmax,
!  sqrt(safmin) <= sqrt(f2 * h2) <= sqrt(safmax).
! also,
!  g2 >> f2, which means that h2 = g2.
                    d = sqrt(f2*h2)
                    c = f2/d
                    if (c >= safmin) then
                        r = f/c
                    else
! f2 / sqrt(f2 * h2) < safmin, then
!  sqrt(safmin) <= f2 * sqrt(safmax) <= h2 / sqrt(f2 * h2) <= h2 * (safmin / f2) <= h2 <= safmax
                        r = f*(h2/d)
                    end if
                    s = conjg(g)*(f/d)
                end if
            else
!
!        use scaled algorithm
!
                u = min(safmax, max(safmin, f1, g1))
                gs = g/u
                g2 = abssq(gs)
                if (f1/u < rtmin) then
!
!           f is not well-scaled when scaled by g1.
!           use a different scaling for f.
!
                    v = min(safmax, max(safmin, f1))
                    w = v/u
                    fs = f/v
                    f2 = abssq(fs)
                    h2 = f2*w**2 + g2
                else
!
!           otherwise use the same scaling for f and g.
!
                    w = one
                    fs = f/u
                    f2 = abssq(fs)
                    h2 = f2 + g2
                end if
! safmin <= f2 <= h2 <= safmax
                if (f2 >= h2*safmin) then
! safmin <= f2/h2 <= 1, and h2/f2 is finite
                    c = sqrt(f2/h2)
                    r = fs/c
                    rtmax = rtmax*2
                    if (f2 > rtmin .and. h2 < rtmax) then
! safmin <= sqrt( f2*h2 ) <= safmax
                        s = conjg(gs)*(fs/sqrt(f2*h2))
                    else
                        s = conjg(gs)*(r/h2)
                    end if
                else
! f2/h2 <= safmin may be subnormal, and h2/f2 may overflow.
! moreover,
!  safmin <= f2*f2 * safmax < f2 * h2 < h2*h2 * safmin <= safmax,
!  sqrt(safmin) <= sqrt(f2 * h2) <= sqrt(safmax).
! also,
!  g2 >> f2, which means that h2 = g2.
                    d = sqrt(f2*h2)
                    c = f2/d
                    if (c >= safmin) then
                        r = fs/c
                    else
! f2 / sqrt(f2 * h2) < safmin, then
!  sqrt(safmin) <= f2 * sqrt(safmax) <= h2 / sqrt(f2 * h2) <= h2 * (safmin / f2) <= h2 <= safmax
                        r = fs*(h2/d)
                    end if
                    s = conjg(gs)*(fs/d)
                end if
! rescale c and r
                c = c*w
                r = r*u
            end if
        end if
        return
    end subroutine zlartg
end module scifort_interpolate_aaa_spectrum
