! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Transfer-function/ZPK/SOS conversions. Coefficients in tf2zpk
! and zpk2tf use descending positive powers of z (or s). Default
! digital SOS rows represent ascending negative powers of z. With
! pairing=minimal the rows may have leading zero coefficients and are
! polynomial factors, not necessarily suitable for direct digital filtering.
module scifort_signal_zpk_conversion
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_signal_analog, only : polynomial_roots
    implicit none
    private

    integer, parameter, public :: zpk_conversion_success = 0
    integer, parameter, public :: zpk_conversion_invalid_input = 1
    integer, parameter, public :: zpk_conversion_numerical_failure = 2

    type, public :: filter_zpk_result
        complex(dp), allocatable :: zeros(:) !! finite roots of the numerator
        complex(dp), allocatable :: poles(:) !! finite roots of the denominator
        complex(dp) :: gain = (0.0_dp, 0.0_dp) !! ratio of leading polynomial coefficients
        integer :: status = zpk_conversion_invalid_input
    end type filter_zpk_result

    type, public :: filter_tf_result
        complex(dp), allocatable :: b(:) !! numerator, highest polynomial power first
        complex(dp), allocatable :: a(:) !! denominator, highest polynomial power first
        integer :: status = zpk_conversion_invalid_input
    end type filter_tf_result

    type, public :: filter_sos_result
        complex(dp), allocatable :: sos(:,:) !! SOS rows (nsections,6); a0 can be zero for minimal/analog
        integer :: status = zpk_conversion_invalid_input
    end type filter_sos_result

    interface tf2sos
        module procedure tf2sos_real, tf2sos_complex
    end interface tf2sos
    interface tf2zpk
        module procedure tf2zpk_real, tf2zpk_complex
    end interface tf2zpk
    interface zpk2tf
        module procedure zpk2tf_complex, zpk2tf_real_gain, zpk2tf_real_roots, &
            zpk2tf_real_roots_complex_gain
    end interface zpk2tf
    interface zpk2sos
        module procedure zpk2sos_complex, zpk2sos_real_gain, zpk2sos_real_roots, &
            zpk2sos_real_roots_complex_gain
    end interface zpk2sos
    interface sos2zpk
        module procedure sos2zpk_real, sos2zpk_complex
    end interface sos2zpk
    public :: tf2sos, tf2zpk, zpk2tf, zpk2sos, sos2zpk

contains

    pure logical function finite_complex(values) result(ok)
        complex(dp), intent(in) :: values(:) !! input coefficients or roots
        ok = all(ieee_is_finite(real(values,dp))) .and. all(ieee_is_finite(aimag(values)))
    end function finite_complex

    function tf2zpk_real(b, a) result(out)
        real(dp), intent(in) :: b(:) !! numerator polynomial, descending powers
        real(dp), intent(in) :: a(:) !! denominator polynomial, descending powers
        type(filter_zpk_result) :: out
        out = tf2zpk_complex(cmplx(b,0.0_dp,dp),cmplx(a,0.0_dp,dp))
    end function tf2zpk_real

    function tf2zpk_complex(b, a) result(out)
        complex(dp), intent(in) :: b(:) !! numerator polynomial, descending powers
        complex(dp), intent(in) :: a(:) !! denominator polynomial, descending powers
        type(filter_zpk_result) :: out
        complex(dp), allocatable :: num(:), den(:)
        integer :: first_a, first_b
        logical :: ok

        if (size(a) == 0 .or. size(b) == 0) return
        if (.not. finite_complex(a) .or. .not. finite_complex(b)) return
        first_a = 1
        do while (first_a <= size(a))
            if (abs(a(first_a)) /= 0.0_dp) exit
            first_a = first_a + 1
        end do
        if (first_a > size(a)) return
        den = a(first_a:)/a(first_a)
        if (.not. finite_complex(den)) then
            out%status = zpk_conversion_numerical_failure
            return
        end if
        call polynomial_roots(den,out%poles,ok)
        if (.not. ok) then
            out%status = zpk_conversion_numerical_failure
            return
        end if
        first_b = 1
        do while (first_b <= size(b))
            if (abs(b(first_b)) /= 0.0_dp) exit
            first_b = first_b + 1
        end do
        if (first_b > size(b)) then
            allocate(out%zeros(0))
            out%gain = (0.0_dp,0.0_dp)
            out%status = zpk_conversion_success
            return
        end if
        num = b(first_b:)/a(first_a)
        out%gain = num(1)
        if (.not. finite_complex(num)) then
            out%status = zpk_conversion_numerical_failure
            return
        end if
        num = num/num(1)
        if (.not. finite_complex(num)) then
            out%status = zpk_conversion_numerical_failure
            return
        end if
        call polynomial_roots(num,out%zeros,ok)
        if (.not. ok) then
            out%status = zpk_conversion_numerical_failure
            return
        end if
        out%status = zpk_conversion_success
    end function tf2zpk_complex

    function tf2sos_real(b, a, pairing, analog) result(out)
        real(dp), intent(in) :: b(:) !! numerator polynomial coefficients
        real(dp), intent(in) :: a(:) !! denominator polynomial coefficients
        character(*), intent(in), optional :: pairing !! nearest, keep_odd, or minimal
        logical, intent(in), optional :: analog !! true for analog s-plane sections
        type(filter_sos_result) :: out

        out = tf2sos_complex(cmplx(b,0.0_dp,dp),cmplx(a,0.0_dp,dp),pairing,analog)
    end function tf2sos_real

    function tf2sos_complex(b, a, pairing, analog) result(out)
        ! SciPy's tf2sos is tf2zpk followed by zpk2sos. Digital SOS
        ! coefficients are in increasing powers of z**(-1). Analog SOS
        ! coefficients are in descending powers of s. No polynomial
        ! expansion is done after finding the roots.
        complex(dp), intent(in) :: b(:) !! numerator polynomial coefficients
        complex(dp), intent(in) :: a(:) !! denominator polynomial coefficients
        character(*), intent(in), optional :: pairing !! nearest, keep_odd, or minimal
        logical, intent(in), optional :: analog !! analog mode requires minimal pairing
        type(filter_sos_result) :: out
        type(filter_zpk_result) :: zp
        character(:), allocatable :: mode
        logical :: analog_mode

        analog_mode = .false.
        if (present(analog)) analog_mode = analog
        mode = 'nearest'
        if (analog_mode) mode = 'minimal'
        if (present(pairing)) mode = trim(pairing)
        if (mode /= 'nearest' .and. mode /= 'keep_odd' .and. mode /= 'minimal') return
        if (analog_mode .and. mode /= 'minimal') return
        zp = tf2zpk_complex(b,a)
        if (zp%status /= zpk_conversion_success) then
            out%status = zp%status
            return
        end if
        out = zpk2sos_complex(zp%zeros,zp%poles,zp%gain,pairing=mode)
        if (out%status /= zpk_conversion_success) return
        if (analog_mode .and. size(zp%poles) == 0) then
            ! For constant analog transfer functions use a degree-zero
            ! numerator and denominator, not artificial s**2 factors.
            out%sos(1,:) = cmplx([0.0_dp,0.0_dp,1.0_dp, &
                0.0_dp,0.0_dp,1.0_dp],0.0_dp,dp)
            out%sos(1,3) = zp%gain
        end if
    end function tf2sos_complex

    pure function monic_from_roots(roots) result(poly)
        complex(dp), intent(in) :: roots(:) !! polynomial roots
        complex(dp), allocatable :: poly(:), next(:)
        integer :: j, n
        poly = [cmplx(1.0_dp,0.0_dp,dp)]
        do j = 1, size(roots)
            n = size(poly)
            allocate(next(n+1))
            next = (0.0_dp,0.0_dp)
            next(:n) = poly
            next(2:) = next(2:) - roots(j)*poly
            call move_alloc(next,poly)
        end do
    end function monic_from_roots

    function zpk2tf_real_roots(z, p, k) result(out)
        real(dp), intent(in) :: z(:) !! real zeros
        real(dp), intent(in) :: p(:) !! real poles
        real(dp), intent(in) :: k !! real gain
        type(filter_tf_result) :: out
        out = zpk2tf_complex(cmplx(z,0.0_dp,dp),cmplx(p,0.0_dp,dp),cmplx(k,0.0_dp,dp))
    end function zpk2tf_real_roots

    function zpk2tf_real_roots_complex_gain(z, p, k) result(out)
        real(dp), intent(in) :: z(:) !! real zeros
        real(dp), intent(in) :: p(:) !! real poles
        complex(dp), intent(in) :: k !! complex gain
        type(filter_tf_result) :: out
        out = zpk2tf_complex(cmplx(z,0.0_dp,dp),cmplx(p,0.0_dp,dp),k)
    end function zpk2tf_real_roots_complex_gain

    function zpk2tf_real_gain(z, p, k) result(out)
        complex(dp), intent(in) :: z(:) !! numerator roots
        complex(dp), intent(in) :: p(:) !! denominator roots
        real(dp), intent(in) :: k !! real gain
        type(filter_tf_result) :: out
        out = zpk2tf_complex(z,p,cmplx(k,0.0_dp,dp))
    end function zpk2tf_real_gain

    function zpk2tf_complex(z, p, k) result(out)
        complex(dp), intent(in) :: z(:) !! numerator roots
        complex(dp), intent(in) :: p(:) !! denominator roots
        complex(dp), intent(in) :: k !! complex gain
        type(filter_tf_result) :: out

        if (.not. finite_complex(z) .or. .not. finite_complex(p)) return
        if (.not. ieee_is_finite(real(k,dp)) .or. .not. ieee_is_finite(aimag(k))) return
        out%b = k*monic_from_roots(z)
        out%a = monic_from_roots(p)
        if (.not. finite_complex(out%b) .or. .not. finite_complex(out%a)) then
            out%status = zpk_conversion_numerical_failure
            return
        end if
        out%status = zpk_conversion_success
    end function zpk2tf_complex

    function zpk2sos_real_roots(z, p, k, pairing) result(out)
        real(dp), intent(in) :: z(:) !! real zeros
        real(dp), intent(in) :: p(:) !! real poles
        real(dp), intent(in) :: k !! real gain
        character(*), intent(in), optional :: pairing !! nearest (default), keep_odd, or minimal
        type(filter_sos_result) :: out
        out = zpk2sos_complex(cmplx(z,0.0_dp,dp),cmplx(p,0.0_dp,dp), &
            cmplx(k,0.0_dp,dp),pairing)
    end function zpk2sos_real_roots

    function zpk2sos_real_roots_complex_gain(z, p, k, pairing) result(out)
        real(dp), intent(in) :: z(:) !! real zeros
        real(dp), intent(in) :: p(:) !! real poles
        complex(dp), intent(in) :: k !! complex gain
        character(*), intent(in), optional :: pairing !! nearest (default), keep_odd, or minimal
        type(filter_sos_result) :: out
        out = zpk2sos_complex(cmplx(z,0.0_dp,dp),cmplx(p,0.0_dp,dp),k,pairing)
    end function zpk2sos_real_roots_complex_gain

    function zpk2sos_real_gain(z, p, k, pairing) result(out)
        complex(dp), intent(in) :: z(:) !! zeros
        complex(dp), intent(in) :: p(:) !! poles
        real(dp), intent(in) :: k !! real gain
        character(*), intent(in), optional :: pairing !! nearest (default), keep_odd, or minimal
        type(filter_sos_result) :: out
        out = zpk2sos_complex(z,p,cmplx(k,0.0_dp,dp),pairing)
    end function zpk2sos_real_gain

    function zpk2sos_complex(z, p, k, pairing) result(out)
        complex(dp), intent(in) :: z(:) !! complex zeros
        complex(dp), intent(in) :: p(:) !! complex poles
        complex(dp), intent(in) :: k !! complex gain
        character(*), intent(in), optional :: pairing !! nearest (default), keep_odd, or minimal
        type(filter_sos_result) :: out
        complex(dp), allocatable :: zz(:), pp(:), zg(:,:), pg(:,:)
        complex(dp), allocatable :: zpair(:), ppair(:)
        logical, allocatable :: assigned(:)
        complex(dp) :: zc, pc
        real(dp) :: cost, best
        logical :: pole_origin, zero_origin
        integer :: n, ns, i, j, match, extra_z, extra_p, count_origin
        character(:), allocatable :: mode

        if (.not. finite_complex(z) .or. .not. finite_complex(p)) return
        if (.not. ieee_is_finite(real(k,dp)) .or. .not. ieee_is_finite(aimag(k))) return
        mode = 'nearest'
        if (present(pairing)) mode = trim(pairing)
        if (mode /= 'nearest' .and. mode /= 'keep_odd' .and. mode /= 'minimal') return
        ! Minimal pairing cannot represent an improper polynomial ratio as
        ! a cascade of sections without introducing additional powers.
        if (mode == 'minimal' .and. size(z) > size(p)) return
        n = max(size(z),size(p))
        if (n > 120) return
        if (n == 0) then
            allocate(out%sos(1,6))
            out%sos(1,:) = cmplx([1.0_dp,0.0_dp,0.0_dp,1.0_dp,0.0_dp,0.0_dp],0.0_dp,dp)
            out%sos(1,1) = k
            out%status = zpk_conversion_success
            return
        end if
        ! SciPy's default digital pairing adds origin roots to balance
        ! numerator/denominator degrees and to make the order even.
        if (mod(n,2) /= 0) n = n+1
        ns = n/2
        extra_z = n-size(z)
        extra_p = n-size(p)
        allocate(zz(n),pp(n),zg(2,ns),pg(2,ns),assigned(ns))
        zz = (0.0_dp,0.0_dp)
        pp = (0.0_dp,0.0_dp)
        zz(:size(z)) = z
        pp(:size(p)) = p
        call group_roots(pp,pg)
        call group_roots(zz,zg)
        allocate(out%sos(ns,6))
        out%sos = (0.0_dp,0.0_dp)
        assigned = .false.
        do i = 1, ns
            pc = sum(pg(:,i))/2.0_dp
            best = huge(1.0_dp)
            match = 0
            do j = 1, ns
                if (assigned(j)) cycle
                zc = sum(zg(:,j))/2.0_dp
                cost = abs(zc-pc)
                if (mode == 'keep_odd' .and. mod(max(size(z),size(p)),2) /= 0) then
                    ! Keep the odd-order first-order factors in the same row.
                    pole_origin = any(abs(pg(:,i)) == 0.0_dp)
                    zero_origin = any(abs(zg(:,j)) == 0.0_dp)
                    if (pole_origin .neqv. zero_origin) cycle
                end if
                if (cost < best) then
                    best = cost
                    match = j
                end if
            end do
            if (match == 0) then
                ! With unequal zero/pole orders, keep_odd can exhaust
                ! sections having the same origin padding. In that case
                ! fall back to the nearest unassigned pair, as a valid
                ! factorization is preferable to an artificial failure.
                best = huge(1.0_dp)
                do j = 1, ns
                    if (assigned(j)) cycle
                    cost = abs(sum(zg(:,j))/2.0_dp - pc)
                    if (cost < best) then
                        best = cost
                        match = j
                    end if
                end do
                if (match == 0) return
            end if
            assigned(match) = .true.
            zpair = zg(:,match)
            ppair = pg(:,i)
            out%sos(i,1:3) = [cmplx(1.0_dp,0.0_dp,dp),-sum(zpair),zpair(1)*zpair(2)]
            out%sos(i,4:6) = [cmplx(1.0_dp,0.0_dp,dp),-sum(ppair),ppair(1)*ppair(2)]
            if (mode == 'minimal') then
                ! The temporary origin roots were inserted only to allow
                ! two-at-a-time pairing. Drop them by shifting polynomial
                ! coefficients to the right, preserving the original roots.
                ! This deliberately allows a0=0 in a minimal SOS section.
                count_origin = min(extra_z,count(zpair == (0.0_dp,0.0_dp)))
                do j = 1, count_origin
                    out%sos(i,1:3) = [cmplx(0.0_dp,0.0_dp,dp),out%sos(i,1:2)]
                end do
                extra_z = extra_z-count_origin
                count_origin = min(extra_p,count(ppair == (0.0_dp,0.0_dp)))
                do j = 1, count_origin
                    out%sos(i,4:6) = [cmplx(0.0_dp,0.0_dp,dp),out%sos(i,4:5)]
                end do
                extra_p = extra_p-count_origin
            end if
        end do
        if (mode == 'minimal' .and. (extra_z /= 0 .or. extra_p /= 0)) then
            out%status = zpk_conversion_numerical_failure
            return
        end if
        ! Distribute overall gain into the first section, matching SciPy.
        out%sos(1,1:3) = k*out%sos(1,1:3)
        if (.not. all(ieee_is_finite(real(out%sos,dp))) .or. &
                .not. all(ieee_is_finite(aimag(out%sos)))) then
            out%status = zpk_conversion_numerical_failure
            return
        end if
        out%status = zpk_conversion_success
    end function zpk2sos_complex

    function sos2zpk_real(sos) result(out)
        real(dp), intent(in) :: sos(:,:) !! SOS coefficients (nsections,6)
        type(filter_zpk_result) :: out
        out = sos2zpk_complex(cmplx(sos,0.0_dp,dp))
    end function sos2zpk_real

    function sos2zpk_complex(sos) result(out)
        complex(dp), intent(in) :: sos(:,:) !! SOS coefficients (nsections,6)
        type(filter_zpk_result) :: out
        type(filter_zpk_result) :: section
        integer :: ns, i, nz, np

        if (size(sos,1) == 0 .or. size(sos,2) /= 6) return
        if (.not. all(ieee_is_finite(real(sos,dp))) .or. &
                .not. all(ieee_is_finite(aimag(sos)))) return
        ns = size(sos,1)
        allocate(out%zeros(2*ns),out%poles(2*ns))
        out%zeros = cmplx(0.0_dp,0.0_dp,dp)
        out%poles = cmplx(0.0_dp,0.0_dp,dp)
        out%gain = cmplx(1.0_dp,0.0_dp,dp)
        do i = 1, ns
            section = tf2zpk_complex(sos(i,1:3),sos(i,4:6))
            if (section%status /= zpk_conversion_success) then
                out%status = section%status
                return
            end if
            nz = size(section%zeros)
            np = size(section%poles)
            out%zeros(2*i-1:2*i-2+nz) = section%zeros
            out%poles(2*i-1:2*i-2+np) = section%poles
            out%gain = out%gain*section%gain
            if (.not. ieee_is_finite(real(out%gain,dp)) .or. &
                    .not. ieee_is_finite(aimag(out%gain))) then
                out%status = zpk_conversion_numerical_failure
                return
            end if
        end do
        ! Like scipy.signal.sos2zpk, the result always includes exactly
        ! two zeros and two poles per section; absent roots are reported
        ! as zeros at the origin (also for minimal sections).
        out%status = zpk_conversion_success
    end function sos2zpk_complex

    subroutine group_roots(roots, pairs)
        complex(dp), intent(in) :: roots(:) !! even-sized list of roots
        complex(dp), intent(out) :: pairs(:,:) !! (2,number of sections)
        logical :: used(size(roots))
        integer :: i, j, pick, ns
        real(dp) :: best, d, tol
        complex(dp) :: r
        used = .false.
        ns = size(pairs,2)
        do i = 1, ns
            pick = 0
            do j = 1, size(roots)
                if (used(j)) cycle
                pick = j
                exit
            end do
            r = roots(pick)
            pairs(1,i) = r
            used(pick) = .true.
            best = huge(1.0_dp)
            pick = 0
            tol = 2.0e-8_dp*max(1.0_dp,abs(r))
            do j = 1, size(roots)
                if (used(j)) cycle
                ! Prefer conjugates for complex filters and real/real pairs.
                if (abs(aimag(r)) <= tol) then
                    if (abs(aimag(roots(j))) > tol) cycle
                    d = abs(roots(j)-r)
                else
                    d = abs(roots(j)-conjg(r))
                end if
                if (d < best) then
                    best = d
                    pick = j
                end if
            end do
            if (pick == 0) then
                do j = 1, size(roots)
                    if (used(j)) cycle
                    pick = j
                    exit
                end do
            end if
            pairs(2,i) = roots(pick)
            used(pick) = .true.
        end do
    end subroutine group_roots

end module scifort_signal_zpk_conversion
