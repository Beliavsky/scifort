! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Analog transfer-function and zero-pole-gain frequency responses.
! Evaluates H(i*w), with TF coefficients in descending powers of s.
module scifort_signal_analog
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal_design, only : frequency_response, design_status_success, &
        design_status_invalid_input, design_status_singular
    implicit none
    private

    real(dp), parameter :: pi = acos(-1.0_dp)

    interface freqs
        module procedure freqs_real
        module procedure freqs_complex
    end interface freqs
    interface freqs_zpk
        module procedure freqs_zpk_real_gain
        module procedure freqs_zpk_complex_gain
        module procedure freqs_zpk_real_roots
        module procedure freqs_zpk_real_roots_complex_gain
        module procedure freqs_zpk_real_zeros_real_gain
        module procedure freqs_zpk_real_zeros_complex_gain
        module procedure freqs_zpk_real_poles_real_gain
        module procedure freqs_zpk_real_poles_complex_gain
    end interface freqs_zpk
    interface findfreqs
        module procedure findfreqs_real
        module procedure findfreqs_complex
    end interface findfreqs
    public :: freqs, freqs_zpk, findfreqs, polynomial_roots

contains

    function findfreqs_real(num, den, n, kind, status) result(w)
        real(dp), intent(in) :: num(:) !! numerator polynomial or real zeros
        real(dp), intent(in) :: den(:) !! denominator polynomial or real poles
        integer, intent(in) :: n !! number of frequencies
        character(*), intent(in), optional :: kind !! 'ba' (default) or 'zp'
        integer, intent(out), optional :: status !! success or invalid input
        real(dp), allocatable :: w(:)

        w = findfreqs_complex(cmplx(num,0.0_dp,dp),cmplx(den,0.0_dp,dp),n,kind,status)
    end function findfreqs_real

    function findfreqs_complex(num, den, n, kind, status) result(w)
        complex(dp), intent(in) :: num(:) !! numerator polynomial or zeros
        complex(dp), intent(in) :: den(:) !! denominator polynomial or poles
        integer, intent(in) :: n !! number of logarithmic frequencies
        character(*), intent(in), optional :: kind !! 'ba' or 'zp'
        integer, intent(out), optional :: status !! zero on success
        real(dp), allocatable :: w(:)
        complex(dp), allocatable :: z(:), p(:)
        character(:), allocatable :: form
        logical :: ok

        allocate(w(0))
        if (present(status)) status = design_status_invalid_input
        if (n < 1 .or. .not. valid_complex(num) .or. .not. valid_complex(den)) return
        form = 'ba'
        if (present(kind)) form = trim(kind)
        select case (form)
        case ('ba')
            if (size(num) < 1 .or. size(den) < 1) return
            if (all(abs(den) == 0.0_dp)) return
            call polynomial_roots(num,z,ok)
            if (.not. ok) return
            call polynomial_roots(den,p,ok)
            if (.not. ok) return
        case ('zp')
            z = num
            p = den
        case default
            return
        end select
        call choose_frequencies(z,p,n,w,ok)
        if (.not. ok) then
            if (allocated(w)) deallocate(w)
            allocate(w(0))
            return
        end if
        if (present(status)) status = design_status_success
    end function findfreqs_complex

    pure logical function valid_complex(x) result(ok)
        complex(dp), intent(in) :: x(:) !! finite complex-valued coefficients or roots
        ok = all(ieee_is_finite(real(x,dp))) .and. all(ieee_is_finite(aimag(x)))
    end function valid_complex

    pure complex(dp) function eval_poly(a, s) result(y)
        complex(dp), intent(in) :: a(:) !! descending polynomial coefficients
        complex(dp), intent(in) :: s !! polynomial argument
        integer :: i
        y = a(1)
        do i = 2, size(a)
            y = y*s + a(i)
        end do
    end function eval_poly

    function freqs_real(b, a, worn, w) result(out)
        real(dp), intent(in) :: b(:) !! numerator coefficients, highest power first
        real(dp), intent(in), optional :: a(:) !! denominator coefficients, default [1]
        integer, intent(in), optional :: worn !! log-spaced frequency count, default 200
        real(dp), intent(in), optional :: w(:) !! explicit angular frequencies in rad/s
        type(frequency_response) :: out
        complex(dp), allocatable :: ac(:)

        ac = [cmplx(1.0_dp,0.0_dp,dp)]
        if (present(a)) ac = cmplx(a,0.0_dp,dp)
        out = freqs_core(cmplx(b,0.0_dp,dp),ac,worn,w)
    end function freqs_real

    function freqs_complex(b, a, worn, w) result(out)
        complex(dp), intent(in) :: b(:) !! numerator coefficients, highest power first
        complex(dp), intent(in), optional :: a(:) !! denominator coefficients, default [1]
        integer, intent(in), optional :: worn !! log-spaced frequency count, default 200
        real(dp), intent(in), optional :: w(:) !! explicit angular frequencies in rad/s
        type(frequency_response) :: out
        complex(dp), allocatable :: ac(:)

        ac = [cmplx(1.0_dp,0.0_dp,dp)]
        if (present(a)) ac = a
        out = freqs_core(b,ac,worn,w)
    end function freqs_complex

    function freqs_core(b, a, worn, w) result(out)
        complex(dp), intent(in) :: b(:) !! numerator polynomial in descending powers
        complex(dp), intent(in) :: a(:) !! denominator polynomial in descending powers
        integer, intent(in), optional :: worn !! number of automatically chosen frequencies
        real(dp), intent(in), optional :: w(:) !! explicit angular frequencies
        type(frequency_response) :: out
        complex(dp), allocatable :: z(:), p(:)
        complex(dp) :: numerator, denominator, s
        integer :: i, n
        logical :: ok

        if (size(b) == 0 .or. size(a) == 0) return
        if (.not. valid_complex(b) .or. .not. valid_complex(a)) return
        if (all(abs(a) == 0.0_dp)) return
        if (present(w)) then
            if (.not. all(ieee_is_finite(w))) return
            out%frequencies = w
        else
            n = 200
            if (present(worn)) n = worn
            if (n < 1) return
            call polynomial_roots(b,z,ok)
            if (.not. ok) return
            call polynomial_roots(a,p,ok)
            if (.not. ok) return
            call choose_frequencies(z,p,n,out%frequencies,ok)
            if (.not. ok) return
        end if

        allocate(out%response(size(out%frequencies)))
        out%status = design_status_success
        do i = 1, size(out%frequencies)
            s = cmplx(0.0_dp,out%frequencies(i),dp)
            numerator = eval_poly(b,s)
            denominator = eval_poly(a,s)
            if (abs(denominator) <= tiny(1.0_dp)) then
                out%response(i) = cmplx(ieee_value(0.0_dp,ieee_quiet_nan), &
                    ieee_value(0.0_dp,ieee_quiet_nan),dp)
                out%status = design_status_singular
            else
                out%response(i) = numerator/denominator
            end if
        end do
    end function freqs_core

    function freqs_zpk_real_gain(z, p, k, worn, w) result(out)
        complex(dp), intent(in) :: z(:) !! zeros of the analog transfer function
        complex(dp), intent(in) :: p(:) !! poles of the analog transfer function
        real(dp), intent(in) :: k !! scalar gain
        integer, intent(in), optional :: worn !! log-spaced frequency count, default 200
        real(dp), intent(in), optional :: w(:) !! explicitly specified angular frequencies
        type(frequency_response) :: out

        out = freqs_zpk_core(z,p,cmplx(k,0.0_dp,dp),worn,w)
    end function freqs_zpk_real_gain

    function freqs_zpk_complex_gain(z, p, k, worn, w) result(out)
        complex(dp), intent(in) :: z(:) !! zeros of the analog transfer function
        complex(dp), intent(in) :: p(:) !! poles of the analog transfer function
        complex(dp), intent(in) :: k !! complex-valued scalar gain
        integer, intent(in), optional :: worn !! log-spaced frequency count, default 200
        real(dp), intent(in), optional :: w(:) !! explicitly specified angular frequencies
        type(frequency_response) :: out

        out = freqs_zpk_core(z,p,k,worn,w)
    end function freqs_zpk_complex_gain

    function freqs_zpk_real_roots(z, p, k, worn, w) result(out)
        real(dp), intent(in) :: z(:) !! real zeros
        real(dp), intent(in) :: p(:) !! real poles
        real(dp), intent(in) :: k !! real scalar gain
        integer, intent(in), optional :: worn !! log-spaced frequency count, default 200
        real(dp), intent(in), optional :: w(:) !! explicitly specified angular frequencies
        type(frequency_response) :: out

        out = freqs_zpk_core(cmplx(z,0.0_dp,dp),cmplx(p,0.0_dp,dp), &
            cmplx(k,0.0_dp,dp),worn,w)
    end function freqs_zpk_real_roots

    function freqs_zpk_real_roots_complex_gain(z, p, k, worn, w) result(out)
        real(dp), intent(in) :: z(:) !! real zeros
        real(dp), intent(in) :: p(:) !! real poles
        complex(dp), intent(in) :: k !! complex-valued scalar gain
        integer, intent(in), optional :: worn !! log-spaced frequency count, default 200
        real(dp), intent(in), optional :: w(:) !! explicitly specified angular frequencies
        type(frequency_response) :: out

        out = freqs_zpk_core(cmplx(z,0.0_dp,dp),cmplx(p,0.0_dp,dp),k,worn,w)
    end function freqs_zpk_real_roots_complex_gain

    function freqs_zpk_real_zeros_real_gain(z, p, k, worn, w) result(out)
        real(dp), intent(in) :: z(:) !! filter zero locations
        complex(dp), intent(in) :: p(:) !! filter pole locations
        real(dp), intent(in) :: k !! scalar gain
        integer, intent(in), optional :: worn !! logarithmic grid count
        real(dp), intent(in), optional :: w(:) !! explicit angular frequencies
        type(frequency_response) :: out

        out = freqs_zpk_core(cmplx(z,0.0_dp,dp),p,cmplx(k,0.0_dp,dp),worn,w)
    end function freqs_zpk_real_zeros_real_gain

    function freqs_zpk_real_zeros_complex_gain(z, p, k, worn, w) result(out)
        real(dp), intent(in) :: z(:) !! filter zero locations
        complex(dp), intent(in) :: p(:) !! filter pole locations
        complex(dp), intent(in) :: k !! scalar gain
        integer, intent(in), optional :: worn !! logarithmic grid count
        real(dp), intent(in), optional :: w(:) !! explicit angular frequencies
        type(frequency_response) :: out

        out = freqs_zpk_core(cmplx(z,0.0_dp,dp),p,k,worn,w)
    end function freqs_zpk_real_zeros_complex_gain

    function freqs_zpk_real_poles_real_gain(z, p, k, worn, w) result(out)
        complex(dp), intent(in) :: z(:) !! filter zero locations
        real(dp), intent(in) :: p(:) !! filter pole locations
        real(dp), intent(in) :: k !! scalar gain
        integer, intent(in), optional :: worn !! logarithmic grid count
        real(dp), intent(in), optional :: w(:) !! explicit angular frequencies
        type(frequency_response) :: out

        out = freqs_zpk_core(z,cmplx(p,0.0_dp,dp),cmplx(k,0.0_dp,dp),worn,w)
    end function freqs_zpk_real_poles_real_gain

    function freqs_zpk_real_poles_complex_gain(z, p, k, worn, w) result(out)
        complex(dp), intent(in) :: z(:) !! filter zero locations
        real(dp), intent(in) :: p(:) !! filter pole locations
        complex(dp), intent(in) :: k !! scalar gain
        integer, intent(in), optional :: worn !! logarithmic grid count
        real(dp), intent(in), optional :: w(:) !! explicit angular frequencies
        type(frequency_response) :: out

        out = freqs_zpk_core(z,cmplx(p,0.0_dp,dp),k,worn,w)
    end function freqs_zpk_real_poles_complex_gain

    function freqs_zpk_core(z, p, k, worn, w) result(out)
        complex(dp), intent(in) :: z(:) !! zeros of numerator transfer polynomial
        complex(dp), intent(in) :: p(:) !! poles of denominator transfer polynomial
        complex(dp), intent(in) :: k !! scalar transfer-function gain
        integer, intent(in), optional :: worn !! number of automatically chosen frequencies
        real(dp), intent(in), optional :: w(:) !! explicit angular frequencies
        type(frequency_response) :: out
        complex(dp) :: numerator, denominator, s
        integer :: i, j, n
        logical :: ok

        if (.not. valid_complex(z) .or. .not. valid_complex(p)) return
        if (.not. valid_complex([k])) return
        if (present(w)) then
            if (.not. all(ieee_is_finite(w))) return
            out%frequencies = w
        else
            n = 200
            if (present(worn)) n = worn
            if (n < 1) return
            call choose_frequencies(z,p,n,out%frequencies,ok)
            if (.not. ok) return
        end if
        allocate(out%response(size(out%frequencies)))
        out%status = design_status_success
        do i = 1, size(out%frequencies)
            s = cmplx(0.0_dp,out%frequencies(i),dp)
            numerator = k
            denominator = cmplx(1.0_dp,0.0_dp,dp)
            do j = 1, size(z)
                numerator = numerator*(s-z(j))
            end do
            do j = 1, size(p)
                denominator = denominator*(s-p(j))
            end do
            if (abs(denominator) <= tiny(1.0_dp)) then
                out%response(i) = cmplx(ieee_value(0.0_dp,ieee_quiet_nan), &
                    ieee_value(0.0_dp,ieee_quiet_nan),dp)
                out%status = design_status_singular
            else
                out%response(i) = numerator/denominator
            end if
        end do
    end function freqs_zpk_core

    subroutine choose_frequencies(z, p, n, frequencies, ok)
        complex(dp), intent(in) :: z(:) !! zeros defining interesting frequency regions
        complex(dp), intent(in) :: p(:) !! poles defining interesting frequency regions
        integer, intent(in) :: n !! number of logarithmically spaced frequencies
        real(dp), allocatable, intent(out) :: frequencies(:) !! selected angular frequencies
        logical, intent(out) :: ok !! success if frequency grid can be determined
        complex(dp), allocatable :: chosen(:), poles(:)
        complex(dp) :: root
        real(dp) :: low, high, vlo, vhi, inte, re, im, mag
        integer :: i, nsel

        ok = .false.
        allocate(frequencies(0))
        if (n < 1) return
        poles = p
        if (size(p) == 0) poles = [cmplx(-1000.0_dp,0.0_dp,dp)]
        allocate(chosen(size(poles)+size(z)))
        nsel = 0
        do i = 1, size(poles)+size(z)
            if (i <= size(poles)) then
                root = poles(i)
            else
                root = z(i-size(poles))
                if (abs(root) >= 1.0e5_dp) cycle
            end if
            if (aimag(root) < -1.0e-10_dp*(1.0_dp+abs(root))) cycle
            nsel = nsel+1
            chosen(nsel) = root
        end do
        if (nsel == 0) return
        vlo = huge(1.0_dp)
        vhi = 0.0_dp
        do i = 1, nsel
            root = chosen(i)
            re = real(root,dp)
            im = max(0.0_dp,aimag(root))
            mag = abs(root)
            inte = 0.0_dp
            if (mag < 1.0e-10_dp) inte = 1.0_dp
            vlo = min(vlo, 0.1_dp*(abs(re+inte)+2.0_dp*im))
            vhi = max(vhi, 3.0_dp*abs(re+inte)+1.5_dp*im)
        end do
        if (vlo <= 0.0_dp .or. vhi <= 0.0_dp) return
        low = nearest_even(log10(vlo)-0.5_dp)
        high = nearest_even(log10(vhi)+0.5_dp)
        deallocate(frequencies)
        allocate(frequencies(n))
        if (n == 1) then
            frequencies(1) = 10.0_dp**low
        else
            do i = 1, n
                frequencies(i) = 10.0_dp**(low+(high-low)*real(i-1,dp)/real(n-1,dp))
            end do
        end if
        ok = all(ieee_is_finite(frequencies))
    end subroutine choose_frequencies

    pure real(dp) function nearest_even(x) result(y)
        real(dp), intent(in) :: x !! logarithmic decade exponent
        integer :: low
        real(dp) :: delta

        low = floor(x)
        delta = x-real(low,dp)
        if (delta < 0.5_dp) then
            y = real(low,dp)
        else if (delta > 0.5_dp) then
            y = real(low+1,dp)
        else if (mod(abs(low),2) == 0) then
            y = real(low,dp)
        else
            y = real(low+1,dp)
        end if
    end function nearest_even

    subroutine polynomial_roots(c, roots, ok)
        complex(dp), intent(in) :: c(:) !! descending polynomial coefficients
        complex(dp), allocatable, intent(out) :: roots(:) !! estimated polynomial roots
        logical, intent(out) :: ok !! true if all roots were resolved
        complex(dp), allocatable :: q(:), previous(:), active(:)
        complex(dp) :: f, df, step, interaction, y
        real(dp) :: scale, value, angle, largest_step, bound
        integer :: istart, n, nzero, i, j, it, degree

        ok = .false.
        allocate(roots(0))
        if (size(c) == 0) return
        istart = 1
        do while (istart <= size(c))
            if (abs(c(istart)) /= 0.0_dp) exit
            istart = istart+1
        end do
        if (istart > size(c)) then
            ok = .true.
            return
        end if
        degree = size(c)-istart
        deallocate(roots)
        allocate(roots(degree))
        if (degree == 0) then
            ok = .true.
            return
        end if
        nzero = 0
        do i = size(c), istart+1, -1
            if (abs(c(i)) /= 0.0_dp) exit
            nzero = nzero+1
        end do
        roots = cmplx(0.0_dp,0.0_dp,dp)
        n = degree-nzero
        if (n == 0) then
            ok = .true.
            return
        end if
        if (n == 1) then
            roots(1) = -c(istart+1)/c(istart)
            ok = valid_complex(roots)
            return
        end if
        scale = tiny(1.0_dp)
        do i = 1, n
            value = abs(c(istart+i)/c(istart))
            if (value > 0.0_dp) scale = max(scale,exp(log(value)/real(i,dp)))
        end do
        if (.not. ieee_is_finite(scale) .or. scale == 0.0_dp) return
        allocate(q(n+1),previous(n),active(n))
        q(1) = cmplx(1.0_dp,0.0_dp,dp)
        do i = 1, n
            q(i+1) = (c(istart+i)/c(istart))/scale**i
        end do
        if (.not. valid_complex(q)) return
        do i = 1, n
            angle = 2.0_dp*pi*(real(i-1,dp)+0.273_dp)/real(n,dp)
            active(i) = 0.87_dp*cmplx(cos(angle),sin(angle),dp)
        end do
        do it = 1, 550
            previous = active
            largest_step = 0.0_dp
            do i = 1, n
                y = previous(i)
                f = q(1)
                df = cmplx(0.0_dp,0.0_dp,dp)
                do j = 2, n+1
                    df = df*y+f
                    f = f*y+q(j)
                end do
                if (abs(df) < tiny(1.0_dp)) then
                    largest_step = huge(1.0_dp)
                    cycle
                end if
                step = f/df
                interaction = cmplx(0.0_dp,0.0_dp,dp)
                do j = 1, n
                    if (j == i) cycle
                    if (abs(y-previous(j)) < tiny(1.0_dp)) cycle
                    interaction = interaction + 1.0_dp/(y-previous(j))
                end do
                if (abs(1.0_dp-step*interaction) > tiny(1.0_dp)) &
                    step = step/(1.0_dp-step*interaction)
                active(i) = y-step
                largest_step = max(largest_step,abs(step)/(1.0_dp+abs(y)))
            end do
            if (.not. valid_complex(active)) return
            if (largest_step < 2.0e-13_dp) exit
        end do
        if (it > 550) then
            ! Repeated roots can converge slowly. For selecting the frequency
            ! band, a small relative polynomial residual is sufficient.
            do i = 1, n
                y = active(i)
                f = q(1)
                bound = abs(q(1))
                do j = 2, n+1
                    f = f*y+q(j)
                    bound = bound*abs(y)+abs(q(j))
                end do
                if (abs(f) > 1.0e-9_dp*max(bound,tiny(1.0_dp))) return
            end do
        end if
        roots(1:n) = scale*active
        ok = valid_complex(roots)
    end subroutine polynomial_roots

end module scifort_signal_analog
