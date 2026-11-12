! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Additional scipy.signal.windows window families. Elementary windows are
! independently implemented from mathematical definitions. The Chebyshev,
! Taylor, and Kaiser-Bessel-derived constructions adapt the SciPy 1.17.0
! scipy/signal/windows/_windows.py algorithms (BSD-3-Clause); SciFort uses a
! log-scaled spectrum and direct DFT for Chebyshev. See CODE_PROVENANCE.md
! and THIRD_PARTY_LICENSES.md for attribution and the full license.
module scifort_signal_windows_extra
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    implicit none
    private
    real(dp), parameter :: pi = acos(-1.0_dp)
    integer, parameter, public :: window_status_success = 0
    integer, parameter, public :: window_status_invalid_input = 1
    public :: window_samples, chebwin, taylor, kaiser_bessel_derived, general_cosine
contains

    function window_samples(name, n, fftbins, parameter, parameter2, status) result(w)
        character(*), intent(in) :: name !! SciPy-style window name
        integer, intent(in) :: n !! requested sample count
        logical, intent(in), optional :: fftbins !! true for periodic, default true
        real(dp), intent(in), optional :: parameter !! shape: beta/std/alpha/tau/p
        real(dp), intent(in), optional :: parameter2 !! general Gaussian width or exponential center
        integer, intent(out), optional :: status !! zero on success; 1 if invalid
        real(dp), allocatable :: w(:)
        real(dp) :: x, a, b, q, h, xx, f, arg
        integer :: m, i, k, nbar, code
        logical :: periodic
        character(:), allocatable :: kind

        allocate(w(0))
        if (present(status)) status = window_status_invalid_input
        if (n < 0 .or. n == huge(n)) return
        if (present(parameter)) then
            if (.not. ieee_is_finite(parameter)) return
        end if
        if (present(parameter2)) then
            if (.not. ieee_is_finite(parameter2)) return
        end if
        kind = trim(name)
        select case (kind)
        case ('boxcar', 'rectangular', 'ones', 'rect', 'box', &
              'hann', 'hanning', 'han', 'hamming', 'hamm', 'ham', &
              'blackman', 'black', 'blk', 'bartlett', 'bart', 'brt', &
              'tukey', 'tuk', 'flattop', 'flat', 'flt', &
              'triang', 'triangle', 'tri', 'parzen', 'parz', 'par', &
              'bohman', 'bman', 'bmn', 'blackmanharris', 'blackharr', 'bkh', &
              'nuttall', 'nutl', 'nut', 'barthann', 'brthan', 'bth', &
              'cosine', 'halfcosine', 'lanczos', 'sinc', &
              'exponential', 'poisson', 'gaussian', 'gauss', 'gss', &
              'kaiser', 'ksr', 'general_gaussian', 'general gauss', &
              'general_gauss', 'ggs', 'general_hamming', 'general hamming', &
              'chebwin', 'cheb', 'taylor', 'taylorwin', &
              'kaiser_bessel_derived', 'kaiser bessel derived', 'kbd')
        case default
            return
        end select
        select case (kind)
        case ('kaiser','ksr','gaussian','gauss','gss', &
              'general_gaussian','general gauss','general_gauss','ggs', &
              'general_hamming','general hamming')
            if (.not. present(parameter)) return
        end select
        select case (kind)
        case ('general_gaussian','general gauss','general_gauss','ggs')
            if (.not. present(parameter2)) return
            if (parameter <= 0.0_dp .or. parameter2 <= 0.0_dp) return
        case ('gaussian','gauss','gss','exponential','poisson')
            if (present(parameter)) then
                if (parameter <= 0.0_dp) return
            end if
        end select
        periodic = .true.
        if (present(fftbins)) periodic = fftbins
        ! scipy.signal.windows.exponential allows an explicitly displaced
        ! center only in periodic mode. In SciFort parameter=tau and
        ! parameter2=center; the existing single-tau API is unchanged.
        if (present(parameter2)) then
            select case (kind)
            case ('exponential','poisson')
                if (.not. periodic) return
            end select
        end if
        ! Some windows have non-cosine-series definitions and need their own
        ! normalization, or (for KBD) reject periodic output altogether.
        select case (kind)
        case ('chebwin','cheb')
            if (.not. present(parameter)) return
            w = chebwin(n, parameter, sym=.not.periodic, status=code)
            if (present(status)) status = code
            return
        case ('taylor','taylorwin')
            nbar = 4
            if (present(parameter)) then
                if (parameter < 1.0_dp .or. parameter > 1000.0_dp) return
                nbar = int(parameter)
                if (real(nbar,dp) /= parameter) return
            end if
            b = 30.0_dp
            if (present(parameter2)) b = parameter2
            w = taylor(n, nbar=nbar, sll=b, sym=.not.periodic, status=code)
            if (present(status)) status = code
            return
        case ('kaiser_bessel_derived','kaiser bessel derived','kbd')
            if (.not. present(parameter)) return
            w = kaiser_bessel_derived(n, parameter, sym=.not.periodic, status=code)
            if (present(status)) status = code
            return
        end select
        m = n
        if (periodic .and. n > 1) m = n + 1
        deallocate(w)
        allocate(w(n))
        if (n < 2) then
            w = 1.0_dp
            if (present(status)) status = window_status_success
            return
        end if
        a = 0.0_dp
        b = 0.0_dp
        if (present(parameter)) a = parameter
        if (present(parameter2)) b = parameter2
        if (kind == 'tukey' .or. kind == 'tuk') then
            a = 0.5_dp
            if (present(parameter)) a = parameter
        end if
        if (kind == 'exponential' .or. kind == 'poisson') then
            a = 1.0_dp
            if (present(parameter)) a = parameter
        end if
        h = real(m - 1, dp)
        do i = 1, n
            x = real(i - 1, dp) / h
            xx = real(i - 1, dp) - 0.5_dp * h
            select case (kind)
            case ('boxcar','rectangular','ones','rect','box')
                f = 1.0_dp
            case ('hann','hanning','han')
                f = 0.5_dp - 0.5_dp * cos(2.0_dp * pi * x)
            case ('hamming','hamm','ham')
                f = 0.54_dp - 0.46_dp * cos(2.0_dp * pi * x)
            case ('blackman','black','blk')
                f = cosine_sum(x, [0.42_dp, 0.5_dp, 0.08_dp])
            case ('bartlett','bart','brt')
                f = 1.0_dp - abs(2.0_dp * x - 1.0_dp)
            case ('flattop','flat','flt')
                f = cosine_sum(x, [0.21557895_dp,0.41663158_dp,0.277263158_dp, &
                                     0.083578947_dp,0.006947368_dp])
            case ('triang','triangle','tri')
                k = min(i, m + 1 - i)
                if (mod(m,2) == 0) then
                    f = real(2*k-1,dp) / real(m,dp)
                else
                    f = 2.0_dp * real(k,dp) / real(m+1,dp)
                end if
            case ('parzen','parz','par')
                q = abs(xx) / (real(m,dp)/2.0_dp)
                if (abs(xx) <= h / 4.0_dp) then
                    f = 1.0_dp - 6.0_dp*q*q + 6.0_dp*q*q*q
                else
                    f = 2.0_dp * (1.0_dp-q)**3
                end if
            case ('bohman','bman','bmn')
                q = abs(2.0_dp*x-1.0_dp)
                if (i == 1 .or. i == m) then
                    f = 0.0_dp
                else
                    f = (1.0_dp-q)*cos(pi*q) + sin(pi*q)/pi
                end if
            case ('blackmanharris','blackharr','bkh')
                f = cosine_sum(x, [0.35875_dp,0.48829_dp,0.14128_dp,0.01168_dp])
            case ('nuttall','nutl','nut')
                f = cosine_sum(x, [0.3635819_dp,0.4891775_dp,0.1365995_dp,0.0106411_dp])
            case ('barthann','brthan','bth')
                q = abs(x - 0.5_dp)
                f = 0.62_dp - 0.48_dp*q + 0.38_dp*cos(2.0_dp*pi*q)
            case ('cosine','halfcosine')
                f = sin(pi*(real(i,dp)-0.5_dp)/real(m,dp))
            case ('lanczos','sinc')
                q = 2.0_dp*x - 1.0_dp
                if (abs(q) < 1.0e-14_dp) then
                    f = 1.0_dp
                else
                    f = sin(pi*q)/(pi*q)
                end if
            case ('tukey','tuk')
                if (a <= 0.0_dp) then
                    f = 1.0_dp
                else if (a >= 1.0_dp) then
                    f = 0.5_dp - 0.5_dp*cos(2.0_dp*pi*x)
                else if (x < a/2.0_dp) then
                    f = 0.5_dp*(1.0_dp+cos(pi*(2.0_dp*x/a-1.0_dp)))
                else if (x > 1.0_dp-a/2.0_dp) then
                    f = 0.5_dp*(1.0_dp+cos(pi*(2.0_dp*x/a-2.0_dp/a+1.0_dp)))
                else
                    f = 1.0_dp
                end if
            case ('exponential','poisson')
                if (present(parameter2)) then
                    f = exp(-abs(real(i-1,dp)-b)/a)
                else
                    f = exp(-abs(xx)/a)
                end if
            case ('gaussian','gauss','gss')
                f = exp(-0.5_dp*(xx/a)**2)
            case ('kaiser','ksr')
                arg = abs(a)*sqrt(max(0.0_dp,1.0_dp-(2.0_dp*x-1.0_dp)**2))
                f = exp(log_besseli0(arg)-log_besseli0(abs(a)))
            case ('general_gaussian','general gauss','general_gauss','ggs')
                f = exp(-0.5_dp*abs(xx/b)**(2.0_dp*a))
            case ('general_hamming','general hamming')
                f = a - (1.0_dp-a)*cos(2.0_dp*pi*x)
            end select
            w(i) = f
        end do
        if (present(status)) status = window_status_success
    end function window_samples

    function general_cosine(n, coefficients, sym, status) result(w)
        integer, intent(in) :: n !! requested window length
        real(dp), intent(in) :: coefficients(:) !! cosine-series weights centered at zero
        logical, intent(in), optional :: sym !! symmetric if true, default true
        integer, intent(out), optional :: status !! zero on success
        real(dp), allocatable :: w(:)
        integer :: i, m
        logical :: symmetric

        allocate(w(0))
        if (present(status)) status = window_status_invalid_input
        if (n < 0 .or. n == huge(n)) return
        if (any(.not. ieee_is_finite(coefficients))) return
        if (n == 0) then
            if (present(status)) status = window_status_success
            return
        end if
        deallocate(w)
        allocate(w(n))
        if (n == 1) then
            ! The short-window convention matches scipy.signal.windows.
            w = 1.0_dp
        else if (size(coefficients) == 0) then
            w = 0.0_dp
        else
            symmetric = .true.
            if (present(sym)) symmetric = sym
            m = n
            if (.not. symmetric) m = n + 1
            do i = 1, n
                w(i) = cosine_sum(real(i-1,dp)/real(m-1,dp), coefficients)
            end do
        end if
        if (present(status)) status = window_status_success
    end function general_cosine

    function chebwin(n, at, sym, status) result(w)
        integer, intent(in) :: n !! window length
        real(dp), intent(in) :: at !! sidelobe attenuation in decibels
        logical, intent(in), optional :: sym !! symmetric if true; default true
        integer, intent(out), optional :: status !! zero on success
        real(dp), allocatable :: w(:)
        real(dp), allocatable :: p(:), u(:)
        real(dp) :: beta, x, degree, ref, t, val, offset, total
        integer :: m, k, i, j
        logical :: symmetric

        allocate(w(0))
        if (present(status)) status = window_status_invalid_input
        if (n < 0 .or. n == huge(n)) return
        if (.not. ieee_is_finite(at)) return
        if (abs(at) > 5000.0_dp) return
        symmetric = .true.
        if (present(sym)) symmetric = sym
        if (n <= 1) then
            deallocate(w)
            allocate(w(n))
            w = 1.0_dp
            if (present(status)) status = window_status_success
            return
        end if
        m = n
        if (.not. symmetric) m = n + 1
        degree = real(m - 1,dp)
        ! Normalizing the frequency samples before the DFT prevents overflow
        ! for large sidelobe attenuations without affecting the final window.
        beta = cosh(acosh(exp(log(10.0_dp)*abs(at)/20.0_dp))/degree)
        ref = log_cosh(degree*acosh(beta))
        allocate(p(0:m-1), u(0:m/2))
        do k = 0, m-1
            x = beta*cos(pi*real(k,dp)/real(m,dp))
            if (abs(x) > 1.0_dp) then
                val = exp(log_cosh(degree*acosh(abs(x))) - ref)
                if (x < -1.0_dp .and. mod(m,2)==0) val = -val
            else
                val = cos(degree*acos(max(-1.0_dp,min(1.0_dp,x)))) * exp(-ref)
            end if
            p(k) = val
        end do
        do j = 0, m/2
            offset = real(j,dp)
            if (mod(m,2)==0) offset = offset - 0.5_dp
            total = 0.0_dp
            do k = 0, m-1
                t = 2.0_dp*pi*offset*real(k,dp)/real(m,dp)
                total = total + p(k)*cos(t)
            end do
            u(j) = total
        end do
        deallocate(w)
        allocate(w(n))
        do i = 1, n
            if (mod(m,2)==0) then
                j = int(abs(real(i,dp)-0.5_dp*real(m+1,dp))+0.5_dp)
            else
                j = abs(i - (m+1)/2)
            end if
            w(i) = u(j)
        end do
        val = maxval(w)
        if (.not. ieee_is_finite(val) .or. val <= 0.0_dp) then
            deallocate(w)
            allocate(w(0))
            return
        end if
        w = w/val
        if (present(status)) status = window_status_success
    end function chebwin

    function taylor(n, nbar, sll, norm, sym, status) result(w)
        integer, intent(in) :: n !! window length
        integer, intent(in), optional :: nbar !! number of approximately constant sidelobes, default 4
        real(dp), intent(in), optional :: sll !! sidelobe suppression in dB, default 30
        logical, intent(in), optional :: norm !! normalize peak, default true
        logical, intent(in), optional :: sym !! symmetric if true; default true
        integer, intent(out), optional :: status !! zero on success
        real(dp), allocatable :: w(:), fm(:)
        real(dp) :: a, sl, s2, numer, denom, t, center
        integer :: m, nb, i, j, k
        logical :: symmetric, normalize

        allocate(w(0))
        if (present(status)) status = window_status_invalid_input
        if (n < 0 .or. n == huge(n)) return
        nb = 4
        if (present(nbar)) nb = nbar
        sl = 30.0_dp
        if (present(sll)) sl = sll
        if (nb < 1 .or. nb > 1000) return
        if (.not. ieee_is_finite(sl) .or. sl <= 0.0_dp .or. sl > 5000.0_dp) return
        symmetric = .true.
        if (present(sym)) symmetric = sym
        normalize = .true.
        if (present(norm)) normalize = norm
        if (n <= 1) then
            deallocate(w)
            allocate(w(n))
            w = 1.0_dp
            if (present(status)) status = window_status_success
            return
        end if
        m = n
        if (.not. symmetric) m = n + 1
        a = acosh(exp(log(10.0_dp)*sl/20.0_dp))/pi
        s2 = real(nb,dp)**2 / (a*a + (real(nb,dp)-0.5_dp)**2)
        allocate(fm(nb-1))
        do i = 1, nb-1
            numer = 1.0_dp
            denom = 2.0_dp
            do j = 1, nb-1
                numer = numer*(1.0_dp - real(i,dp)**2/ &
                    (s2*(a*a + (real(j,dp)-0.5_dp)**2)))
                if (i /= j) denom = denom*(1.0_dp - real(i,dp)**2/real(j,dp)**2)
            end do
            if (mod(i,2)==0) numer = -numer
            fm(i) = numer/denom
        end do
        if (any(.not. ieee_is_finite(fm))) return
        center = 1.0_dp + 2.0_dp*sum(fm)
        if (normalize .and. (center <= tiny(1.0_dp) .or. .not. ieee_is_finite(center))) return
        deallocate(w)
        allocate(w(n))
        do k = 0, n-1
            t = 1.0_dp
            do i = 1, nb-1
                t = t + 2.0_dp*fm(i)*cos(2.0_dp*pi*real(i,dp)* &
                    (real(k,dp)-0.5_dp*real(m-1,dp))/real(m,dp))
            end do
            if (normalize) t = t/center
            w(k+1) = t
        end do
        if (any(.not. ieee_is_finite(w))) then
            deallocate(w)
            allocate(w(0))
            return
        end if
        if (present(status)) status = window_status_success
    end function taylor

    function kaiser_bessel_derived(n, beta, sym, status) result(w)
        integer, intent(in) :: n !! even window length
        real(dp), intent(in) :: beta !! Kaiser shape parameter
        logical, intent(in), optional :: sym !! must be true (default)
        integer, intent(out), optional :: status !! zero on success
        real(dp), allocatable :: w(:)
        real(dp) :: a, x, s, total
        integer :: i, m

        allocate(w(0))
        if (present(status)) status = window_status_invalid_input
        if (n < 0 .or. .not. ieee_is_finite(beta)) return
        if (abs(beta) > 5000.0_dp) return
        if (present(sym)) then
            if (.not. sym) return
        end if
        if (n == 0) then
            if (present(status)) status = window_status_success
            return
        end if
        if (mod(n,2) /= 0) return
        m = n/2
        a = log_besseli0(abs(beta))
        total = 0.0_dp
        do i = 0, m
            x = 2.0_dp*real(i,dp)/real(m,dp) - 1.0_dp
            total = total + exp(log_besseli0(abs(beta)*sqrt(max(0.0_dp,1.0_dp-x*x))) - a)
        end do
        if (total <= 0.0_dp .or. .not. ieee_is_finite(total)) return
        deallocate(w)
        allocate(w(n))
        s = 0.0_dp
        do i = 0, m-1
            x = 2.0_dp*real(i,dp)/real(m,dp) - 1.0_dp
            s = s + exp(log_besseli0(abs(beta)*sqrt(max(0.0_dp,1.0_dp-x*x))) - a)
            w(i+1) = sqrt(s/total)
            w(n-i) = w(i+1)
        end do
        if (present(status)) status = window_status_success
    end function kaiser_bessel_derived

    pure function log_cosh(x) result(y)
        real(dp), intent(in) :: x !! nonnegative argument
        real(dp) :: y
        y = x - log(2.0_dp) + log(1.0_dp + exp(-2.0_dp*x))
    end function log_cosh

    pure function cosine_sum(x, coeff) result(y)
        real(dp), intent(in) :: x !! coordinate divided by (M-1)
        real(dp), intent(in) :: coeff(:) !! cosine-series coefficients
        real(dp) :: y
        integer :: j
        y = coeff(1)
        do j = 2, size(coeff)
            if (mod(j,2)==0) then
                y = y - coeff(j)*cos(2.0_dp*pi*real(j-1,dp)*x)
            else
                y = y + coeff(j)*cos(2.0_dp*pi*real(j-1,dp)*x)
            end if
        end do
    end function cosine_sum

    pure function log_besseli0(x) result(y)
        real(dp), intent(in) :: x !! absolute Bessel argument
        real(dp) :: y, term, total, t, z
        integer :: k
        if (x < 45.0_dp) then
            term = 1.0_dp
            total = 1.0_dp
            do k = 1, 400
                term = term * (0.5_dp*x/real(k,dp))**2
                total = total + term
                if (term <= epsilon(total)*total) exit
            end do
            y = log(total)
        else
            z = 1.0_dp/x
            t = 1.0_dp + z*(1.0_dp/8.0_dp + z*(9.0_dp/128.0_dp + &
                z*(225.0_dp/3072.0_dp + z*(11025.0_dp/98304.0_dp + &
                z*(893025.0_dp/3932160.0_dp + z*108056025.0_dp/188743680.0_dp)))))
            y = x - 0.5_dp*log(2.0_dp*pi*x) + log(t)
        end if
    end function log_besseli0
end module scifort_signal_windows_extra
