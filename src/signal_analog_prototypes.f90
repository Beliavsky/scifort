! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Normalized analog lowpass prototypes in zero-pole-gain form.
! Butterworth and Chebyshev poles follow the classical circle/ellipse
! constructions. Bessel poles are roots of a scaled reverse Bessel
! polynomial, refined simultaneously with the Aberth-Ehrlich iteration.
! The scaled half-integer Bessel K kernel avoids high-order polynomial cancellation.
! No external root-finding or linear-algebra dependency is needed.
module scifort_signal_analog_prototypes
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_bessel_complex, only : kve_complex
    implicit none
    private

    real(dp), parameter :: pi = acos(-1.0_dp)
    integer, parameter, public :: prototype_success = 0
    integer, parameter, public :: prototype_invalid_input = 1
    integer, parameter, public :: prototype_no_convergence = 2
    integer, parameter, public :: prototype_too_expensive = 3
    integer, parameter :: max_trigonometric_order = 64
    integer, parameter :: max_bessel_order = 20

    type, public :: analog_prototype_result
        complex(dp), allocatable :: zeros(:) !! finite zeros of H(s), empty for all but cheb2ap
        complex(dp), allocatable :: poles(:) !! left-half-plane poles of H(s)
        real(dp) :: gain = 0.0_dp !! multiplier in gain*product(s-z)/product(s-p)
        integer :: status = prototype_invalid_input !! prototype_success on success
    end type analog_prototype_result

    public :: buttap, cheb1ap, cheb2ap, besselap

contains

    function buttap(n) result(out)
        integer, intent(in) :: n !! prototype order, 0 through 64
        type(analog_prototype_result) :: out
        integer :: i
        real(dp) :: theta

        if (n < 0 .or. n > max_trigonometric_order) return
        allocate(out%zeros(0), out%poles(n))
        do i = 1, n
            theta = pi*real(2*i-n-1, dp)/real(2*n, dp)
            out%poles(i) = -cmplx(cos(theta), sin(theta), dp)
        end do
        out%gain = 1.0_dp
        out%status = prototype_success
    end function buttap

    function cheb1ap(n, rp) result(out)
        integer, intent(in) :: n !! prototype order, 0 through 64
        real(dp), intent(in) :: rp !! positive passband ripple (decibels)
        type(analog_prototype_result) :: out
        real(dp) :: epsilon, mu, theta, log_gain
        integer :: i

        if (n < 0 .or. n > max_trigonometric_order) return
        if (.not. ieee_is_finite(rp)) return
        if (rp <= 0.0_dp .or. rp > 300.0_dp) return
        allocate(out%zeros(0), out%poles(n))
        if (n == 0) then
            out%gain = 10.0_dp**(-rp/20.0_dp)
        else
            epsilon = sqrt(ripple_energy(rp))
            if (epsilon <= 0.0_dp) return
            mu = asinh(1.0_dp/epsilon)/real(n, dp)
            do i = 1, n
                theta = pi*real(2*i-n-1, dp)/real(2*n, dp)
                out%poles(i) = -cmplx(sinh(mu)*cos(theta), cosh(mu)*sin(theta), dp)
            end do
            log_gain = sum(log(abs(out%poles)))
            if (mod(n,2) == 0) log_gain = log_gain - 0.5_dp*log(1.0_dp+epsilon**2)
            out%gain = exp(log_gain)
        end if
        if (.not. ieee_is_finite(out%gain)) return
        out%status = prototype_success
    end function cheb1ap

    function cheb2ap(n, rs) result(out)
        integer, intent(in) :: n !! prototype order, 0 through 64
        real(dp), intent(in) :: rs !! positive stopband attenuation (decibels)
        type(analog_prototype_result) :: out
        real(dp) :: mu, theta, epsilon, log_gain, m
        integer :: i, nz
        complex(dp) :: q

        if (n < 0 .or. n > max_trigonometric_order) return
        if (.not. ieee_is_finite(rs)) return
        if (rs <= 0.0_dp .or. rs > 300.0_dp) return
        allocate(out%zeros(n-mod(n,2)), out%poles(n))
        if (n == 0) then
            out%gain = 1.0_dp
        else
            epsilon = sqrt(ripple_energy(rs))
            if (epsilon <= 0.0_dp) return
            mu = asinh(epsilon)/real(n, dp)
            nz = 0
            do i = 1, n
                m = real(2*i-n-1, dp)
                theta = pi*m/real(2*n, dp)
                q = cmplx(sinh(mu)*cos(theta), cosh(mu)*sin(theta), dp)
                out%poles(i) = -1.0_dp/q
                if (abs(m) > 0.5_dp) then
                    nz = nz+1
                    out%zeros(nz) = cmplx(0.0_dp,1.0_dp/sin(theta),dp)
                end if
            end do
            log_gain = sum(log(abs(out%poles))) - sum(log(abs(out%zeros)))
            out%gain = exp(log_gain)
        end if
        if (.not. ieee_is_finite(out%gain)) return
        out%status = prototype_success
    end function cheb2ap

    ! Evaluate 10**(db/10)-1 accurately for very small ripple levels.
    pure real(dp) function ripple_energy(db) result(energy)
        real(dp), intent(in) :: db !! positive attenuation in decibels
        real(dp) :: x

        x = db*log(10.0_dp)/10.0_dp
        if (x < 1.0e-5_dp) then
            energy = x*(1.0_dp+x*(0.5_dp+x*(1.0_dp/6.0_dp+x/24.0_dp)))
        else
            energy = exp(x)-1.0_dp
        end if
    end function ripple_energy

    function besselap(n, norm) result(out)
        integer, intent(in) :: n !! order, 0 to 20 for the DP root solver
        character(*), intent(in), optional :: norm !! 'phase' (default), 'delay', or 'mag'
        type(analog_prototype_result) :: out
        character(:), allocatable :: choice
        real(dp) :: log_last, frequency_scale, log_k
        integer :: i
        logical :: ok

        choice = 'phase'
        if (present(norm)) choice = trim(norm)
        if (choice /= 'phase' .and. choice /= 'delay' .and. choice /= 'mag') return
        if (n < 0) return
        if (n > max_bessel_order) then
            out%status = prototype_too_expensive
            return
        end if
        allocate(out%zeros(0), out%poles(n))
        if (n == 0) then
            out%gain = 1.0_dp
            out%status = prototype_success
            return
        end if
        call reverse_bessel_roots(n, out%poles, ok)
        if (.not. ok) then
            out%status = prototype_no_convergence
            return
        end if
        log_last = 0.0_dp
        do i = 1, n
            log_last = log_last + log(real(n+i,dp)/2.0_dp)
        end do
        select case (choice)
        case ('phase')
            frequency_scale = exp(-log_last/real(n,dp))
            out%poles = out%poles*frequency_scale
            out%gain = 1.0_dp
        case ('delay')
            out%gain = exp(log_last)
        case ('mag')
            frequency_scale = bessel_half_power_frequency(out%poles)
            if (.not. ieee_is_finite(frequency_scale) .or. frequency_scale <= 0.0_dp) then
                out%status = prototype_no_convergence
                return
            end if
            out%poles = out%poles/frequency_scale
            log_k = log_last-real(n,dp)*log(frequency_scale)
            out%gain = exp(log_k)
        end select
        if (.not. ieee_is_finite(out%gain)) return
        out%status = prototype_success
    end function besselap

    subroutine reverse_bessel_roots(n, poles, ok)
        integer, intent(in) :: n !! positive polynomial degree, at most 20
        complex(dp), intent(out) :: poles(:) !! poles in delay normalization
        logical, intent(out) :: ok !! true for distinct, stable, validated poles
        real(dp) :: a1, a2, b0, b1, b2, b3, s, r, x, k
        real(dp) :: last_step, threshold
        complex(dp) :: roots(n), next_roots(n), f, deriv, delta, repulsion, step, tmp
        integer :: i, j, iteration, pair

        ok = .false.
        ! The ordinary Bessel polynomial y_n(x) has zeros x = 1/p,
        ! where p are the poles of its reverse polynomial theta_n(p).
        ! Seed ordinary zeros using Campos-Calderon's polynomial approximations.
        ! This avoids loss of significance in reverse-polynomial coefficients.
        if (n == 1) then
            roots(1) = cmplx(-1.0_dp,0.0_dp,dp)
        else
            s = 2.0_dp*n**2-3.0_dp*n**4+real(n,dp)**5
            r = 2.0_dp*n**2+real(n,dp)**3
            a1 = (-6.0_dp-6.0_dp*n)/r
            a2 = 6.0_dp/r
            b3 = (16.0_dp-8.0_dp*n)/s
            b2 = (-24.0_dp-12.0_dp*n+12.0_dp*n**2)/s
            b1 = (8.0_dp+24.0_dp*n-12.0_dp*n**2-2.0_dp*n**3)/s
            b0 = (-6.0_dp*n+5.0_dp*n**3-real(n,dp)**4)/s
            do i = 1, n
                k = real(i,dp)
                x = a1*k+a2*k*k
                roots(i) = cmplx(x,b0+b1*k+b2*k*k+b3*k**3,dp)
            end do
        end if

        ! The zeros of y_n(x) coincide with those of exp(1/x)*K_(n+1/2)(1/x).
        ! SciFort's scaled complex K implementation avoids catastrophic
        ! cancellation in direct polynomial Horner evaluation at higher orders.
        do iteration = 1, 100
            last_step = 0.0_dp
            do i = 1, n
                call bessel_zero_function(n, roots(i), f, deriv)
                if (.not. finite_complex_scalar(f) .or. &
                    .not. finite_complex_scalar(deriv)) return
                if (abs(deriv) <= tiny(1.0_dp)) return
                delta = f/deriv
                repulsion = (0.0_dp,0.0_dp)
                do j = 1, n
                    if (i == j) cycle
                    if (abs(roots(i)-roots(j)) <= 1.0e-14_dp) return
                    repulsion = repulsion+1.0_dp/(roots(i)-roots(j))
                end do
                if (abs(1.0_dp-delta*repulsion) <= tiny(1.0_dp)) return
                step = delta/(1.0_dp-delta*repulsion)
                next_roots(i) = roots(i)-step
                last_step = max(last_step,abs(step))
            end do
            roots = next_roots
            if (.not. all(ieee_is_finite(real(roots,dp)))) return
            if (.not. all(ieee_is_finite(aimag(roots)))) return
            threshold = 3.0e-12_dp*max(1.0_dp,maxval(abs(roots)))
            if (last_step <= threshold) exit
        end do
        if (iteration > 100) return
        if (any(abs(roots) <= tiny(1.0_dp))) return
        poles = 1.0_dp/roots

        ! Sort by imaginary part, then enforce exactly conjugate poles.
        do i = 2, n
            tmp = poles(i)
            j = i-1
            do while (j >= 1)
                if (aimag(poles(j)) <= aimag(tmp)) exit
                poles(j+1) = poles(j)
                j = j-1
            end do
            poles(j+1) = tmp
        end do
        do i = 1, n/2
            pair = n+1-i
            tmp = 0.5_dp*(poles(i)+conjg(poles(pair)))
            poles(i) = tmp
            poles(pair) = conjg(tmp)
        end do
        if (mod(n,2) == 1) poles((n+1)/2) = cmplx(real(poles((n+1)/2),dp),0.0_dp,dp)
        if (any(real(poles,dp) >= 0.0_dp)) return
        ! Vieta's sum catches lost roots and spurious converged iterations.
        if (abs(sum(real(poles,dp))+0.5_dp*n*(n+1)) > 1.0e-8_dp*real(n,dp)) return
        do i = 1, n
            do j = 1, i-1
                if (abs(poles(i)-poles(j)) < 1.0e-8_dp) return
            end do
        end do
        ok = .true.
    end subroutine reverse_bessel_roots

    pure subroutine bessel_zero_function(n, x, f, df)
        integer, intent(in) :: n !! degree of the ordinary Bessel polynomial
        complex(dp), intent(in) :: x !! complex nonzero evaluation point
        complex(dp), intent(out) :: f !! scaled modified-Bessel-K value
        complex(dp), intent(out) :: df !! derivative of the scaled value with respect to x
        complex(dp) :: z, km, kp

        z = 1.0_dp/x
        f = kve_complex(real(n,dp)+0.5_dp,z)
        km = kve_complex(real(n,dp)-0.5_dp,z)
        kp = kve_complex(real(n,dp)+1.5_dp,z)
        ! K'_nu(z) = -(K_(nu-1)(z)+K_(nu+1)(z))/2, and z=1/x.
        df = (0.5_dp*(km+kp)-f)/(x*x)
    end subroutine bessel_zero_function

    pure logical function finite_complex_scalar(z) result(ok)
        complex(dp), intent(in) :: z !! value checked for finite parts
        ok = ieee_is_finite(real(z,dp)) .and. ieee_is_finite(aimag(z))
    end function finite_complex_scalar

    pure real(dp) function bessel_half_power_frequency(poles) result(w)
        complex(dp), intent(in) :: poles(:) !! delay-normalized poles
        real(dp) :: lower, upper, midpoint, log_magnitude
        integer :: i

        lower = 0.0_dp
        upper = real(size(poles)*4,dp)
        do i = 1, 90
            midpoint = 0.5_dp*(lower+upper)
            log_magnitude = sum(log(abs(poles)) - &
                log(abs(cmplx(0.0_dp,midpoint,dp)-poles)))
            if (log_magnitude > -0.5_dp*log(2.0_dp)) then
                lower = midpoint
            else
                upper = midpoint
            end if
        end do
        w = 0.5_dp*(lower+upper)
    end function bessel_half_power_frequency

end module scifort_signal_analog_prototypes
