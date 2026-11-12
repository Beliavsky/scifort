! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Frequency-sampling FIR design, weighted least-squares FIR design, and
! FFT-cepstrum minimum-phase conversion. Independent implementations of
! standard filter-design formulas. See CODE_PROVENANCE.md.
module scifort_signal_fir_advanced
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_signal_spectral, only : forward_fft, spectral_named_window, spectral_status_success
    implicit none
    private

    integer, parameter, public :: fir_advanced_status_success = 0
    integer, parameter, public :: fir_advanced_status_invalid_input = 1
    integer, parameter, public :: fir_advanced_status_singular = 2
    real(dp), parameter :: pi = acos(-1.0_dp)
    public :: firwin2, firls, minimum_phase

contains

    ! Forward FFT is exp(-i*omega*t); its inverse is the conjugated transform.
    subroutine inverse_fft(values)
        complex(dp), intent(inout) :: values(:) !! complex samples transformed in place
        values = conjg(values)
        call forward_fft(values)
        values = conjg(values) / real(size(values), dp)
    end subroutine inverse_fft

    integer function next_power_two(n) result(out)
        integer, intent(in) :: n !! requested minimum transform length
        out = 1
        if (n < 1) then
            out = 0
            return
        end if
        do while (out < n)
            if (out > huge(out)/2) then
                out = 0
                return
            end if
            out = 2*out
        end do
    end function next_power_two

    function firwin2(numtaps, freq, gain, nfreqs, window, antisymmetric, fs, status, &
            window_parameter, window_parameter2, window_coefficients, window_values) result(taps)
        integer, intent(in) :: numtaps !! FIR coefficient count, positive and smaller than nfreqs
        real(dp), intent(in) :: freq(:) !! nondecreasing frequencies, including 0 and fs/2
        real(dp), intent(in) :: gain(:) !! gain at each specified frequency
        integer, intent(in), optional :: nfreqs !! one-sided frequency grid size; default next power of two plus one
        character(*), intent(in), optional :: window !! symmetric window, default hamming; 'none' for no window
        logical, intent(in), optional :: antisymmetric !! generate type III or type IV FIR
        real(dp), intent(in), optional :: fs !! sample rate, default 2
        integer, intent(out), optional :: status !! success or invalid input
        real(dp), intent(in), optional :: window_parameter !! first named-window shape parameter
        real(dp), intent(in), optional :: window_parameter2 !! second shape parameter
        real(dp), intent(in), optional :: window_coefficients(:) !! general-cosine window weights
        real(dp), intent(in), optional :: window_values(:) !! explicit symmetric window, size numtaps
        real(dp), allocatable :: taps(:)
        real(dp), allocatable :: frequencies(:), win(:)
        complex(dp), allocatable :: spectrum(:)
        real(dp) :: rate, nyquist, eps_freq, x, y, phase
        integer :: npoints, ngrid, nfft, i, j, seg, st
        logical :: anti
        character(:), allocatable :: wname

        allocate(taps(0))
        if (present(status)) status = fir_advanced_status_invalid_input
        rate = 2.0_dp
        if (present(fs)) rate = fs
        if (.not. ieee_is_finite(rate)) return
        if (rate <= 0.0_dp .or. numtaps < 1) return
        nyquist = rate/2.0_dp
        npoints = size(freq)
        if (npoints < 2 .or. size(gain) /= npoints) return
        if (.not. all(ieee_is_finite(freq)) .or. .not. all(ieee_is_finite(gain))) return
        if (freq(1) /= 0.0_dp .or. freq(npoints) /= nyquist) return
        if (any(freq(2:) < freq(:npoints-1))) return
        if (freq(2) == 0.0_dp .or. freq(npoints-1) == nyquist) return
        anti = .false.
        if (present(antisymmetric)) anti = antisymmetric
        if (anti) then
            if (gain(1) /= 0.0_dp) return
            if (mod(numtaps, 2) == 1 .and. gain(npoints) /= 0.0_dp) return
        else
            if (mod(numtaps, 2) == 0 .and. gain(npoints) /= 0.0_dp) return
        end if

        ngrid = next_power_two(numtaps)
        if (ngrid == 0 .or. ngrid == huge(ngrid)) return
        ngrid = ngrid + 1
        if (present(nfreqs)) ngrid = nfreqs
        if (ngrid <= numtaps .or. ngrid < 2 .or. ngrid-1 > huge(ngrid)/2) return
        nfft = 2*(ngrid-1)

        frequencies = freq
        eps_freq = epsilon(1.0_dp)*nyquist
        do i = 1, npoints-1
            if (freq(i) == freq(i+1)) then
                frequencies(i) = freq(i)-eps_freq
                frequencies(i+1) = freq(i+1)+eps_freq
            end if
        end do
        if (any(frequencies(2:) <= frequencies(:npoints-1))) return
        wname = 'hamming'
        if (present(window)) wname = trim(window)
        if (present(window_values)) then
            ! Unlike SciPy get_window, SciFort also permits a precomputed
            ! symmetric window, provided that no named window is requested.
            if (present(window) .or. present(window_parameter) .or. &
                    present(window_parameter2) .or. present(window_coefficients)) return
            if (size(window_values) /= numtaps) return
            if (.not. all(ieee_is_finite(window_values))) return
            win = window_values
        else if (wname == 'none') then
            if (present(window_parameter) .or. present(window_parameter2) .or. &
                    present(window_coefficients)) return
            allocate(win(numtaps))
            win = 1.0_dp
        else
            win = spectral_named_window(wname, numtaps, st, &
                parameter=window_parameter, parameter2=window_parameter2, &
                coefficients=window_coefficients, fftbins=.false.)
            if (st /= spectral_status_success) return
        end if

        allocate(spectrum(nfft))
        spectrum = (0.0_dp, 0.0_dp)
        seg = 1
        do j = 0, ngrid-1
            x = nyquist * real(j,dp) / real(ngrid-1,dp)
            do while (seg < npoints-1)
                if (x <= frequencies(seg+1)) exit
                seg = seg + 1
            end do
            y = gain(seg) + (gain(seg+1)-gain(seg)) * &
                (x-frequencies(seg))/(frequencies(seg+1)-frequencies(seg))
            phase = -0.5_dp * real(numtaps-1,dp)*pi*x/nyquist
            spectrum(j+1) = y*cmplx(cos(phase), sin(phase), dp)
            if (anti) spectrum(j+1) = cmplx(0.0_dp,1.0_dp,dp)*spectrum(j+1)
        end do
        ! The DC and Nyquist bins of a real inverse transform must be real.
        spectrum(1) = cmplx(real(spectrum(1),dp),0.0_dp,dp)
        spectrum(ngrid) = cmplx(real(spectrum(ngrid),dp),0.0_dp,dp)
        do j = 1, ngrid-2
            spectrum(nfft-j+1) = conjg(spectrum(j+1))
        end do
        call inverse_fft(spectrum)
        taps = real(spectrum(1:numtaps),dp)*win
        if (anti .and. mod(numtaps,2) == 1) taps((numtaps+1)/2) = 0.0_dp
        if (.not. all(ieee_is_finite(taps))) then
            deallocate(taps)
            allocate(taps(0))
            return
        end if
        if (present(status)) status = fir_advanced_status_success
    end function firwin2

    logical function cholesky_solve(a, b, x) result(ok)
        real(dp), intent(in) :: a(:,:) !! symmetric positive definite system matrix
        real(dp), intent(in) :: b(:) !! system right-hand side
        real(dp), intent(out) :: x(:) !! solution vector
        real(dp) :: ell(size(b),size(b)), y(size(b)), t, floor_diag
        integer :: n, i, j, k

        ok = .false.
        n = size(b)
        if (size(a,1) /= n .or. size(a,2) /= n .or. size(x) /= n) return
        ell = 0.0_dp
        floor_diag = epsilon(1.0_dp)*maxval(abs(a))
        do i = 1, n
            do j = 1, i
                t = a(i,j)
                do k = 1, j-1
                    t = t - ell(i,k)*ell(j,k)
                end do
                if (i == j) then
                    if (.not. ieee_is_finite(t)) return
                    if (t <= floor_diag) return
                    ell(i,j) = sqrt(t)
                else
                    ell(i,j) = t/ell(j,j)
                end if
            end do
        end do
        do i = 1, n
            t = b(i)
            do k = 1, i-1
                t = t-ell(i,k)*y(k)
            end do
            y(i) = t/ell(i,i)
        end do
        do i = n, 1, -1
            t = y(i)
            do k = i+1, n
                t = t-ell(k,i)*x(k)
            end do
            x(i) = t/ell(i,i)
        end do
        ok = all(ieee_is_finite(x))
    end function cholesky_solve

    pure real(dp) function sinc_small(x) result(y)
        real(dp), intent(in) :: x !! real-valued argument
        real(dp) :: xx
        if (abs(x) < 0.01_dp) then
            xx = x*x
            y = 1.0_dp + xx*(-1.0_dp/6.0_dp + xx*(1.0_dp/120.0_dp - xx/5040.0_dp))
        else
            y = sin(x)/x
        end if
    end function sinc_small

    pure real(dp) function sine_moment(x) result(y)
        ! (sin(x)-x*cos(x))/x**2; avoid cancellation for small x.
        real(dp), intent(in) :: x !! real-valued argument
        real(dp) :: xx
        if (abs(x) < 0.01_dp) then
            xx = x*x
            y = x*(1.0_dp/3.0_dp + xx*(-1.0_dp/30.0_dp + &
                xx*(1.0_dp/840.0_dp - xx/45360.0_dp)))
        else
            y = (sin(x)-x*cos(x))/(x*x)
        end if
    end function sine_moment

    logical function symmetric_rank_solve(a, b, x) result(ok)
        ! Rank-revealing symmetric Jacobi eigensolver. Use the pseudoinverse
        ! of the normal-equation matrix only when fast Cholesky fails.
        ! This is an independent fallback, requiring no BLAS/LAPACK library.
        real(dp), intent(in) :: a(:,:) !! symmetric positive semidefinite matrix
        real(dp), intent(in) :: b(:) !! right-hand side
        real(dp), intent(out) :: x(:) !! minimum-norm solution to retained rank
        real(dp) :: mat(size(b),size(b)), v(size(b),size(b))
        real(dp) :: scale, off, diagmax, eig_tol, tau, t, c, sn
        real(dp) :: app, aqq, apq, aip, aiq, vip, viq, lam, projection
        integer :: n, p, q, i, sweep, rank
        logical :: converged

        ok = .false.
        n = size(b)
        x = 0.0_dp
        if (n < 1 .or. size(a,1) /= n .or. size(a,2) /= n .or. size(x) /= n) return
        if (.not. all(ieee_is_finite(a)) .or. .not. all(ieee_is_finite(b))) return
        scale = maxval(abs(a))
        if (scale <= 0.0_dp) return
        mat = a/scale
        v = 0.0_dp
        do i = 1, n
            v(i,i) = 1.0_dp
        end do
        ! Cyclic Jacobi diagonalization of the symmetric Gram matrix.
        ! It resolves eigenvalues down to the precision of the matrix itself.
        converged = .false.
        do sweep = 1, 60
            diagmax = maxval(abs([(mat(i,i), i=1,n)]))
            off = 0.0_dp
            do p = 1, n-1
                do q = p+1, n
                    off = max(off,abs(mat(p,q)))
                end do
            end do
            if (off <= 4.0_dp*epsilon(1.0_dp)*diagmax) then
                converged = .true.
                exit
            end if
            do p = 1, n-1
                do q = p+1, n
                    apq = mat(p,q)
                    if (abs(apq) <= 2.0_dp*epsilon(1.0_dp)*diagmax) cycle
                    app = mat(p,p)
                    aqq = mat(q,q)
                    tau = (aqq-app)/(2.0_dp*apq)
                    if (tau == 0.0_dp) then
                        t = 1.0_dp
                    else
                        t = sign(1.0_dp,tau)/(abs(tau)+sqrt(1.0_dp+tau*tau))
                    end if
                    c = 1.0_dp/sqrt(1.0_dp+t*t)
                    sn = t*c
                    mat(p,p) = app-t*apq
                    mat(q,q) = aqq+t*apq
                    mat(p,q) = 0.0_dp
                    mat(q,p) = 0.0_dp
                    do i = 1, n
                        if (i /= p .and. i /= q) then
                            aip = mat(i,p)
                            aiq = mat(i,q)
                            mat(i,p) = c*aip-sn*aiq
                            mat(p,i) = mat(i,p)
                            mat(i,q) = sn*aip+c*aiq
                            mat(q,i) = mat(i,q)
                        end if
                        vip = v(i,p)
                        viq = v(i,q)
                        v(i,p) = c*vip-sn*viq
                        v(i,q) = sn*vip+c*viq
                    end do
                end do
            end do
        end do
        if (.not. converged) return
        eig_tol = real(n,dp)*epsilon(1.0_dp)*maxval(abs([(mat(i,i),i=1,n)]))
        rank = 0
        do i = 1, n
            lam = mat(i,i)
            if (lam <= eig_tol) cycle
            rank = rank+1
            projection = dot_product(v(:,i), b/scale)/lam
            x = x+projection*v(:,i)
        end do
        ok = rank > 0 .and. all(ieee_is_finite(x))
    end function symmetric_rank_solve

    function firls(numtaps, bands, desired, weight, fs, status) result(taps)
        integer, intent(in) :: numtaps !! positive odd filter length
        real(dp), intent(in) :: bands(:) !! frequency edge pairs in [0, fs/2], monotone nonoverlapping
        real(dp), intent(in) :: desired(:) !! gain at each band edge, nonnegative
        real(dp), intent(in), optional :: weight(:) !! nonnegative per-band weights, default 1
        real(dp), intent(in), optional :: fs !! sample rate, default 2
        integer, intent(out), optional :: status !! success, invalid input, or singular
        real(dp), allocatable :: taps(:)
        real(dp), allocatable :: q(:), rhs(:), a(:,:), sol(:), w(:)
        real(dp) :: rate, nyquist, lo, hi, width, mid, dmid, dd, z, x, base
        integer :: nb, m, i, k, j
        logical :: ok

        allocate(taps(0))
        if (present(status)) status = fir_advanced_status_invalid_input
        rate = 2.0_dp
        if (present(fs)) rate = fs
        if (.not. ieee_is_finite(rate)) return
        if (rate <= 0.0_dp .or. numtaps < 1 .or. mod(numtaps,2) == 0) return
        if (size(bands) < 2 .or. mod(size(bands),2) /= 0) return
        if (size(desired) /= size(bands)) return
        if (.not. all(ieee_is_finite(bands)) .or. .not. all(ieee_is_finite(desired))) return
        if (any(desired < 0.0_dp)) return
        nyquist = rate/2.0_dp
        if (any(bands < 0.0_dp) .or. any(bands > nyquist)) return
        nb = size(bands)/2
        do i = 1, nb
            if (bands(2*i) <= bands(2*i-1)) return
            if (i > 1) then
                if (bands(2*i-1) < bands(2*i-2)) return
            end if
        end do
        allocate(w(nb))
        w = 1.0_dp
        if (present(weight)) then
            if (size(weight) /= nb) return
            if (.not. all(ieee_is_finite(weight))) return
            if (any(weight < 0.0_dp)) return
            w = weight
        end if
        if (all(w == 0.0_dp)) then
            if (present(status)) status = fir_advanced_status_singular
            return
        end if
        ! Scaling all band weights does not change the minimizer. This also
        ! protects the Gram matrix against overflow for very large weights.
        w = w/maxval(w)

        m = (numtaps-1)/2
        allocate(q(numtaps),rhs(m+1),a(m+1,m+1),sol(m+1))
        q = 0.0_dp
        rhs = 0.0_dp
        do i = 1, nb
            if (w(i) == 0.0_dp) cycle
            lo = bands(2*i-1)/nyquist
            hi = bands(2*i)/nyquist
            width = hi-lo
            mid = lo + 0.5_dp*width
            dmid = desired(2*i-1)+0.5_dp*(desired(2*i)-desired(2*i-1))
            dd = desired(2*i)-desired(2*i-1)
            q(1) = q(1) + w(i)*width
            rhs(1) = rhs(1) + w(i)*width*dmid
            do k = 1, numtaps-1
                z = pi*real(k,dp)
                x = 0.5_dp*z*width
                base = width*cos(z*mid)*sinc_small(x)
                q(k+1) = q(k+1) + w(i)*base
                if (k <= m) then
                    rhs(k+1) = rhs(k+1) + w(i)*width*( &
                        dmid*cos(z*mid)*sinc_small(x) - &
                        0.5_dp*dd*sin(z*mid)*sine_moment(x))
                end if
            end do
        end do
        if (.not. all(ieee_is_finite(q)) .or. .not. all(ieee_is_finite(rhs))) return
        do j = 0, m
            do k = 0, m
                a(j+1,k+1) = q(abs(j-k)+1)+q(j+k+1)
            end do
        end do
        ok = cholesky_solve(a,rhs,sol)
        if (.not. ok) ok = symmetric_rank_solve(a,rhs,sol)
        if (.not. ok) then
            if (present(status)) status = fir_advanced_status_singular
            return
        end if
        deallocate(taps)
        allocate(taps(numtaps))
        taps(m+1) = 2.0_dp*sol(1)
        do j = 1, m
            taps(m+1-j) = sol(j+1)
            taps(m+1+j) = sol(j+1)
        end do
        if (.not. all(ieee_is_finite(taps))) then
            deallocate(taps)
            allocate(taps(0))
            if (present(status)) status = fir_advanced_status_singular
            return
        end if
        if (present(status)) status = fir_advanced_status_success
    end function firls

    function minimum_phase(h, method, n_fft, half, status) result(taps)
        real(dp), intent(in) :: h(:) !! real-valued FIR coefficients, at least 3
        character(*), intent(in), optional :: method !! 'homomorphic' (default) or 'hilbert'
        integer, intent(in), optional :: n_fft !! FFT length, default power of 2 >= 200*(size(h)-1)
        logical, intent(in), optional :: half !! half-length square-root magnitude, default true
        integer, intent(out), optional :: status !! success, invalid input, or singular
        real(dp), allocatable :: taps(:)
        complex(dp), allocatable :: spectrum(:), temp(:)
        real(dp), allocatable :: amplitude(:)
        real(dp) :: a_min, delta_p, delta_s, gain_scale, theta
        integer :: n, nf, i, nout, mid
        logical :: half_output
        character(:), allocatable :: mode

        allocate(taps(0))
        if (present(status)) status = fir_advanced_status_invalid_input
        n = size(h)
        if (n < 3 .or. .not. all(ieee_is_finite(h))) return
        half_output = .true.
        if (present(half)) half_output = half
        mode = 'homomorphic'
        if (present(method)) mode = trim(method)
        if (mode /= 'homomorphic' .and. mode /= 'hilbert') return
        if (mode == 'hilbert' .and. .not. half_output) return
        if (n-1 > huge(n)/200) return
        nf = next_power_two(200*(n-1))
        if (present(n_fft)) nf = n_fft
        if (nf < n) return
        nout = n
        if (half_output) nout = (n+1)/2
        allocate(spectrum(nf),temp(nf),amplitude(nf))
        spectrum = (0.0_dp,0.0_dp)
        spectrum(1:n) = cmplx(h,0.0_dp,dp)
        call forward_fft(spectrum)

        if (mode == 'homomorphic') then
            amplitude = abs(spectrum)
            if (.not. any(amplitude > 0.0_dp)) then
                if (present(status)) status = fir_advanced_status_singular
                return
            end if
            a_min = minval(amplitude, mask=amplitude>0.0_dp)
            amplitude = log(amplitude + 1.0e-7_dp*a_min)
            if (half_output) amplitude = 0.5_dp*amplitude
            temp = cmplx(amplitude,0.0_dp,dp)
            call inverse_fft(temp)
            ! The one-sided complex cepstrum has twice the strictly
            ! positive quefrencies, with DC and (if present) Nyquist single.
            temp = cmplx(real(temp,dp),0.0_dp,dp)
            mid = nf/2
            do i = 2, mid
                temp(i) = 2.0_dp*temp(i)
            end do
            if (mod(nf,2) == 1) temp(mid+1) = cmplx(real(temp(mid+1),dp),0.0_dp,dp)
            if (mod(nf,2) == 0) temp(mid+1) = (0.0_dp,0.0_dp)
            if (mid+2 <= nf) temp(mid+2:nf) = (0.0_dp,0.0_dp)
            call forward_fft(temp)
            temp = exp(temp)
            call inverse_fft(temp)
        else
            mid = n/2
            do i = 0, nf-1
                theta = 2.0_dp*pi*real(mid,dp)*real(i,dp)/real(nf,dp)
                amplitude(i+1) = real(spectrum(i+1)*cmplx(cos(theta),sin(theta),dp),dp)
            end do
            delta_p = maxval(amplitude)-1.0_dp
            delta_s = -minval(amplitude)
            if (1.0_dp+delta_p+delta_s < 0.0_dp .or. &
                1.0_dp-delta_p+delta_s < 0.0_dp) then
                if (present(status)) status = fir_advanced_status_singular
                return
            end if
            gain_scale = sqrt(1.0_dp+delta_p+delta_s)+sqrt(1.0_dp-delta_p+delta_s)
            if (gain_scale <= 0.0_dp) then
                if (present(status)) status = fir_advanced_status_singular
                return
            end if
            amplitude = max(0.0_dp,(amplitude+delta_s)*(4.0_dp/gain_scale**2))
            amplitude = sqrt(amplitude)+1.0e-10_dp
            temp = cmplx(log(amplitude),0.0_dp,dp)
            call inverse_fft(temp)
            do i = 1, nf
                if (i == 1 .or. i == nf/2+1) then
                    temp(i) = (0.0_dp,0.0_dp)
                else if (i >= 2 .and. i <= nf/2) then
                    temp(i) = temp(i)
                else
                    temp(i) = -temp(i)
                end if
            end do
            call forward_fft(temp)
            temp = cmplx(amplitude,0.0_dp,dp)*exp(temp)
            call inverse_fft(temp)
        end if
        if (.not. all(ieee_is_finite(real(temp(:nout),dp)))) then
            if (present(status)) status = fir_advanced_status_singular
            return
        end if
        taps = real(temp(:nout),dp)
        if (present(status)) status = fir_advanced_status_success
    end function minimum_phase
end module scifort_signal_fir_advanced
