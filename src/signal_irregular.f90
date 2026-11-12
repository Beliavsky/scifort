! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Independent weighted least-squares Lomb-Scargle spectral analysis and
! event phase synchrony, following the published model specifications.
module scifort_signal_irregular
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    implicit none
    private

    integer, parameter, public :: irregular_status_success = 0
    integer, parameter, public :: irregular_status_invalid_input = 1
    integer, parameter, public :: irregular_status_singular = 2
    real(dp), parameter :: twopi = 2.0_dp * acos(-1.0_dp)

    type, public :: vector_strength_result
        real(dp), allocatable :: strength(:) !! length equal to number of periods
        real(dp), allocatable :: phase(:) !! mean phase, in radians (-pi, pi]
        integer :: status = irregular_status_invalid_input !! zero on success
    end type vector_strength_result

    ! At each angular frequency, return the actual weighted Fourier-series
    ! coefficients as well as the periodogram. Frequency bins can fail
    ! independently if the harmonic design is rank deficient.
    type, public :: harmonic_fit_result
        real(dp), allocatable :: frequencies(:) !! angular trial frequencies
        real(dp), allocatable :: offset(:) !! fitted intercept on the original observation scale
        real(dp), allocatable :: cosine(:,:) !! cosine coefficients (nterms,nfreq)
        real(dp), allocatable :: sine(:,:) !! sine coefficients (nterms,nfreq)
        real(dp), allocatable :: power(:) !! one-half N times weighted explained sum of squares
        real(dp), allocatable :: normalized_power(:) !! fractional weighted variance explained
        real(dp), allocatable :: rms(:) !! square root of weighted mean squared residual
        logical, allocatable :: valid(:) !! .true. at full-rank trial frequencies
        ! Only allocated when estimate_covariance=.true.  Order is
        ! (intercept, cos(1), sin(1), ..., cos(nterms), sin(nterms)).
        real(dp), allocatable :: covariance(:,:,:) !! parameter covariance (nparam+1,nparam+1,nfreq)
        real(dp), allocatable :: standard_error(:,:) !! parameter SE (nparam+1,nfreq)
        real(dp), allocatable :: residual_variance(:) !! estimated weighted noise scale per frequency
        integer :: degrees_of_freedom = 0 !! number of positive weights minus fitted parameter count
        logical :: intercept_fitted = .false. !! true if floating_mean fitted an intercept
        integer :: nterms = 0 !! number of fitted harmonics
        integer :: status = irregular_status_invalid_input !! success, invalid input, or singular fit
    end type harmonic_fit_result

    ! Reuse the frequency-dependent weighted harmonic QR factorization for
    ! many observation series sampled on the same irregular time grid.
    ! The cached basis is private: callers can read status and valid bins but
    ! cannot modify the factorization or its normalized weights.
    type, public :: harmonic_fit_plan
        integer :: status = irregular_status_invalid_input !! validity of prepared frequency grid
        integer :: nterms = 0 !! number of harmonics in each trial fit
        logical :: floating_mean = .false. !! whether to fit an intercept
        logical, allocatable :: valid(:) !! full-rank frequency bins
        real(dp), allocatable :: frequencies(:) !! angular frequencies (radians/time)
        integer, private :: nobs = 0 !! number of time samples
        integer, private :: dof = 0 !! positive-weight observations minus parameter count
        real(dp), allocatable, private :: weights(:) !! sum-one normalized observation weights
        real(dp), allocatable, private :: root_weights(:) !! square roots of normalized weights
        real(dp), allocatable, private :: q(:,:,:) !! orthonormal weighted design columns (n,p,m)
        real(dp), allocatable, private :: r(:,:,:) !! upper triangular factors (p,p,m)
        real(dp), allocatable, private :: column_means(:,:) !! harmonic weighted means (p,m)
    end type harmonic_fit_plan

    interface vectorstrength
        module procedure vectorstrength_array
        module procedure vectorstrength_scalar
    end interface vectorstrength

    public :: lombscargle, lombscargle_amplitude, lombscargle_multiharmonic, vectorstrength
    public :: lombscargle_harmonic_fit, lombscargle_harmonic_predict
    public :: lombscargle_harmonic_predict_std
    public :: lombscargle_harmonic_prepare, lombscargle_harmonic_fit_prepared
    public :: lombscargle_harmonic_fit_batch

contains

    function lombscargle(x, y, freqs, precenter, normalize, weights, floating_mean, status, nterms) result(pgram)
        real(dp), intent(in) :: x(:) !! sample locations, usually observation times
        real(dp), intent(in) :: y(:) !! real observations at those locations
        real(dp), intent(in) :: freqs(:) !! angular frequencies in radians/time
        logical, intent(in), optional :: precenter !! subtract unweighted arithmetic mean (legacy SciPy mode)
        character(*), intent(in), optional :: normalize !! 'power' (default) or 'normalize'
        real(dp), intent(in), optional :: weights(:) !! nonnegative observation weights
        logical, intent(in), optional :: floating_mean !! fit frequency-specific offset; default false
        integer, intent(out), optional :: status !! zero on success
        integer, intent(in), optional :: nterms !! optional Fourier harmonics (SciFort extension), default one
        real(dp), allocatable :: pgram(:)

        real(dp), allocatable :: normalized(:)
        complex(dp), allocatable :: amplitude(:)
        character(:), allocatable :: mode
        integer :: ierr

        if (present(nterms)) then
            if (nterms /= 1) then
                pgram = lombscargle_multiharmonic(x, y, freqs, nterms, precenter, &
                    normalize, weights, floating_mean, status)
                return
            end if
        end if
        mode = 'power'
        if (present(normalize)) mode = trim(normalize)
        allocate(pgram(size(freqs)))
        pgram = ieee_value(0.0_dp, ieee_quiet_nan)
        if (present(status)) status = irregular_status_invalid_input
        if (mode /= 'power' .and. mode /= 'normalize') return

        call lomb_kernel(x, y, freqs, precenter, weights, floating_mean, &
            pgram, normalized, amplitude, ierr)
        if (mode == 'normalize') pgram = normalized
        if (present(status)) status = ierr
    end function lombscargle

    function lombscargle_amplitude(x, y, freqs, precenter, weights, floating_mean, status) result(amplitude)
        real(dp), intent(in) :: x(:) !! sample locations
        real(dp), intent(in) :: y(:) !! observed real signal
        real(dp), intent(in) :: freqs(:) !! angular frequencies in radians/time
        logical, intent(in), optional :: precenter !! subtract unweighted arithmetic mean
        real(dp), intent(in), optional :: weights(:) !! nonnegative observation weights
        logical, intent(in), optional :: floating_mean !! fit frequency-specific offset; default false
        integer, intent(out), optional :: status !! zero on success
        complex(dp), allocatable :: amplitude(:)

        real(dp), allocatable :: power(:), normalized(:)
        integer :: ierr

        call lomb_kernel(x, y, freqs, precenter, weights, floating_mean, &
            power, normalized, amplitude, ierr)
        if (present(status)) status = ierr
    end function lombscargle_amplitude

    subroutine lomb_kernel(x, y, freqs, precenter, weights, floating_mean, power, normalized, amplitude, status)
        real(dp), intent(in) :: x(:) !! observation times
        real(dp), intent(in) :: y(:) !! observations
        real(dp), intent(in) :: freqs(:) !! angular frequencies
        logical, intent(in), optional :: precenter !! subtract arithmetic mean
        real(dp), intent(in), optional :: weights(:) !! nonnegative weights
        logical, intent(in), optional :: floating_mean !! fit frequency-specific offset
        real(dp), allocatable, intent(out) :: power(:) !! unnormalized spectral power
        real(dp), allocatable, intent(out) :: normalized(:) !! fractional variance explained
        complex(dp), allocatable, intent(out) :: amplitude(:) !! complex best-fit amplitude
        integer, intent(out) :: status !! success, invalid input, or singular reference signal

        real(dp), allocatable :: w(:), yv(:)
        real(dp) :: c, s, cc, ss, cs, yc, ys, tau, angle, a, b
        real(dp) :: mean_y, yy, maxweight, improvement, wsum
        real(dp), parameter :: min_denom = 0.5_dp * epsilon(1.0_dp)
        logical :: floating, center
        integer :: i, j, n, m

        n = size(x)
        m = size(freqs)
        allocate(power(m), normalized(m), amplitude(m))
        power = ieee_value(0.0_dp, ieee_quiet_nan)
        normalized = ieee_value(0.0_dp, ieee_quiet_nan)
        amplitude = cmplx(ieee_value(0.0_dp, ieee_quiet_nan), &
            ieee_value(0.0_dp, ieee_quiet_nan), dp)
        status = irregular_status_invalid_input

        if (n == 0 .or. size(y) /= n .or. m == 0) return
        if (any(.not. ieee_is_finite(x)) .or. any(.not. ieee_is_finite(y))) return
        if (any(.not. ieee_is_finite(freqs))) return
        allocate(w(n), yv(n))
        if (present(weights)) then
            if (size(weights) /= n) return
            if (any(.not. ieee_is_finite(weights))) return
            if (any(weights < 0.0_dp)) return
            maxweight = maxval(weights)
            if (maxweight <= 0.0_dp) return
            w = weights / maxweight
        else
            w = 1.0_dp
        end if
        wsum = sum(w)
        if (wsum <= 0.0_dp) return
        w = w / wsum

        yv = y
        center = .false.
        if (present(precenter)) center = precenter
        if (center) yv = yv - sum(y) / real(n, dp)
        mean_y = sum(w * yv)
        floating = .false.
        if (present(floating_mean)) floating = floating_mean
        if (floating) then
            yy = sum(w * (yv - mean_y)**2)
            if (maxval(abs(yv - yv(1))) <= tiny(1.0_dp)) yy = 0.0_dp
        else
            yy = sum(w * yv**2)
        end if

        status = irregular_status_success
        if (yy <= 0.0_dp) status = irregular_status_singular
        do j = 1, m
            ! Rotate the sinusoid basis to eliminate the weighted cosine/sine
            ! cross term. A fitted offset uses centered Gram products.
            c = 0.0_dp
            s = 0.0_dp
            cc = 0.0_dp
            cs = 0.0_dp
            do i = 1, n
                angle = freqs(j) * x(i)
                c = c + w(i) * cos(angle)
                s = s + w(i) * sin(angle)
                cc = cc + w(i) * cos(angle)**2
                cs = cs + w(i) * cos(angle) * sin(angle)
            end do
            ss = 1.0_dp - cc
            if (floating) then
                cc = cc - c*c
                ss = ss - s*s
                cs = cs - c*s
            end if
            tau = 0.5_dp * atan2(2.0_dp*cs, cc-ss)

            cc = 0.0_dp
            c = 0.0_dp
            s = 0.0_dp
            yc = 0.0_dp
            ys = 0.0_dp
            do i = 1, n
                angle = freqs(j) * x(i) - tau
                a = cos(angle)
                b = sin(angle)
                c = c + w(i) * a
                s = s + w(i) * b
                cc = cc + w(i) * a*a
                yc = yc + w(i) * yv(i)*a
                ys = ys + w(i) * yv(i)*b
            end do
            ss = 1.0_dp - cc
            if (floating) then
                yc = yc - mean_y*c
                ys = ys - mean_y*s
                cc = cc - c*c
                ss = ss - s*s
            end if
            ! SciPy uses a half-machine-epsilon floor at singular bins.
            a = yc / max(cc, min_denom)
            b = ys / max(ss, min_denom)
            improvement = a*yc + b*ys
            power(j) = 0.5_dp * real(n, dp) * improvement
            if (yy > 0.0_dp) normalized(j) = improvement / yy
            amplitude(j) = cmplx(a,b,dp) * cmplx(cos(tau), sin(tau),dp)
        end do
    end subroutine lomb_kernel

    ! A SciFort extension of SciPy's generalized Lomb-Scargle periodogram:
    ! jointly fit a Fourier series with nterms harmonics at each trial
    ! frequency.  Weighted, twice-reorthogonalized QR avoids forming the
    ! ill-conditioned normal equations of multi-harmonic regressions.
    function lombscargle_multiharmonic(x, y, freqs, nterms, precenter, normalize, &
            weights, floating_mean, status) result(pgram)
        real(dp), intent(in) :: x(:) !! irregular sample locations (time units)
        real(dp), intent(in) :: y(:) !! real observations, same length as x
        real(dp), intent(in) :: freqs(:) !! angular trial frequencies (radians/time)
        integer, intent(in) :: nterms !! positive number of jointly fitted Fourier harmonics
        logical, intent(in), optional :: precenter !! subtract unweighted sample mean before fitting
        character(*), intent(in), optional :: normalize !! 'power' (default) or 'normalize'
        real(dp), intent(in), optional :: weights(:) !! nonnegative observation weights, same length as x
        logical, intent(in), optional :: floating_mean !! fit a weighted intercept at each frequency
        integer, intent(out), optional :: status !! success, invalid input, or rank-deficient trial frequency
        real(dp), allocatable :: pgram(:)

        real(dp), allocatable :: w(:), rootw(:), centered(:), target(:), q(:,:), v(:)
        real(dp) :: maxweight, mean_y, yy, theta, mean_col, original_norm, norm_v
        real(dp) :: projection, improvement, tol
        character(:), allocatable :: mode
        integer :: n, m, k, col, h, j, pass, st
        logical :: floating, rank_ok

        n = size(x)
        m = size(freqs)
        allocate(pgram(m))
        pgram = ieee_value(0.0_dp, ieee_quiet_nan)
        st = irregular_status_invalid_input
        if (present(status)) status = st
        if (nterms < 1 .or. nterms > (huge(nterms)-1)/2) return
        if (n < 3 .or. size(y) /= n .or. m < 1) return
        if (any(.not. ieee_is_finite(x)) .or. any(.not. ieee_is_finite(y))) return
        if (any(.not. ieee_is_finite(freqs))) return
        mode = 'power'
        if (present(normalize)) mode = trim(normalize)
        if (mode /= 'power' .and. mode /= 'normalize') return

        floating = .false.
        if (present(floating_mean)) floating = floating_mean
        allocate(w(n))
        if (present(weights)) then
            if (size(weights) /= n) return
            if (any(.not. ieee_is_finite(weights))) return
            if (any(weights < 0.0_dp)) return
            maxweight = maxval(weights)
            if (maxweight <= 0.0_dp) return
            w = weights/maxweight
        else
            w = 1.0_dp
        end if
        ! The positive-weight sample count must exceed the number of
        ! fitted parameters; otherwise the residual is not identifiable.
        if (count(w > 0.0_dp) <= 2*nterms + merge(1, 0, floating)) return
        w = w/sum(w)
        rootw = sqrt(w)
        centered = y
        if (present(precenter)) then
            if (precenter) centered = centered - sum(y)/real(n, dp)
        end if
        if (floating) then
            mean_y = sum(w*centered)
            centered = centered - mean_y
        end if
        target = rootw*centered
        yy = dot_product(target, target)
        ! Centering a constant array can leave a roundoff-sized residual.
        if (floating .and. maxval(abs(y-y(1))) == 0.0_dp) yy = 0.0_dp
        if (.not. ieee_is_finite(yy)) return
        if (yy <= 0.0_dp) then
            if (mode == 'power') pgram = 0.0_dp
            if (present(status)) status = irregular_status_singular
            return
        end if

        allocate(q(n, 2*nterms), v(n))
        st = irregular_status_success
        tol = 64.0_dp*epsilon(1.0_dp)*real(2*nterms, dp)
        do k = 1, m
            rank_ok = .true.
            do col = 1, 2*nterms
                h = (col+1)/2
                do j = 1, n
                    theta = real(h, dp)*freqs(k)*x(j)
                    if (.not. ieee_is_finite(theta)) then
                        rank_ok = .false.
                        exit
                    end if
                    if (mod(col, 2) == 1) then
                        v(j) = cos(theta)
                    else
                        v(j) = sin(theta)
                    end if
                end do
                if (.not. rank_ok) exit
                if (floating) then
                    mean_col = sum(w*v)
                    v = v - mean_col
                end if
                v = rootw*v
                original_norm = sqrt(dot_product(v, v))
                ! Twice-reorthogonalized modified Gram-Schmidt (QR).
                do pass = 1, 2
                    do j = 1, col-1
                        projection = dot_product(q(:,j), v)
                        v = v - projection*q(:,j)
                    end do
                end do
                norm_v = sqrt(dot_product(v, v))
                if (norm_v <= tol .or. norm_v <= tol*original_norm) then
                    rank_ok = .false.
                    exit
                end if
                q(:,col) = v/norm_v
            end do
            if (.not. rank_ok) then
                st = irregular_status_singular
                cycle
            end if
            improvement = 0.0_dp
            do col = 1, 2*nterms
                projection = dot_product(q(:,col), target)
                improvement = improvement + projection*projection
            end do
            improvement = min(yy, max(0.0_dp, improvement))
            if (mode == 'power') then
                pgram(k) = 0.5_dp*real(n, dp)*improvement
            else
                pgram(k) = improvement/yy
            end if
        end do
        if (present(status)) status = st
    end function lombscargle_multiharmonic

    ! Recover the actual weighted Fourier-series coefficients from a
    ! multiharmonic Lomb-Scargle fit.  Unlike the periodogram, this solves
    ! the QR triangular system to obtain physical cosine/sine amplitudes.
    ! The intercept is returned on the ORIGINAL y scale, even when
    ! precenter=.true. (the data are then fitted after subtracting their
    ! unweighted mean).
    function lombscargle_harmonic_fit(x, y, freqs, nterms, precenter, weights, &
        floating_mean, estimate_covariance) result(fit)
        real(dp), intent(in) :: x(:) !! irregular observation times
        real(dp), intent(in) :: y(:) !! real observations, same size as x
        real(dp), intent(in) :: freqs(:) !! angular trial frequencies (radians/time)
        integer, intent(in) :: nterms !! number of harmonics, positive
        logical, intent(in), optional :: precenter !! subtract unweighted mean before fitting
        real(dp), intent(in), optional :: weights(:) !! finite nonnegative observation weights
        logical, intent(in), optional :: floating_mean !! jointly estimate frequency-dependent intercept
        logical, intent(in), optional :: estimate_covariance !! estimate regression covariance and standard errors
        type(harmonic_fit_result) :: fit

        real(dp), allocatable :: w(:), rootw(:), centered(:), target(:)
        real(dp), allocatable :: q(:,:), r(:,:), v(:), means(:), z(:), coefficients(:)
        real(dp), allocatable :: inverse_r(:,:), harmonic_covariance(:,:), temp(:)
        real(dp) :: maxweight, shift, mean_y, yy, theta, original_norm, norm_v
        real(dp) :: projection, improvement, tol, intercept, model_value, error_sum
        real(dp) :: residual_scale, value
        integer :: n, m, k, col, h, j, pass, nparam, nfit, i
        logical :: floating, rank_ok, calculate_covariance
        real(dp) :: nan

        n = size(x)
        m = size(freqs)
        nan = ieee_value(0.0_dp, ieee_quiet_nan)
        nfit = 0
        if (nterms >= 1 .and. nterms <= n/2) nfit = nterms
        fit%nterms = nfit
        allocate(fit%frequencies(m))
        fit%frequencies = freqs
        allocate(fit%offset(m), fit%cosine(nfit,m), fit%sine(nfit,m))
        allocate(fit%power(m), fit%normalized_power(m), fit%rms(m), fit%valid(m))
        fit%offset = nan
        fit%cosine = nan
        fit%sine = nan
        fit%power = nan
        fit%normalized_power = nan
        fit%rms = nan
        fit%valid = .false.
        fit%status = irregular_status_invalid_input
        calculate_covariance = .false.
        if (present(estimate_covariance)) calculate_covariance = estimate_covariance
        if (calculate_covariance) then
            allocate(fit%covariance(1+2*nfit,1+2*nfit,m), &
                fit%standard_error(1+2*nfit,m),fit%residual_variance(m))
            fit%covariance = nan
            fit%standard_error = nan
            fit%residual_variance = nan
        end if
        if (nterms < 1 .or. nterms > (huge(nterms)-1)/2) return
        nparam = 2*nterms
        if (n < 3 .or. size(y) /= n .or. m < 1) return
        if (any(.not. ieee_is_finite(x)) .or. any(.not. ieee_is_finite(y))) return
        if (any(.not. ieee_is_finite(freqs))) return
        floating = .false.
        if (present(floating_mean)) floating = floating_mean
        fit%intercept_fitted = floating
        if (nparam + merge(1,0,floating) >= n) return

        allocate(w(n))
        if (present(weights)) then
            if (size(weights) /= n) return
            if (any(.not. ieee_is_finite(weights)) .or. any(weights < 0.0_dp)) return
            maxweight = maxval(weights)
            if (maxweight <= 0.0_dp) return
            w = weights/maxweight
        else
            w = 1.0_dp
        end if
        fit%degrees_of_freedom = count(w > 0.0_dp) - nparam - merge(1,0,floating)
        if (fit%degrees_of_freedom <= 0) return
        w = w/sum(w)
        rootw = sqrt(w)
        centered = y
        shift = 0.0_dp
        if (present(precenter)) then
            if (precenter) then
                shift = sum(y)/real(n,dp)
                centered = centered-shift
            end if
        end if
        mean_y = 0.0_dp
        if (floating) then
            mean_y = sum(w*centered)
            centered = centered-mean_y
        end if
        target = rootw*centered
        yy = dot_product(target,target)
        if (floating .and. maxval(abs(y-y(1))) == 0.0_dp) yy = 0.0_dp
        if (.not. ieee_is_finite(yy)) return
        if (yy <= 0.0_dp) then
            fit%status = irregular_status_singular
            return
        end if

        allocate(q(n,nparam),r(nparam,nparam),v(n),means(nparam),z(nparam),coefficients(nparam))
        if (calculate_covariance) then
            allocate(inverse_r(nparam,nparam),harmonic_covariance(nparam,nparam),temp(nparam))
        end if
        fit%status = irregular_status_success
        tol = 64.0_dp*epsilon(1.0_dp)*real(nparam,dp)
        do k = 1, m
            rank_ok = .true.
            r = 0.0_dp
            means = 0.0_dp
            do col = 1, nparam
                h = (col+1)/2
                do j = 1, n
                    theta = real(h,dp)*freqs(k)*x(j)
                    if (.not. ieee_is_finite(theta)) then
                        rank_ok = .false.
                        exit
                    end if
                    if (mod(col,2) == 1) then
                        v(j) = cos(theta)
                    else
                        v(j) = sin(theta)
                    end if
                end do
                if (.not. rank_ok) exit
                if (floating) then
                    means(col) = sum(w*v)
                    v = v-means(col)
                end if
                v = rootw*v
                original_norm = sqrt(dot_product(v,v))
                do pass = 1, 2
                    do j = 1, col-1
                        projection = dot_product(q(:,j),v)
                        r(j,col) = r(j,col)+projection
                        v = v-projection*q(:,j)
                    end do
                end do
                norm_v = sqrt(dot_product(v,v))
                if (norm_v <= tol .or. norm_v <= tol*original_norm) then
                    rank_ok = .false.
                    exit
                end if
                r(col,col) = norm_v
                q(:,col) = v/norm_v
            end do
            if (.not. rank_ok) then
                fit%status = irregular_status_singular
                cycle
            end if
            do col = 1, nparam
                z(col) = dot_product(q(:,col),target)
            end do
            coefficients = z
            do col = nparam, 1, -1
                if (col < nparam) coefficients(col) = coefficients(col) - &
                    dot_product(r(col,col+1:nparam),coefficients(col+1:nparam))
                coefficients(col) = coefficients(col)/r(col,col)
            end do
            if (any(.not. ieee_is_finite(coefficients))) then
                fit%status = irregular_status_singular
                cycle
            end if
            intercept = shift+mean_y
            if (floating) intercept = intercept-dot_product(means,coefficients)
            if (.not. ieee_is_finite(intercept)) then
                fit%status = irregular_status_singular
                cycle
            end if
            do h = 1, nterms
                fit%cosine(h,k) = coefficients(2*h-1)
                fit%sine(h,k) = coefficients(2*h)
            end do
            fit%offset(k) = intercept
            improvement = min(yy,max(0.0_dp,dot_product(z,z)))
            fit%power(k) = 0.5_dp*real(n,dp)*improvement
            fit%normalized_power(k) = improvement/yy
            error_sum = 0.0_dp
            do j = 1, n
                model_value = intercept
                do h = 1, nterms
                    theta = real(h,dp)*freqs(k)*x(j)
                    model_value = model_value + fit%cosine(h,k)*cos(theta) + &
                        fit%sine(h,k)*sin(theta)
                end do
                error_sum = error_sum+w(j)*(y(j)-model_value)**2
            end do
            if (.not. ieee_is_finite(error_sum)) then
                fit%status = irregular_status_singular
                fit%cosine(:,k) = nan
                fit%sine(:,k) = nan
                fit%offset(k) = nan
                fit%power(k) = nan
                fit%normalized_power(k) = nan
                cycle
            end if
            if (calculate_covariance) then
                ! w is normalized to sum one.  For independent errors whose
                ! variances are inversely proportional to the ORIGINAL
                ! observation weights, the common noise scale is the
                ! weighted RSS divided by the residual degrees of freedom.
                residual_scale = max(0.0_dp,error_sum)/real(fit%degrees_of_freedom,dp)
                if (.not. ieee_is_finite(residual_scale)) then
                    fit%status = irregular_status_singular
                    cycle
                end if
                ! Solve R * inverse_r = I by backward substitution.
                ! Unlike a normal-equations inverse, this uses the same
                ! stable QR factor that supplied the coefficients.
                inverse_r = 0.0_dp
                do col = 1, nparam
                    do i = nparam, 1, -1
                        value = 0.0_dp
                        if (i == col) value = 1.0_dp
                        if (i < nparam) value = value - &
                            dot_product(r(i,i+1:nparam),inverse_r(i+1:nparam,col))
                        inverse_r(i,col) = value/r(i,i)
                    end do
                end do
                harmonic_covariance = residual_scale * matmul(inverse_r,transpose(inverse_r))
                if (any(.not. ieee_is_finite(harmonic_covariance))) then
                    fit%status = irregular_status_singular
                    cycle
                end if
                fit%residual_variance(k) = residual_scale
                fit%covariance(2:,2:,k) = harmonic_covariance
                if (floating) then
                    ! Orthogonality of the intercept to the centered design
                    ! gives the usual full covariance without forming X'WX.
                    temp = matmul(harmonic_covariance,means)
                    fit%covariance(1,1,k) = residual_scale + dot_product(means,temp)
                    fit%covariance(1,2:,k) = -temp
                    fit%covariance(2:,1,k) = -temp
                    fit%standard_error(1,k) = sqrt(max(0.0_dp,fit%covariance(1,1,k)))
                end if
                do i = 1, nparam
                    fit%standard_error(i+1,k) = sqrt(max(0.0_dp,harmonic_covariance(i,i)))
                end do
            end if
            fit%rms(k) = sqrt(error_sum)
            fit%valid(k) = .true.
        end do
    end function lombscargle_harmonic_fit

    ! Evaluate one frequency's fitted harmonic series at arbitrary times.
    ! The returned NaNs preserve invalidity without forcing a STOP.
    function lombscargle_harmonic_predict(fit, xnew, frequency_index, status) result(prediction)
        type(harmonic_fit_result), intent(in) :: fit !! fitted harmonic coefficients
        real(dp), intent(in) :: xnew(:) !! new observation times, including unsorted times
        integer, intent(in) :: frequency_index !! one-based trial frequency index
        integer, intent(out), optional :: status !! success or invalid input
        real(dp), allocatable :: prediction(:)
        integer :: h
        real(dp) :: theta

        allocate(prediction(size(xnew)))
        prediction = ieee_value(0.0_dp, ieee_quiet_nan)
        if (present(status)) status = irregular_status_invalid_input
        if (.not. allocated(fit%valid)) return
        if (.not. allocated(fit%frequencies)) return
        if (.not. allocated(fit%offset)) return
        if (.not. allocated(fit%cosine)) return
        if (.not. allocated(fit%sine)) return
        if (fit%nterms < 1) return
        if (frequency_index < 1 .or. frequency_index > size(fit%valid)) return
        if (.not. fit%valid(frequency_index)) return
        if (size(fit%frequencies) /= size(fit%valid)) return
        if (size(fit%offset) /= size(fit%valid)) return
        if (size(fit%cosine,1) /= fit%nterms .or. size(fit%sine,1) /= fit%nterms) return
        if (size(fit%cosine,2) /= size(fit%valid)) return
        if (size(fit%sine,2) /= size(fit%valid)) return
        if (any(.not. ieee_is_finite(xnew))) return
        prediction = fit%offset(frequency_index)
        do h = 1, fit%nterms
            theta = real(h,dp)*fit%frequencies(frequency_index)
            prediction = prediction + fit%cosine(h,frequency_index)*cos(theta*xnew) + &
                fit%sine(h,frequency_index)*sin(theta*xnew)
        end do
        if (any(.not. ieee_is_finite(prediction))) then
            prediction = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if
        if (present(status)) status = irregular_status_success
    end function lombscargle_harmonic_predict

    ! Standard error of the fitted MEAN curve, not a prediction interval for
    ! a new noisy observation.  Errors are conditional on fixed frequencies
    ! and, when precenter is used without a floating mean, on its fixed shift.
    function lombscargle_harmonic_predict_std(fit, xnew, frequency_index, status) result(standard_error)
        type(harmonic_fit_result), intent(in) :: fit !! result fitted with estimate_covariance=.true.
        real(dp), intent(in) :: xnew(:) !! finite prediction times
        integer, intent(in) :: frequency_index !! one-based frequency index
        integer, intent(out), optional :: status !! zero for successful uncertainty evaluation
        real(dp), allocatable :: standard_error(:)

        real(dp), allocatable :: basis(:), product(:), evaluated(:)
        real(dp) :: theta, variance
        integer :: h, j, p, k

        allocate(standard_error(size(xnew)))
        standard_error = ieee_value(0.0_dp, ieee_quiet_nan)
        if (present(status)) status = irregular_status_invalid_input
        if (.not. allocated(fit%covariance) .or. .not. allocated(fit%valid)) return
        if (.not. allocated(fit%frequencies)) return
        if (fit%nterms < 1) return
        k = frequency_index
        if (k < 1 .or. k > size(fit%valid)) return
        if (.not. fit%valid(k)) return
        if (size(fit%frequencies) /= size(fit%valid)) return
        p = 1+2*fit%nterms
        if (any(shape(fit%covariance) /= [p,p,size(fit%valid)])) return
        if (any(.not. ieee_is_finite(xnew))) return
        allocate(basis(p),product(p),evaluated(size(xnew)))
        do j = 1, size(xnew)
            basis(1) = 0.0_dp
            if (fit%intercept_fitted) basis(1) = 1.0_dp
            do h = 1, fit%nterms
                theta = real(h,dp)*fit%frequencies(k)*xnew(j)
                if (.not. ieee_is_finite(theta)) return
                basis(2*h) = cos(theta)
                basis(2*h+1) = sin(theta)
            end do
            if (fit%intercept_fitted) then
                if (any(.not. ieee_is_finite(fit%covariance(:,:,k)))) return
                product = matmul(fit%covariance(:,:,k),basis)
                variance = dot_product(basis,product)
            else
                if (any(.not. ieee_is_finite(fit%covariance(2:,2:,k)))) return
                product(2:) = matmul(fit%covariance(2:,2:,k),basis(2:))
                variance = dot_product(basis(2:),product(2:))
            end if
            if (.not. ieee_is_finite(variance)) return
            ! Allow tiny roundoff below zero for a positive semidefinite
            ! covariance matrix; larger negatives indicate a failed fit.
            if (variance < -64.0_dp*epsilon(1.0_dp)* &
                max(1.0_dp,maxval(abs(fit%covariance(2:,2:,k))))) return
            evaluated(j) = sqrt(max(0.0_dp,variance))
        end do
        standard_error = evaluated
        if (present(status)) status = irregular_status_success
    end function lombscargle_harmonic_predict_std

    ! Build the time- and frequency-dependent QR factors exactly once.
    ! Subsequent series require only weighted projections and triangular
    ! solves, rather than new harmonic basis and QR factorizations.
    function lombscargle_harmonic_prepare(x, freqs, nterms, weights, floating_mean) result(plan)
        real(dp), intent(in) :: x(:) !! finite irregular sample times
        real(dp), intent(in) :: freqs(:) !! finite angular trial frequencies
        integer, intent(in) :: nterms !! positive number of harmonics
        real(dp), intent(in), optional :: weights(:) !! nonnegative weights, common across channels
        logical, intent(in), optional :: floating_mean !! fit a weighted intercept, default false
        type(harmonic_fit_plan) :: plan

        integer :: n, m, p, k, col, h, j, pass
        real(dp) :: maxweight, theta, original_norm, norm_v, projection, tol
        real(dp), allocatable :: v(:)
        logical :: rank_ok

        n = size(x)
        m = size(freqs)
        plan%nobs = n
        allocate(plan%frequencies(m), plan%valid(m))
        plan%frequencies = freqs
        plan%valid = .false.
        if (present(floating_mean)) plan%floating_mean = floating_mean
        if (nterms >= 1 .and. nterms <= n/2) plan%nterms = nterms
        if (nterms < 1 .or. nterms > (huge(nterms)-1)/2) return
        p = 2*nterms
        if (n < 3 .or. m == 0) return
        if (p+merge(1,0,plan%floating_mean) >= n) return
        if (any(.not. ieee_is_finite(x)) .or. any(.not. ieee_is_finite(freqs))) return
        allocate(plan%weights(n))
        if (present(weights)) then
            if (size(weights) /= n) return
            if (any(.not. ieee_is_finite(weights)) .or. any(weights < 0.0_dp)) return
            maxweight = maxval(weights)
            if (maxweight <= 0.0_dp) return
            plan%weights = weights/maxweight
        else
            plan%weights = 1.0_dp
        end if
        plan%dof = count(plan%weights > 0.0_dp)-p-merge(1,0,plan%floating_mean)
        if (plan%dof <= 0) return
        plan%weights = plan%weights/sum(plan%weights)
        allocate(plan%root_weights(n))
        plan%root_weights = sqrt(plan%weights)
        allocate(plan%q(n,p,m), plan%r(p,p,m), plan%column_means(p,m), v(n))
        plan%q = 0.0_dp
        plan%r = 0.0_dp
        plan%column_means = 0.0_dp
        plan%status = irregular_status_success
        tol = 64.0_dp*epsilon(1.0_dp)*real(p,dp)
        do k = 1, m
            rank_ok = .true.
            do col = 1, p
                h = (col+1)/2
                do j = 1, n
                    theta = real(h,dp)*freqs(k)*x(j)
                    if (.not. ieee_is_finite(theta)) then
                        rank_ok = .false.
                        exit
                    end if
                    if (mod(col,2) == 1) then
                        v(j) = cos(theta)
                    else
                        v(j) = sin(theta)
                    end if
                end do
                if (.not. rank_ok) exit
                if (plan%floating_mean) then
                    plan%column_means(col,k) = sum(plan%weights*v)
                    v = v-plan%column_means(col,k)
                end if
                v = plan%root_weights*v
                original_norm = sqrt(dot_product(v,v))
                do pass = 1, 2
                    do j = 1, col-1
                        projection = dot_product(plan%q(:,j,k),v)
                        plan%r(j,col,k) = plan%r(j,col,k)+projection
                        v = v-projection*plan%q(:,j,k)
                    end do
                end do
                norm_v = sqrt(dot_product(v,v))
                if (.not. ieee_is_finite(norm_v)) then
                    rank_ok = .false.
                    exit
                end if
                if (norm_v <= tol .or. norm_v <= tol*original_norm) then
                    rank_ok = .false.
                    exit
                end if
                plan%r(col,col,k) = norm_v
                plan%q(:,col,k) = v/norm_v
            end do
            plan%valid(k) = rank_ok
            if (.not. rank_ok) plan%status = irregular_status_singular
        end do
    end function lombscargle_harmonic_prepare

    ! Fit one channel using the prepared QR factors, returning the same
    ! coefficient, power, and uncertainty fields as the ordinary fit.
    function lombscargle_harmonic_fit_prepared(plan, y, precenter, estimate_covariance) result(fit)
        type(harmonic_fit_plan), intent(in) :: plan !! previously prepared frequency grid
        real(dp), intent(in) :: y(:) !! observations at the plan's sample times
        logical, intent(in), optional :: precenter !! subtract unweighted mean before fitting
        logical, intent(in), optional :: estimate_covariance !! compute covariance and mean-curve SE
        type(harmonic_fit_result) :: fit

        real(dp), allocatable :: target(:), centered(:), coefficients(:), z(:), residual(:)
        real(dp), allocatable :: inverse_r(:,:), harmonic_covariance(:,:), temp(:)
        real(dp) :: nan, shift, mean_y, yy, improvement, error_sum, intercept
        real(dp) :: residual_scale, value
        integer :: m, n, p, k, col, h, i
        logical :: calculate_covariance

        m = 0
        if (allocated(plan%frequencies)) m = size(plan%frequencies)
        fit%nterms = plan%nterms
        fit%intercept_fitted = plan%floating_mean
        fit%degrees_of_freedom = plan%dof
        nan = ieee_value(0.0_dp, ieee_quiet_nan)
        allocate(fit%frequencies(m), fit%offset(m), fit%cosine(fit%nterms,m), &
            fit%sine(fit%nterms,m), fit%power(m), fit%normalized_power(m), &
            fit%rms(m), fit%valid(m))
        if (m > 0) fit%frequencies = plan%frequencies
        fit%offset = nan
        fit%cosine = nan
        fit%sine = nan
        fit%power = nan
        fit%normalized_power = nan
        fit%rms = nan
        fit%valid = .false.
        fit%status = irregular_status_invalid_input
        calculate_covariance = .false.
        if (present(estimate_covariance)) calculate_covariance = estimate_covariance
        p = 2*fit%nterms
        if (calculate_covariance) then
            allocate(fit%covariance(p+1,p+1,m), fit%standard_error(p+1,m), &
                fit%residual_variance(m))
            fit%covariance = nan
            fit%standard_error = nan
            fit%residual_variance = nan
        end if
        if (plan%status /= irregular_status_success .and. &
            plan%status /= irregular_status_singular) return
        n = plan%nobs
        if (size(y) /= n .or. any(.not. ieee_is_finite(y))) return
        centered = y
        shift = 0.0_dp
        if (present(precenter)) then
            if (precenter) then
                shift = sum(y)/real(n,dp)
                centered = centered-shift
            end if
        end if
        mean_y = 0.0_dp
        if (plan%floating_mean) then
            mean_y = sum(plan%weights*centered)
            centered = centered-mean_y
        end if
        target = plan%root_weights*centered
        yy = dot_product(target,target)
        if (plan%floating_mean .and. maxval(abs(y-y(1))) == 0.0_dp) yy = 0.0_dp
        if (.not. ieee_is_finite(yy)) return
        if (yy <= 0.0_dp) then
            fit%status = irregular_status_singular
            return
        end if
        allocate(coefficients(p),z(p),residual(n))
        if (calculate_covariance) then
            allocate(inverse_r(p,p),harmonic_covariance(p,p),temp(p))
        end if
        fit%status = plan%status
        do k = 1, m
            if (.not. plan%valid(k)) cycle
            do col = 1, p
                z(col) = dot_product(plan%q(:,col,k),target)
            end do
            coefficients = z
            do col = p, 1, -1
                if (col < p) coefficients(col) = coefficients(col)- &
                    dot_product(plan%r(col,col+1:p,k),coefficients(col+1:p))
                coefficients(col) = coefficients(col)/plan%r(col,col,k)
            end do
            if (any(.not. ieee_is_finite(coefficients))) then
                fit%status = irregular_status_singular
                cycle
            end if
            intercept = shift+mean_y
            if (plan%floating_mean) intercept = intercept- &
                dot_product(plan%column_means(:,k),coefficients)
            if (.not. ieee_is_finite(intercept)) then
                fit%status = irregular_status_singular
                cycle
            end if
            ! Q spans the weighted, mean-adjusted harmonic design.  Computing
            ! residuals in Q space avoids repeated trigonometric evaluations.
            residual = target-matmul(plan%q(:,:,k),z)
            error_sum = dot_product(residual,residual)
            if (.not. ieee_is_finite(error_sum)) then
                fit%status = irregular_status_singular
                cycle
            end if
            if (calculate_covariance) then
                residual_scale = max(0.0_dp,error_sum)/real(fit%degrees_of_freedom,dp)
                if (.not. ieee_is_finite(residual_scale)) then
                    fit%status = irregular_status_singular
                    cycle
                end if
                inverse_r = 0.0_dp
                do col = 1, p
                    do i = p, 1, -1
                        value = 0.0_dp
                        if (i == col) value = 1.0_dp
                        if (i < p) value = value- &
                            dot_product(plan%r(i,i+1:p,k),inverse_r(i+1:p,col))
                        inverse_r(i,col) = value/plan%r(i,i,k)
                    end do
                end do
                harmonic_covariance = residual_scale*matmul(inverse_r,transpose(inverse_r))
                if (any(.not. ieee_is_finite(harmonic_covariance))) then
                    fit%status = irregular_status_singular
                    cycle
                end if
                fit%residual_variance(k) = residual_scale
                fit%covariance(2:,2:,k) = harmonic_covariance
                if (plan%floating_mean) then
                    temp = matmul(harmonic_covariance,plan%column_means(:,k))
                    fit%covariance(1,1,k) = residual_scale+dot_product(plan%column_means(:,k),temp)
                    fit%covariance(1,2:,k) = -temp
                    fit%covariance(2:,1,k) = -temp
                    fit%standard_error(1,k) = sqrt(max(0.0_dp,fit%covariance(1,1,k)))
                end if
                do i = 1, p
                    fit%standard_error(i+1,k) = sqrt(max(0.0_dp,harmonic_covariance(i,i)))
                end do
            end if
            do h = 1, fit%nterms
                fit%cosine(h,k) = coefficients(2*h-1)
                fit%sine(h,k) = coefficients(2*h)
            end do
            fit%offset(k) = intercept
            improvement = min(yy,max(0.0_dp,dot_product(z,z)))
            fit%power(k) = 0.5_dp*real(n,dp)*improvement
            fit%normalized_power(k) = improvement/yy
            fit%rms(k) = sqrt(error_sum)
            fit%valid(k) = .true.
        end do
    end function lombscargle_harmonic_fit_prepared

    ! Channel j occupies y(:,j).  The shared plan is never modified, so it
    ! can be reused for any number of calls or independently held channels.
    function lombscargle_harmonic_fit_batch(plan, y, precenter, estimate_covariance) result(fits)
        type(harmonic_fit_plan), intent(in) :: plan !! reusable weighted harmonic design
        real(dp), intent(in) :: y(:,:) !! input observations (nsamples,nchannels)
        logical, intent(in), optional :: precenter !! center each channel independently
        logical, intent(in), optional :: estimate_covariance !! estimate each channel's covariance
        type(harmonic_fit_result), allocatable :: fits(:)
        integer :: j

        allocate(fits(size(y,2)))
        do j = 1, size(y,2)
            fits(j) = lombscargle_harmonic_fit_prepared(plan,y(:,j),precenter,estimate_covariance)
        end do
    end function lombscargle_harmonic_fit_batch

    function vectorstrength_array(events, periods) result(out)
        real(dp), intent(in) :: events(:) !! event times, nonempty
        real(dp), intent(in) :: periods(:) !! positive test periods, nonempty
        type(vector_strength_result) :: out

        real(dp) :: cs, sn, theta
        integer :: i, j

        allocate(out%strength(size(periods)), out%phase(size(periods)))
        out%strength = ieee_value(0.0_dp, ieee_quiet_nan)
        out%phase = ieee_value(0.0_dp, ieee_quiet_nan)
        if (size(events) == 0 .or. size(periods) == 0) return
        if (any(.not. ieee_is_finite(events))) return
        if (any(.not. ieee_is_finite(periods))) return
        if (any(periods <= 0.0_dp)) return

        do j = 1, size(periods)
            cs = 0.0_dp
            sn = 0.0_dp
            do i = 1, size(events)
                theta = twopi * (events(i) / periods(j))
                cs = cs + cos(theta)
                sn = sn + sin(theta)
            end do
            cs = cs / real(size(events), dp)
            sn = sn / real(size(events), dp)
            out%strength(j) = sqrt(cs*cs+sn*sn)
            out%phase(j) = atan2(sn, cs)
        end do
        out%status = irregular_status_success
    end function vectorstrength_array

    function vectorstrength_scalar(events, period) result(out)
        real(dp), intent(in) :: events(:) !! event times
        real(dp), intent(in) :: period !! single positive test period
        type(vector_strength_result) :: out

        out = vectorstrength_array(events, [period])
    end function vectorstrength_scalar

end module scifort_signal_irregular
