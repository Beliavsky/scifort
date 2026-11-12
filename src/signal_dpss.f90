! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Discrete prolate spheroidal sequences (Slepian tapers).
! Independent implementation of the tridiagonal eigenvalue formulation
! described by Slepian (1978) and Percival and Walden (1993).
! Uses either implicit-shift QL or selected eigenpairs from tridiagonal
! Sturm bisection and inverse iteration, without external LAPACK.
module scifort_signal_dpss
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    implicit none
    private

    integer, parameter, public :: dpss_status_success = 0
    integer, parameter, public :: dpss_status_invalid_input = 1
    integer, parameter, public :: dpss_status_no_convergence = 2
    real(dp), parameter :: pi = acos(-1.0_dp)

    type, public :: dpss_result
        real(dp), allocatable :: tapers(:,:) !! (sample index, taper index), one-based
        real(dp), allocatable :: ratios(:) !! spectral concentration, one per taper
        integer :: status = dpss_status_invalid_input
    end type dpss_result

    public :: dpss
contains

    function dpss(m, nw, kmax, sym, norm, solver, ratio_method) result(out)
        integer, intent(in) :: m !! output window length
        real(dp), intent(in) :: nw !! positive time-half-bandwidth product; NW < M/2
        integer, intent(in), optional :: kmax !! number of tapers, default one
        logical, intent(in), optional :: sym !! symmetric (default) or periodic
        character(*), intent(in), optional :: norm !! 'l2', 'approximate', or 'subsample'; default l2 with kmax
        character(*), intent(in), optional :: solver !! 'auto' (default), 'ql', or 'partial'
        character(*), intent(in), optional :: ratio_method !! 'auto' (default), 'direct', or 'fft'
        type(dpss_result) :: out

        real(dp), allocatable :: d(:), e(:), eigenvectors(:,:), taper(:)
        complex(dp), allocatable :: spectrum(:)
        real(dp) :: w, idx, correction, vmax, vfirst, ratio, center
        real(dp) :: tmp
        integer :: n, k, i, j, source, it, ierr, fft_length
        logical :: symmetric, use_partial, use_fft_ratio
        character(:), allocatable :: mode, eigensolver, ratio_algorithm

        allocate(out%tapers(max(0,m),0), out%ratios(0))
        if (m < 0 .or. .not. ieee_is_finite(nw)) return
        if (m == 0) then
            out%status = dpss_status_success
            return
        end if
        k = 1
        if (present(kmax)) k = kmax
        if (k < 1 .or. k > m) return
        mode = 'approximate'
        if (present(kmax)) mode = 'l2'
        if (present(norm)) mode = trim(norm)
        if (mode /= 'l2' .and. mode /= 'approximate' .and. mode /= 'subsample') return
        eigensolver = 'auto'
        if (present(solver)) eigensolver = trim(solver)
        if (eigensolver /= 'auto' .and. eigensolver /= 'ql' .and. &
            eigensolver /= 'partial') return
        ratio_algorithm = 'auto'
        if (present(ratio_method)) ratio_algorithm = trim(ratio_method)
        if (ratio_algorithm /= 'auto' .and. ratio_algorithm /= 'direct' .and. &
            ratio_algorithm /= 'fft') return
        if (m > 1) then
            if (nw <= 0.0_dp .or. nw >= 0.5_dp * real(m,dp)) return
        end if
        symmetric = .true.
        if (present(sym)) symmetric = sym
        n = m
        if (.not. symmetric .and. m > 1) n = m + 1
        deallocate(out%tapers, out%ratios)
        allocate(out%tapers(m,k), out%ratios(k))
        out%tapers = 0.0_dp
        out%ratios = 0.0_dp
        if (m == 1) then
            out%tapers = 1.0_dp
            out%ratios = 1.0_dp
            out%status = dpss_status_success
            return
        end if

        ! For short vectors direct correlation avoids FFT setup; on long
        ! tapers the fast linear convolution avoids quadratic work.
        use_fft_ratio = ratio_algorithm == 'fft'
        if (ratio_algorithm == 'auto') use_fft_ratio = n >= 256
        if (use_fft_ratio) then
            ! An exact linear autocorrelation requires at least 2*n-1 bins.
            if (n > (huge(n)-1)/2) then
                out%status = dpss_status_no_convergence
                return
            end if
            fft_length = 1
            do while (fft_length < 2*n-1)
                if (fft_length > huge(fft_length)/2) then
                    out%status = dpss_status_no_convergence
                    return
                end if
                fft_length = 2*fft_length
            end do
            allocate(spectrum(fft_length), stat=ierr)
            if (ierr /= 0) then
                out%status = dpss_status_no_convergence
                return
            end if
        end if
        use_partial = eigensolver == 'partial'
        if (eigensolver == 'auto') use_partial = n >= 256 .and. k <= 32
        allocate(d(n), e(n), taper(n))
        w = nw / real(n, dp)
        do i = 1, n
            idx = 0.5_dp * real(n - 2*i + 1,dp)
            d(i) = idx * idx * cos(2.0_dp*pi*w)
            e(i) = 0.0_dp
            if (i < n) e(i) = 0.5_dp * real(i,dp) * real(n-i,dp)
        end do
        if (use_partial) then
            allocate(eigenvectors(n,k))
            call tridiagonal_selected_eigenvectors(d, e, eigenvectors, ierr)
        else
            allocate(eigenvectors(n,n))
            eigenvectors = 0.0_dp
            do i = 1, n
                eigenvectors(i,i) = 1.0_dp
            end do
            call tridiagonal_ql(d, e, eigenvectors, ierr)
        end if
        if (ierr /= 0) then
            out%status = dpss_status_no_convergence
            out%tapers = 0.0_dp
            out%ratios = 0.0_dp
            return
        end if
        ! Extract descending eigenvalue order. The QL algorithm returns
        ! unsorted, l2-normalized eigenvectors.
        do j = 1, k
            if (.not. use_partial) then
                source = maxloc(d(j:n),dim=1) + j - 1
                tmp = d(j)
                d(j) = d(source)
                d(source) = tmp
                if (source /= j) then
                    taper = eigenvectors(:,j)
                    eigenvectors(:,j) = eigenvectors(:,source)
                    eigenvectors(:,source) = taper
                end if
            end if
            taper = eigenvectors(:,j)
            ! SciPy/Percival-Walden sign conventions: even orders have
            ! positive sum; odd orders start with a positive main lobe.
            if (mod(j,2) == 1) then
                if (sum(taper) < 0.0_dp) taper = -taper
            else
                vfirst = max(1.0e-7_dp, 1.0_dp / real(n,dp))
                it = maxloc(abs(taper),dim=1)
                do i = 1, n
                    if (taper(i)**2 > vfirst) then
                        it = i
                        exit
                    end if
                end do
                if (taper(it) < 0.0_dp) taper = -taper
            end if
            ! The Toeplitz sinc Rayleigh quotient uses the *full*
            ! symmetric taper, even when the returned window is periodic.
            if (use_fft_ratio) then
                ratio = concentration_ratio_fft(taper, w, spectrum)
            else
                ratio = concentration_ratio_direct(taper, w)
            end if
            if (.not. ieee_is_finite(ratio)) then
                out%status = dpss_status_no_convergence
                return
            end if
            out%ratios(j) = min(1.0_dp, max(0.0_dp, ratio))
            eigenvectors(:,j) = taper
        end do

        correction = 1.0_dp
        if (mode /= 'l2') then
            vmax = maxval(eigenvectors(:,1:k))
            if (vmax <= 0.0_dp) then
                out%status = dpss_status_no_convergence
                return
            end if
            correction = 1.0_dp / vmax
            if (mod(n,2) == 0) then
                if (mode == 'approximate') then
                    correction = correction * real(n,dp)**2 / (real(n,dp)**2 + nw)
                else
                    ! Peak of the bandlimited periodic interpolation halfway
                    ! between the two central samples.  The first taper
                    ! defines SciPy's common 'subsample' scale for all k.
                    center = subsample_center_value(eigenvectors(:,1))
                    if (.not. ieee_is_finite(center) .or. center <= 0.0_dp) then
                        out%status = dpss_status_no_convergence
                        return
                    end if
                    correction = 1.0_dp / center
                end if
            end if
        end if
        out%tapers = eigenvectors(1:m,1:k) * correction
        out%status = dpss_status_success
    end function dpss


    ! Spectral concentration is the Rayleigh quotient of the symmetric
    ! Toeplitz matrix C(i,j) = sin(2*pi*W*(i-j))/(pi*(i-j)), diagonal 2*W.
    ! The direct O(n*n) calculation is retained as a small-window reference.
    pure real(dp) function concentration_ratio_direct(taper, w) result(ratio)
        real(dp), intent(in) :: taper(:) !! full l2-normalized DPSS eigenvector
        real(dp), intent(in) :: w !! normalized half bandwidth, 0 < W < 1/2
        real(dp) :: rxx, rlag
        integer :: n, lag

        n = size(taper)
        ratio = 2.0_dp*w*sum(taper*taper)
        do lag = 1, n-1
            rxx = sum(taper(1:n-lag)*taper(1+lag:n))
            rlag = 2.0_dp*sin(2.0_dp*pi*w*real(lag,dp))/(pi*real(lag,dp))
            ratio = ratio + rlag*rxx
        end do
    end function concentration_ratio_direct

    ! FFT autocorrelation R(lag) = sum_i taper(i)*taper(i+lag).
    ! Zero-padding to >= 2*n-1 prevents circular-aliasing of the correlation.
    ! The workspace is reused for successive tapers. A private radix-2 kernel
    ! avoids the signal_spectral -> signal_dpss module dependency cycle.
    real(dp) function concentration_ratio_fft(taper, w, work) result(ratio)
        real(dp), intent(in) :: taper(:) !! full l2-normalized DPSS eigenvector
        real(dp), intent(in) :: w !! normalized half bandwidth, 0 < W < 1/2
        complex(dp), intent(inout) :: work(:) !! reusable power-of-two Fourier scratch
        real(dp) :: rlag, term, corrected, compensation, next_sum
        integer :: n, lag

        n = size(taper)
        work = (0.0_dp,0.0_dp)
        work(1:n) = cmplx(taper,0.0_dp,dp)
        call dpss_radix2_fft(work,.false.)
        work = cmplx(real(work,dp)**2 + aimag(work)**2,0.0_dp,dp)
        call dpss_radix2_fft(work,.true.)
        ratio = 2.0_dp*w*sum(taper*taper)
        compensation = 0.0_dp
        do lag = 1, n-1
            rlag = 2.0_dp*sin(2.0_dp*pi*w*real(lag,dp))/(pi*real(lag,dp))
            term = rlag*real(work(lag+1),dp)
            corrected = term - compensation
            next_sum = ratio + corrected
            compensation = (next_sum - ratio) - corrected
            ratio = next_sum
        end do
    end function concentration_ratio_fft

    ! In-place radix-2 FFT, with a normalized inverse. A twiddle recurrence
    ! refreshed every 64 butterflies bounds long-stage phase drift.
    subroutine dpss_radix2_fft(values, inverse)
        complex(dp), intent(inout) :: values(:) !! samples in; transform out
        logical, intent(in) :: inverse !! true for normalized inverse transform
        complex(dp) :: tmp, twiddle, root
        real(dp) :: angle
        integer :: n, i, j, bit, stage, half, first, t

        n = size(values)
        j = 1
        do i = 1, n-1
            if (i < j) then
                tmp = values(i)
                values(i) = values(j)
                values(j) = tmp
            end if
            bit = n/2
            do while (bit > 0)
                if (j <= bit) exit
                j = j-bit
                bit = bit/2
            end do
            j = j+bit
        end do
        stage = 2
        do while (stage <= n)
            half = stage/2
            angle = -2.0_dp*pi/real(stage,dp)
            if (inverse) angle = -angle
            root = cmplx(cos(angle),sin(angle),dp)
            do first = 1, n, stage
                twiddle = (1.0_dp,0.0_dp)
                do t = 0, half-1
                    if (mod(t,64) == 0) then
                        twiddle = cmplx(cos(angle*real(t,dp)), &
                            sin(angle*real(t,dp)),dp)
                    end if
                    tmp = twiddle*values(first+t+half)
                    values(first+t+half) = values(first+t)-tmp
                    values(first+t) = values(first+t)+tmp
                    twiddle = twiddle*root
                end do
            end do
            if (stage == n) exit
            stage = 2*stage
        end do
        if (inverse) values = values/real(n,dp)
    end subroutine dpss_radix2_fft

    ! Exact finite-Fourier-series interpolation at index (n-1)/2 of an
    ! even-length, real, symmetric taper. The Dirichlet-kernel identity
    !   1 + 2 * sum_{k=1}^{n/2} cos(k*theta)
    !      = sin((n+1)*theta/2) / sin(theta/2)
    ! eliminates the FFT and evaluates SciPy's half-sample correction in
    ! O(n) arithmetic. For even n the denominator is never zero because
    ! n-2*j+1 is odd. The Nyquist cosine vanishes by symmetry.
    pure real(dp) function subsample_center_value(taper) result(value)
        real(dp), intent(in) :: taper(:) !! positive-centered first DPSS taper, even length
        real(dp) :: theta, weight, accum, compensation, term, next_sum
        integer :: n, j

        n = size(taper)
        accum = 0.0_dp
        compensation = 0.0_dp
        do j = 1, n
            theta = pi * real(n - 2*j + 1,dp) / real(n,dp)
            weight = sin(0.5_dp*(real(n,dp)+1.0_dp)*theta) / sin(0.5_dp*theta)
            term = taper(j)*weight - compensation
            next_sum = accum + term
            compensation = (next_sum - accum) - term
            accum = next_sum
        end do
        value = accum / real(n,dp)
    end function subsample_center_value

    ! Compute only the requested *largest* eigenvectors of a symmetric
    ! tridiagonal matrix. This avoids the O(n*n) storage used by QL.
    ! Bisection finds each eigenvalue using Sturm inertia counts, followed by
    ! inverse iteration on a pivoted tridiagonal LU factorization. Two rounds
    ! of Gram-Schmidt protect the ordering of closely spaced eigenvectors.
    subroutine tridiagonal_selected_eigenvectors(diagonal, offdiag, vectors, status)
        real(dp), intent(in) :: diagonal(:) !! tridiagonal main diagonal, length n
        real(dp), intent(in) :: offdiag(:) !! upper diagonal, first n-1 elements used
        real(dp), intent(out) :: vectors(:,:) !! n by k selected eigenvectors, descending eigenvalues
        integer, intent(out) :: status !! zero on convergence, one otherwise
        real(dp), allocatable :: a(:), b(:), rhs(:), residual(:)
        real(dp), allocatable :: ld(:), ud(:), ud2(:), pivdiag(:)
        integer, allocatable :: piv(:)
        real(dp) :: scale, low, high, left, right, middle, eigenvalue
        real(dp) :: row_sum, norm, proj, shift, err
        integer :: n, k, i, j, h, iteration, sweep, target, bisect_step
        logical :: converged

        status = 1
        n = size(diagonal)
        k = size(vectors,2)
        if (n < 2 .or. k < 1 .or. k > n) return
        if (size(offdiag) < n-1 .or. size(vectors,1) /= n) return
        scale = max(maxval(abs(diagonal)), maxval(abs(offdiag(1:n-1))))
        if (.not. ieee_is_finite(scale) .or. scale <= 0.0_dp) return
        allocate(a(n), b(n-1), rhs(n), residual(n))
        allocate(ld(n-1), ud(n-1), ud2(max(0,n-2)), pivdiag(n), piv(n-1))
        a = diagonal / scale
        b = offdiag(1:n-1) / scale
        low = huge(1.0_dp)
        high = -huge(1.0_dp)
        do i = 1, n
            row_sum = 0.0_dp
            if (i > 1) row_sum = row_sum + abs(b(i-1))
            if (i < n) row_sum = row_sum + abs(b(i))
            low = min(low, a(i)-row_sum)
            high = max(high, a(i)+row_sum)
        end do
        low = low - 16.0_dp*epsilon(1.0_dp)
        high = high + 16.0_dp*epsilon(1.0_dp)

        do j = 1, k
            target = n-j+1
            left = low
            right = high
            do bisect_step = 1, 80
                middle = left + 0.5_dp*(right-left)
                if (middle <= left .or. middle >= right) exit
                if (sturm_negative_pivots(a, b, middle) < target) then
                    left = middle
                else
                    right = middle
                end if
            end do
            eigenvalue = left + 0.5_dp*(right-left)
            ! Slightly offset the shift so singular LU factors do not divide
            ! by zero. The perturbation is below the requested accuracy.
            shift = eigenvalue + 8.0_dp*epsilon(1.0_dp)
            call tridiagonal_shift_lu(a, b, shift, ld, ud, ud2, pivdiag, piv)

            do i = 1, n
                rhs(i) = sin(0.037_dp*real(i,dp)*real(j+1,dp)) + &
                    cos(0.071_dp*real(i,dp)*real(2*j+5,dp))
            end do
            converged = .false.
            do iteration = 1, 24
                call tridiagonal_shift_solve(ld, ud, ud2, pivdiag, piv, rhs)
                if (.not. all(ieee_is_finite(rhs))) exit
                do sweep = 1, 2
                    do h = 1, j-1
                        proj = dot_product(rhs, vectors(:,h))
                        rhs = rhs - proj*vectors(:,h)
                    end do
                end do
                norm = maxval(abs(rhs))
                if (.not. ieee_is_finite(norm) .or. norm <= tiny(norm)) exit
                rhs = rhs/norm
                norm = sqrt(sum(rhs*rhs))
                if (norm <= tiny(norm)) exit
                rhs = rhs/norm
                residual = (a-eigenvalue)*rhs
                residual(1) = residual(1) + b(1)*rhs(2)
                do i = 2, n-1
                    residual(i) = residual(i) + b(i-1)*rhs(i-1) + b(i)*rhs(i+1)
                end do
                residual(n) = residual(n) + b(n-1)*rhs(n-1)
                err = sqrt(sum(residual*residual))
                if (err <= 2.0e-13_dp) then
                    converged = .true.
                    exit
                end if
            end do
            if (.not. converged) return
            vectors(:,j) = rhs
        end do
        status = 0
    end subroutine tridiagonal_selected_eigenvectors

    ! Number of eigenvalues <= the trial point, by Sturm's theorem.
    pure integer function sturm_negative_pivots(a, b, shift) result(count)
        real(dp), intent(in) :: a(:) !! scaled diagonal
        real(dp), intent(in) :: b(:) !! scaled offdiagonal
        real(dp), intent(in) :: shift !! trial eigenvalue in scaled units
        real(dp) :: pivot, threshold
        integer :: i

        threshold = 16.0_dp*tiny(1.0_dp)
        count = 0
        pivot = a(1) - shift
        if (pivot <= 0.0_dp) count = 1
        do i = 2, size(a)
            if (abs(pivot) < threshold) pivot = -threshold
            pivot = a(i) - shift - b(i-1)*(b(i-1)/pivot)
            if (pivot <= 0.0_dp) count = count + 1
        end do
    end function sturm_negative_pivots

    ! Pivoted LU factorization of a shifted tridiagonal matrix; stores a
    ! second superdiagonal only when adjacent rows must be interchanged.
    subroutine tridiagonal_shift_lu(a, b, shift, ld, ud, ud2, pivdiag, piv)
        real(dp), intent(in) :: a(:) !! scaled main diagonal
        real(dp), intent(in) :: b(:) !! scaled symmetric offdiagonal
        real(dp), intent(in) :: shift !! inverse-iteration shift
        real(dp), intent(out) :: ld(:) !! L multipliers
        real(dp), intent(out) :: ud(:) !! U first superdiagonal
        real(dp), intent(out) :: ud2(:) !! U second superdiagonal
        real(dp), intent(out) :: pivdiag(:) !! U diagonal
        integer, intent(out) :: piv(:) !! chosen row for each elimination step
        real(dp) :: factor, temp, floor_pivot
        integer :: n, i

        n = size(a)
        ld = b
        ud = b
        ud2 = 0.0_dp
        pivdiag = a - shift
        floor_pivot = 16.0_dp*epsilon(1.0_dp)
        do i = 1, n-1
            piv(i) = i
            if (abs(pivdiag(i)) >= abs(ld(i))) then
                if (abs(pivdiag(i)) < floor_pivot) then
                    pivdiag(i) = sign(floor_pivot, pivdiag(i))
                end if
                factor = ld(i)/pivdiag(i)
                pivdiag(i+1) = pivdiag(i+1) - factor*ud(i)
                ld(i) = factor
            else
                factor = pivdiag(i)/ld(i)
                pivdiag(i) = ld(i)
                temp = pivdiag(i+1)
                pivdiag(i+1) = ud(i) - factor*temp
                ld(i) = factor
                ud(i) = temp
                if (i < n-1) then
                    ud2(i) = ud(i+1)
                    ud(i+1) = -factor*ud(i+1)
                end if
                piv(i) = i+1
            end if
        end do
        if (abs(pivdiag(n)) < floor_pivot) then
            pivdiag(n) = sign(floor_pivot, pivdiag(n))
        end if
    end subroutine tridiagonal_shift_lu

    subroutine tridiagonal_shift_solve(ld, ud, ud2, pivdiag, piv, rhs)
        real(dp), intent(in) :: ld(:) !! L multipliers
        real(dp), intent(in) :: ud(:) !! U first superdiagonal
        real(dp), intent(in) :: ud2(:) !! U second superdiagonal
        real(dp), intent(in) :: pivdiag(:) !! diagonal pivots
        integer, intent(in) :: piv(:) !! adjacent-row permutations
        real(dp), intent(inout) :: rhs(:) !! input vector, overwritten by solution
        real(dp) :: tmp
        integer :: i, n

        n = size(rhs)
        do i = 1, n-1
            if (piv(i) == i) then
                rhs(i+1) = rhs(i+1) - ld(i)*rhs(i)
            else
                tmp = rhs(i)
                rhs(i) = rhs(i+1)
                rhs(i+1) = tmp - ld(i)*rhs(i)
            end if
        end do
        rhs(n) = rhs(n)/pivdiag(n)
        rhs(n-1) = (rhs(n-1)-ud(n-1)*rhs(n))/pivdiag(n-1)
        do i = n-2, 1, -1
            rhs(i) = (rhs(i)-ud(i)*rhs(i+1)-ud2(i)*rhs(i+2))/pivdiag(i)
        end do
    end subroutine tridiagonal_shift_solve

    subroutine tridiagonal_ql(d, e, z, status)
        real(dp), intent(inout) :: d(:) !! diagonal, overwritten by eigenvalues
        real(dp), intent(inout) :: e(:) !! upper offdiagonal, e(n) initially zero
        real(dp), intent(inout) :: z(:,:) !! initial identity matrix, then eigenvectors
        integer, intent(out) :: status !! zero on convergence
        integer :: n, l, m, iter, i, row
        real(dp) :: g, r, s, c, p, f, b, tmp, scale

        status = 0
        n = size(d)
        do l = 1, n
            iter = 0
            do
                m = l
                do while (m < n)
                    scale = abs(d(m)) + abs(d(m+1))
                    if (abs(e(m)) <= epsilon(1.0_dp)*scale) exit
                    m = m + 1
                end do
                if (m == l) exit
                iter = iter + 1
                if (iter > 100) then
                    status = 1
                    return
                end if
                g = (d(l+1)-d(l)) / (2.0_dp*e(l))
                r = hypot(g,1.0_dp)
                g = d(m) - d(l) + e(l)/(g + sign(r,g))
                s = 1.0_dp
                c = 1.0_dp
                p = 0.0_dp
                do i = m-1, l, -1
                    f = s * e(i)
                    b = c * e(i)
                    if (abs(f) >= abs(g)) then
                        c = g / f
                        r = hypot(c,1.0_dp)
                        e(i+1) = f*r
                        s = 1.0_dp/r
                        c = c*s
                    else
                        s = f/g
                        r = hypot(s,1.0_dp)
                        e(i+1) = g*r
                        c = 1.0_dp/r
                        s = s*c
                    end if
                    g = d(i+1) - p
                    r = (d(i) - g)*s + 2.0_dp*c*b
                    p = s*r
                    d(i+1) = g + p
                    g = c*r-b
                    do row = 1, n
                        tmp = z(row,i+1)
                        z(row,i+1) = s*z(row,i) + c*tmp
                        z(row,i) = c*z(row,i) - s*tmp
                    end do
                end do
                d(l) = d(l)-p
                e(l) = g
                e(m) = 0.0_dp
            end do
        end do
    end subroutine tridiagonal_ql
end module scifort_signal_dpss
