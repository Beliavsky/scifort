! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Equiripple FIR design by a barycentric Remez exchange iteration.
! The implementation is independent of SciPy's C/C++ implementation.
module scifort_signal_remez
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    implicit none
    private

    integer, parameter, public :: remez_status_success = 0
    integer, parameter, public :: remez_status_invalid_input = 1
    integer, parameter, public :: remez_status_no_convergence = 2
    ! Extended intermediates avoid cancellation in high-order barycentric
    ! interpolation. Fall back to dp on compilers without this real kind.
    integer, parameter :: qp = selected_real_kind(30, 400)
    ! When binary128 is absent, x87-style extended precision (if supported)
    ! improves high-order exchanges substantially over plain real64.
    integer, parameter :: ep = selected_real_kind(18, 300)
    integer, parameter :: wp = merge(qp, merge(ep, dp, ep > 0), qp > 0)
    real(wp), parameter :: pi = acos(-1.0_wp)
    public :: remez

contains

    function remez(numtaps, bands, desired, weight, filter_type, maxiter, grid_density, fs, status, &
        iterations, ripple) result(taps)
        integer, intent(in) :: numtaps !! number of FIR coefficients
        real(dp), intent(in) :: bands(:) !! ascending band edge pairs in [0, fs/2]
        real(dp), intent(in) :: desired(:) !! constant gain or differentiation slope per band
        real(dp), intent(in), optional :: weight(:) !! positive error weight per band
        character(*), intent(in), optional :: filter_type !! 'bandpass' (default), 'hilbert', or 'differentiator'
        integer, intent(in), optional :: maxiter !! maximum exchanges, default 25
        integer, intent(in), optional :: grid_density !! grid points per tap per Nyquist interval, default 16
        real(dp), intent(in), optional :: fs !! sample frequency, default 1 (as in SciPy)
        integer, intent(out), optional :: status !! success, invalid input, or no convergence
        integer, intent(out), optional :: iterations !! number of completed exchange iterations
        real(dp), intent(out), optional :: ripple !! maximum weighted response error on the design grid
        real(dp), allocatable :: taps(:)
        real(wp), allocatable :: xx(:), target(:), ww(:), envelope(:), err(:)
        real(wp), allocatable :: nodes(:), values(:), lambdas(:), coef(:), wband(:)
        integer, allocatable :: ext(:), new_ext(:)
        real(wp) :: rate, nyquist, df, f, phase, v, scale, delta, minerr, maxerr, max_all
        integer :: nb, m, nex, gd, mit, total, b, j, nstep, idx, loop, st, ktype, nsize
        logical :: anti, converged, ok
        character(:), allocatable :: kind

        allocate(taps(0))
        if (present(status)) status = remez_status_invalid_input
        if (present(iterations)) iterations = 0
        if (present(ripple)) ripple = 0.0_wp
        rate = 1.0_wp
        if (present(fs)) rate = fs
        if (.not. ieee_is_finite(rate) .or. rate <= 0.0_wp) return
        nyquist = rate/2.0_wp
        nb = size(desired)
        if (numtaps < 2 .or. nb < 1 .or. size(bands) /= 2*nb) return
        if (.not. all(ieee_is_finite(bands)) .or. .not. all(ieee_is_finite(desired))) return
        if (any(bands < 0.0_wp) .or. any(bands > nyquist)) return
        do b = 1, nb
            if (bands(2*b) <= bands(2*b-1)) return
        end do
        do b = 2, nb
            if (bands(2*b-1) <= bands(2*b-2)) return
        end do
        allocate(wband(nb))
        wband = 1.0_wp
        if (present(weight)) then
            if (size(weight) /= nb) return
            if (.not. all(ieee_is_finite(weight))) return
            if (any(weight <= 0.0_wp)) return
            wband = weight
        end if
        mit = 25
        if (present(maxiter)) mit = maxiter
        gd = 16
        if (present(grid_density)) gd = grid_density
        if (mit < 1 .or. gd < 2) return
        if (numtaps > 2000000) return
        if (gd > (huge(gd)-100)/(2*(numtaps+1))) return

        kind = 'bandpass'
        if (present(filter_type)) kind = trim(filter_type)
        select case (kind)
        case ('bandpass')
            anti = .false.
        case ('hilbert', 'differentiator')
            anti = .true.
        case default
            return
        end select
        if (anti) then
            if (mod(numtaps,2) == 1) then
                ktype = 3
                m = (numtaps-1)/2
            else
                ktype = 4
                m = numtaps/2
            end if
        else
            if (mod(numtaps,2) == 1) then
                ktype = 1
                m = (numtaps+1)/2
            else
                ktype = 2
                m = numtaps/2
            end if
        end if
        nex = m+1
        total = 0
        do b = 1, nb
            nstep = max(3,ceiling((bands(2*b)-bands(2*b-1))/nyquist * real(gd*(numtaps+1),wp)))
            total = total+nstep+1
        end do
        if (total < nex .or. total > 2000000) return
        allocate(xx(total), target(total), ww(total), envelope(total), err(total), &
            ext(nex), new_ext(nex), nodes(nex), values(nex), lambdas(nex), coef(m))
        df = nyquist/real(gd*(numtaps+1),wp)
        idx = 0
        do b = 1, nb
            nstep = max(3,ceiling((bands(2*b)-bands(2*b-1))/nyquist * real(gd*(numtaps+1),wp)))
            do j = 0, nstep
                f = bands(2*b-1) + (bands(2*b)-bands(2*b-1))*real(j,wp)/real(nstep,wp)
                phase = pi*f/nyquist
                select case(ktype)
                case(1)
                    v = 1.0_wp
                case(2)
                    v = cos(phase/2.0_wp)
                case(3)
                    v = sin(phase)
                case(4)
                    v = sin(phase/2.0_wp)
                end select
                ! The antisymmetric response vanishes at DC; types II/III also
                ! vanish at Nyquist. Remove these constrained bins from the grid.
                if (abs(v) < 2.0e-12_wp) cycle
                idx = idx+1
                xx(idx) = cos(phase)
                envelope(idx) = v
                target(idx) = desired(b)
                scale = wband(b)
                if (kind == 'differentiator') then
                    target(idx) = desired(b)*f/rate
                    if (abs(desired(b)) > 1.0e-14_wp) scale = scale*rate/max(f,df/100.0_wp)
                end if
                target(idx) = target(idx)/v
                ww(idx) = scale*abs(v)
            end do
        end do
        nsize = idx
        if (nsize < nex) then
            if (present(status)) status = remez_status_no_convergence
            return
        end if
        do j = 1, nex
            ext(j) = 1 + nint(real(j-1,wp)*real(nsize-1,wp)/real(nex-1,wp))
        end do

        st = remez_status_no_convergence
        if (present(status)) status = st
        converged = .false.
        do loop = 1, mit
            call exchange_coefficients(xx(ext),target(ext),ww(ext),nodes,values,lambdas,delta,ok)
            if (.not. ok) exit
            do j = 1, nsize
                v = barycentric_value(xx(j),nodes,values,lambdas)
                err(j) = (v-target(j))*ww(j)
            end do
            if (.not. all(ieee_is_finite(err(:nsize)))) exit
            call choose_extrema(err(:nsize),nex,new_ext,ok)
            if (.not. ok) exit
            minerr = minval(abs(err(new_ext)))
            maxerr = maxval(abs(err(new_ext)))
            max_all = maxval(abs(err(:nsize)))
            if ((all(ext == new_ext) .or. maxerr-minerr <= 1.0e-8_wp*maxerr) .and. &
                max_all <= 1.005_wp*max(abs(delta),maxerr)) then
                ext = new_ext
                converged = .true.
                exit
            end if
            ext = new_ext
        end do
        if (present(iterations)) iterations = min(loop,mit)
        if (.not. converged) then
            if (present(status)) status = st
            return
        end if
        call exchange_coefficients(xx(ext),target(ext),ww(ext),nodes,values,lambdas,delta,ok)
        if (.not. ok) return
        call extract_cosine_coefficients(nodes,values,lambdas,coef)
        call reconstruct_impulse(coef,ktype,taps)
        if (.not. all(ieee_is_finite(taps))) then
            deallocate(taps)
            allocate(taps(0))
            return
        end if
        ! Detect numerically ill-conditioned polynomial reconstruction.
        ! In those cases an apparently converged exchange can produce taps
        ! whose *actual* response does not satisfy the minimax certificate.
        max_all = 0.0_wp
        do j = 1, nsize
            v = realized_amplitude(taps, ktype, acos(xx(j)))
            max_all = max(max_all, abs((v-target(j)*envelope(j))*ww(j)/envelope(j)))
        end do
        ! Allow only a rounding-scale floor, not a fixed absolute 1e-8
        ! threshold that can hide a poor reconstruction of a tiny ripple.
        scale = max(1.0_wp,real(size(taps),wp)*real(maxval(abs(taps)),wp))
        scale = scale*maxval(ww(:nsize)/abs(envelope(:nsize)))
        if (.not. ieee_is_finite(max_all) .or. &
            max_all > max(1.05_wp*abs(delta),40.0_wp*real(epsilon(1.0_dp),wp)*scale)) then
            deallocate(taps)
            allocate(taps(0))
            return
        end if
        if (present(ripple)) ripple = max_all
        if (present(status)) status = remez_status_success
    end function remez

    subroutine exchange_coefficients(x,d,w,nodes,values,lambdas,delta,ok)
        real(wp), intent(in) :: x(:) !! extremal cosines, strictly decreasing
        real(wp), intent(in) :: d(:) !! transformed desired response
        real(wp), intent(in) :: w(:) !! transformed positive weights
        real(wp), intent(out) :: nodes(:) !! interpolation nodes
        real(wp), intent(out) :: values(:) !! interpolation values
        real(wp), intent(out) :: lambdas(:) !! normalized barycentric weights
        real(wp), intent(out) :: delta !! signed equiripple level
        logical, intent(out) :: ok !! valid exchange state
        real(wp) :: logs(size(x)), largest, numerator, denominator, dx
        integer :: i, j, n

        n = size(x)
        nodes = x
        ok = .false.
        if (any(w <= 0.0_wp)) return
        do i = 1, n
            logs(i) = 0.0_wp
            do j = 1, n
                if (j == i) cycle
                dx = abs(x(i)-x(j))
                if (dx <= 0.0_wp) return
                logs(i) = logs(i)-log(dx)
            end do
        end do
        largest = maxval(logs)
        do i = 1, n
            lambdas(i) = exp(logs(i)-largest)
            if (mod(i,2) == 0) lambdas(i) = -lambdas(i)
        end do
        numerator = sum(lambdas*d)
        denominator = sum(abs(lambdas)/w)
        if (denominator <= 0.0_wp .or. .not. ieee_is_finite(denominator)) return
        delta = numerator/denominator
        do i = 1, n
            values(i) = d(i)-sign(1.0_wp,lambdas(i))*delta/w(i)
        end do
        ok = all(ieee_is_finite(values))
    end subroutine exchange_coefficients

    real(wp) function barycentric_value(x,nodes,values,lambdas) result(y)
        real(wp), intent(in) :: x !! evaluation point
        real(wp), intent(in) :: nodes(:) !! interpolation nodes
        real(wp), intent(in) :: values(:) !! polynomial values at nodes
        real(wp), intent(in) :: lambdas(:) !! barycentric weights
        real(wp) :: a, numerator, denominator
        integer :: k

        numerator = 0.0_wp
        denominator = 0.0_wp
        do k = 1, size(nodes)
            if (x == nodes(k)) then
                y = values(k)
                return
            end if
            a = lambdas(k)/(x-nodes(k))
            numerator = numerator+a*values(k)
            denominator = denominator+a
        end do
        y = numerator/denominator
    end function barycentric_value

    subroutine choose_extrema(error, ne, result, ok)
        real(wp), intent(in) :: error(:) !! current weighted approximation error
        integer, intent(in) :: ne !! number of required alternation extrema
        integer, intent(out) :: result(:) !! selected grid indices
        logical, intent(out) :: ok !! whether selection succeeded
        integer :: candidate(size(error)), nc, i, j, k, best_i, s, last_sign
        integer :: even_i, odd_i
        integer, allocatable :: prev(:,:)
        real(wp), allocatable :: score(:,:), logamp(:)
        real(wp) :: best, even_best, odd_best, val

        ok = .false.
        nc = 0
        last_sign = 0
        do i = 1, size(error)
            if (error(i) == 0.0_wp) cycle
            s = 1
            if (error(i) < 0.0_wp) s = -1
            if (s /= last_sign) then
                nc = nc+1
                candidate(nc) = i
                last_sign = s
            else if (abs(error(i)) >= abs(error(candidate(nc)))) then
                candidate(nc) = i
            end if
        end do
        if (nc < ne) return
        allocate(score(nc,ne),prev(nc,ne),logamp(nc))
        score = -huge(1.0_wp)
        prev = 0
        do i = 1, nc
            logamp(i) = log(max(abs(error(candidate(i))),tiny(1.0_wp)))
            score(i,1) = logamp(i)
        end do
        ! Candidate signs alternate, so the previous candidate must have
        ! opposite index parity. Prefix maxima reduce the selection from
        ! O(nc**2*ne) to O(nc*ne), preserving the original tie ordering.
        do k = 2, ne
            even_best = -huge(1.0_wp)
            odd_best = -huge(1.0_wp)
            even_i = 0
            odd_i = 0
            do i = k, nc
                j = i-1
                val = score(j,k-1)
                if (mod(j,2) == 0) then
                    if (val > even_best) then
                        even_best = val
                        even_i = j
                    end if
                else
                    if (val > odd_best) then
                        odd_best = val
                        odd_i = j
                    end if
                end if
                if (mod(i,2) == 0) then
                    if (odd_i == 0) cycle
                    score(i,k) = odd_best+logamp(i)
                    prev(i,k) = odd_i
                else
                    if (even_i == 0) cycle
                    score(i,k) = even_best+logamp(i)
                    prev(i,k) = even_i
                end if
            end do
        end do
        best = -huge(1.0_wp)
        best_i = 0
        do i = ne, nc
            if (score(i,ne) > best) then
                best = score(i,ne)
                best_i = i
            end if
        end do
        if (best_i == 0) return
        i = best_i
        do k = ne, 1, -1
            result(k) = candidate(i)
            i = prev(i,k)
        end do
        ok = .true.
    end subroutine choose_extrema

    subroutine extract_cosine_coefficients(nodes,values,lambdas,c)
        real(wp), intent(in) :: nodes(:) !! extremal cosines
        real(wp), intent(in) :: values(:) !! polynomial values at nodes
        real(wp), intent(in) :: lambdas(:) !! barycentric weights
        real(wp), intent(out) :: c(:) !! Chebyshev polynomial coefficients
        real(wp) :: ang, v
        integer :: i, k, m

        m = size(c)
        c = 0.0_wp
        do i = 1, m
            ang = pi*(real(i,wp)-0.5_wp)/real(m,wp)
            v = barycentric_value(cos(ang),nodes,values,lambdas)
            c(1) = c(1)+v
            do k = 2, m
                c(k) = c(k)+2.0_wp*v*cos(real(k-1,wp)*ang)
            end do
        end do
        c = c/real(m,wp)
    end subroutine extract_cosine_coefficients

    real(wp) function realized_amplitude(taps,ktype,omega) result(amp)
        real(dp), intent(in) :: taps(:) !! real FIR impulse response
        integer, intent(in) :: ktype !! linear-phase filter type
        real(wp), intent(in) :: omega !! angular frequency
        integer :: i, m

        amp = 0.0_wp
        select case (ktype)
        case (1)
            m = (size(taps)+1)/2
            amp = taps(m)
            do i = 1, m-1
                amp = amp+2.0_wp*taps(m-i)*cos(real(i,wp)*omega)
            end do
        case (2)
            m = size(taps)/2
            do i = 1, m
                amp = amp+2.0_wp*taps(m+1-i)*cos((real(i,wp)-0.5_wp)*omega)
            end do
        case (3)
            m = (size(taps)-1)/2
            do i = 1, m
                amp = amp+2.0_wp*taps(m+1-i)*sin(real(i,wp)*omega)
            end do
        case (4)
            m = size(taps)/2
            do i = 1, m
                amp = amp+2.0_wp*taps(m+1-i)*sin((real(i,wp)-0.5_wp)*omega)
            end do
        end select
    end function realized_amplitude

    subroutine reconstruct_impulse(c,ktype,taps)
        real(wp), intent(in) :: c(:) !! cosine polynomial coefficients, index 1 is constant
        integer, intent(in) :: ktype !! linear-phase FIR type I, II, III or IV
        real(dp), allocatable, intent(out) :: taps(:) !! symmetric/antisymmetric FIR taps
        real(wp), allocatable :: b(:)
        integer :: m, i, n

        m = size(c)
        select case(ktype)
        case(1)
            n = 2*m-1
            allocate(taps(n))
            taps(m) = c(1)
            do i = 1, m-1
                taps(m-i) = 0.5_wp*c(i+1)
                taps(m+i) = taps(m-i)
            end do
        case(2)
            n = 2*m
            allocate(taps(n),b(m))
            b = 0.0_wp
            b(1) = c(1)
            do i = 2, m
                b(i-1) = b(i-1)+0.5_wp*c(i)
                b(i) = b(i)+0.5_wp*c(i)
            end do
            do i = 1, m
                taps(m+1-i) = 0.5_wp*b(i)
                taps(m+i) = taps(m+1-i)
            end do
        case(3)
            n = 2*m+1
            allocate(taps(n),b(m))
            b = 0.0_wp
            b(1) = c(1)
            do i = 2, m
                b(i) = b(i)+0.5_wp*c(i)
                if (i > 2) b(i-2) = b(i-2)-0.5_wp*c(i)
            end do
            taps(m+1) = 0.0_wp
            do i = 1, m
                taps(m+1-i) = 0.5_wp*b(i)
                taps(m+1+i) = -taps(m+1-i)
            end do
        case(4)
            n = 2*m
            allocate(taps(n),b(m))
            b = 0.0_wp
            b(1) = c(1)
            do i = 2, m
                b(i) = b(i)+0.5_wp*c(i)
                b(i-1) = b(i-1)-0.5_wp*c(i)
            end do
            do i = 1, m
                taps(m+1-i) = 0.5_wp*b(i)
                taps(m+i) = -taps(m+1-i)
            end do
        end select
    end subroutine reconstruct_impulse
end module scifort_signal_remez
