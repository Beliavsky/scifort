! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Gustafsson initial conditions for forward-backward digital filtering.
! The zero-state forward/backward and backward/forward responses are
! reconciled by a small complex least-squares problem.  Real and complex
! signals use the same kernel; no BLAS or LAPACK dependency is required.
module scifort_signal_filtfilt_gust
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal_lfilter, only : lfilter, lfilter_status_success
    use scifort_signal_sos, only : sosfilt, sos_status_success
    implicit none
    private

    integer, parameter, public :: gust_status_success = 0
    integer, parameter, public :: gust_status_invalid_input = 1
    integer, parameter, public :: gust_status_numerical_failure = 2

    interface filtfilt_gust
        module procedure gust_real, gust_complex, gust_real_coeff_complex_signal
        module procedure gust_complex_coeff_real_signal
    end interface filtfilt_gust
    interface sosfiltfilt_gust
        module procedure sos_gust_real, sos_gust_complex
        module procedure sos_gust_real_coeff_complex_signal, sos_gust_complex_coeff_real_signal
    end interface sosfiltfilt_gust
    public :: filtfilt_gust, sosfiltfilt_gust

contains

    function gust_real(b, a, x, irlen, status) result(y)
        real(dp), intent(in) :: b(:) !! numerator, in increasing lag order
        real(dp), intent(in) :: a(:) !! denominator, nonzero first coefficient
        real(dp), intent(in) :: x(:) !! real-valued input signal
        integer, intent(in), optional :: irlen !! positive impulse-response truncation length
        integer, intent(out), optional :: status !! zero on success
        real(dp), allocatable :: y(:)
        complex(dp), allocatable :: yc(:)
        integer :: code

        yc = gust_kernel(cmplx(b, 0.0_dp, dp), cmplx(a, 0.0_dp, dp), &
            cmplx(x, 0.0_dp, dp), irlen, code)
        y = real(yc, dp)
        if (present(status)) status = code
    end function gust_real

    function gust_complex(b, a, x, irlen, status) result(y)
        complex(dp), intent(in) :: b(:) !! numerator, in increasing lag order
        complex(dp), intent(in) :: a(:) !! denominator, nonzero first coefficient
        complex(dp), intent(in) :: x(:) !! complex-valued input signal
        integer, intent(in), optional :: irlen !! positive impulse-response truncation length
        integer, intent(out), optional :: status !! zero on success
        complex(dp), allocatable :: y(:)
        integer :: code

        y = gust_kernel(b, a, x, irlen, code)
        if (present(status)) status = code
    end function gust_complex

    function gust_real_coeff_complex_signal(b, a, x, irlen, status) result(y)
        real(dp), intent(in) :: b(:) !! real numerator
        real(dp), intent(in) :: a(:) !! real denominator
        complex(dp), intent(in) :: x(:) !! complex input
        integer, intent(in), optional :: irlen !! impulse-response truncation length
        integer, intent(out), optional :: status !! zero on success
        complex(dp), allocatable :: y(:)
        integer :: code

        y = gust_kernel(cmplx(b, 0.0_dp, dp), cmplx(a, 0.0_dp, dp), x, irlen, code)
        if (present(status)) status = code
    end function gust_real_coeff_complex_signal

    function gust_complex_coeff_real_signal(b, a, x, irlen, status) result(y)
        complex(dp), intent(in) :: b(:) !! complex numerator
        complex(dp), intent(in) :: a(:) !! complex denominator
        real(dp), intent(in) :: x(:) !! real-valued input
        integer, intent(in), optional :: irlen !! impulse-response truncation length
        integer, intent(out), optional :: status !! zero on success
        complex(dp), allocatable :: y(:)
        integer :: code

        y = gust_kernel(b, a, cmplx(x, 0.0_dp, dp), irlen, code)
        if (present(status)) status = code
    end function gust_complex_coeff_real_signal

    function gust_kernel(b, a, x, irlen, status) result(y)
        complex(dp), intent(in) :: b(:) !! numerator coefficients
        complex(dp), intent(in) :: a(:) !! denominator coefficients
        complex(dp), intent(in) :: x(:) !! input samples
        integer, intent(in), optional :: irlen !! positive effective impulse-response length
        integer, intent(out) :: status !! status code
        complex(dp), allocatable :: y(:)
        complex(dp), allocatable :: obs(:,:), s(:,:), matrix(:,:), correction(:,:), &
            delta(:), ic(:), zeros(:), zi(:), v(:), forward(:), backward(:), bf(:)
        real(dp) :: nanval
        integer :: n, order, m, k, code
        logical :: truncated, solved

        nanval = ieee_value(0.0_dp, ieee_quiet_nan)
        n = size(x)
        allocate(y(n))
        y = cmplx(nanval, nanval, dp)
        status = gust_status_invalid_input
        if (size(b) < 1 .or. size(a) < 1 .or. n < 1) return
        if (abs(a(1)) <= tiny(1.0_dp)) return
        if (any(.not. ieee_is_finite(real(b, dp))) .or. any(.not. ieee_is_finite(aimag(b)))) return
        if (any(.not. ieee_is_finite(real(a, dp))) .or. any(.not. ieee_is_finite(aimag(a)))) return
        if (any(.not. ieee_is_finite(real(x, dp))) .or. any(.not. ieee_is_finite(aimag(x)))) return
        if (present(irlen)) then
            if (irlen < 1) return
        end if

        order = max(size(b), size(a)) - 1
        if (order == 0) then
            y = (b(1) / a(1))**2 * x
            status = gust_status_success
            return
        end if
        m = n
        if (present(irlen)) then
            if (irlen < (n+1)/2) m = irlen
        end if
        truncated = m < n

        ! Columns of the observability matrix: output from a single unit
        ! state with zero signal. S filters the reversed state response.
        allocate(obs(m,order), s(m,order), zeros(m), zi(order))
        zeros = (0.0_dp, 0.0_dp)
        do k = 1, order
            zi = (0.0_dp, 0.0_dp)
            zi(k) = (1.0_dp, 0.0_dp)
            obs(:,k) = lfilter(b, a, zeros, zi=zi, status=code)
            if (code /= lfilter_status_success) return
            s(:,k) = lfilter(b, a, obs(m:1:-1,k), status=code)
            if (code /= lfilter_status_success) return
        end do

        ! The discrepancy of the two zero-state filtering orders supplies
        ! the least-squares right-hand side.
        forward = lfilter(b, a, x, status=code)
        if (code /= lfilter_status_success) return
        y = lfilter(b, a, forward(n:1:-1), status=code)
        if (code /= lfilter_status_success) return
        y = y(n:1:-1)
        backward = lfilter(b, a, x(n:1:-1), status=code)
        if (code /= lfilter_status_success) return
        bf = lfilter(b, a, backward(n:1:-1), status=code)
        if (code /= lfilter_status_success) return
        if (truncated) then
            allocate(matrix(2*m,2*order), correction(2*m,2*order), delta(2*m))
            matrix = (0.0_dp, 0.0_dp)
            correction = (0.0_dp, 0.0_dp)
            matrix(1:m,1:order) = s(m:1:-1,:) - obs
            matrix(m+1:2*m,order+1:2*order) = obs(m:1:-1,:) - s
            correction(1:m,1:order) = s(m:1:-1,:)
            correction(m+1:2*m,order+1:2*order) = obs(m:1:-1,:)
            delta(1:m) = bf(1:m) - y(1:m)
            delta(m+1:2*m) = bf(n-m+1:n) - y(n-m+1:n)
        else
            allocate(matrix(n,2*order), correction(n,2*order), delta(n))
            matrix(:,1:order) = s(m:1:-1,:) - obs
            matrix(:,order+1:2*order) = obs(m:1:-1,:) - s
            correction(:,1:order) = s(m:1:-1,:)
            correction(:,order+1:2*order) = obs(m:1:-1,:)
            delta = bf - y
        end if

        solved = .false.
        if (size(matrix,1) <= size(matrix,2)) then
            call householder_minimum_norm(matrix,delta,ic,solved)
        end if
        if (.not. solved) then
            if (size(matrix,1) <= 1200 .and. size(matrix,2) <= 32) then
                call jacobi_least_squares(matrix,delta,ic,solved)
                if (.not. solved) call pivoted_least_squares(matrix,delta,ic,solved)
            else
                call pivoted_least_squares(matrix,delta,ic,solved)
            end if
        end if
        if (.not. solved) then
            y = cmplx(nanval, nanval, dp)
            status = gust_status_numerical_failure
            return
        end if
        v = matmul(correction,ic)
        if (truncated) then
            y(1:m) = y(1:m) + v(1:m)
            y(n-m+1:n) = y(n-m+1:n) + v(m+1:2*m)
        else
            y = y + v
        end if
        if (.not. all(ieee_is_finite(real(y, dp))) .or. &
            .not. all(ieee_is_finite(aimag(y)))) then
            status = gust_status_numerical_failure
            return
        end if
        status = gust_status_success
    end function gust_kernel


    ! SOS variant of the Gustafsson construction.  Unlike conversion via
    ! sos2tf, every state-response column and filtering pass uses a cascade
    ! of second-order sections, avoiding high-order polynomial coefficients.
    function sos_gust_real(sos, x, irlen, status) result(y)
        real(dp), intent(in) :: sos(:,:) !! SOS rows (nsections,6), a0=1
        real(dp), intent(in) :: x(:) !! input samples
        integer, intent(in), optional :: irlen !! positive impulse response cutoff
        integer, intent(out), optional :: status !! zero on success
        real(dp), allocatable :: y(:)
        complex(dp), allocatable :: yc(:)
        integer :: code

        yc = sos_gust_kernel(cmplx(sos,0.0_dp,dp),cmplx(x,0.0_dp,dp),irlen,code)
        y = real(yc,dp)
        if (present(status)) status = code
    end function sos_gust_real

    function sos_gust_complex(sos, x, irlen, status) result(y)
        complex(dp), intent(in) :: sos(:,:) !! SOS rows (nsections,6), a0=1
        complex(dp), intent(in) :: x(:) !! input samples
        integer, intent(in), optional :: irlen !! positive impulse response cutoff
        integer, intent(out), optional :: status !! zero on success
        complex(dp), allocatable :: y(:)
        integer :: code

        y = sos_gust_kernel(sos,x,irlen,code)
        if (present(status)) status = code
    end function sos_gust_complex

    function sos_gust_real_coeff_complex_signal(sos,x,irlen,status) result(y)
        real(dp), intent(in) :: sos(:,:) !! real SOS rows
        complex(dp), intent(in) :: x(:) !! complex input signal
        integer, intent(in), optional :: irlen !! positive impulse response cutoff
        integer, intent(out), optional :: status !! zero on success
        complex(dp), allocatable :: y(:)
        integer :: code

        y = sos_gust_kernel(cmplx(sos,0.0_dp,dp),x,irlen,code)
        if (present(status)) status = code
    end function sos_gust_real_coeff_complex_signal

    function sos_gust_complex_coeff_real_signal(sos,x,irlen,status) result(y)
        complex(dp), intent(in) :: sos(:,:) !! complex SOS rows
        real(dp), intent(in) :: x(:) !! real input signal
        integer, intent(in), optional :: irlen !! positive impulse response cutoff
        integer, intent(out), optional :: status !! zero on success
        complex(dp), allocatable :: y(:)
        integer :: code

        y = sos_gust_kernel(sos,cmplx(x,0.0_dp,dp),irlen,code)
        if (present(status)) status = code
    end function sos_gust_complex_coeff_real_signal

    function sos_gust_kernel(sos,x,irlen,status) result(y)
        complex(dp), intent(in) :: sos(:,:) !! SOS rows, with a0=1
        complex(dp), intent(in) :: x(:) !! input samples
        integer, intent(in), optional :: irlen !! positive impulse-response cutoff
        integer, intent(out) :: status !! success or error
        complex(dp), allocatable :: y(:)
        complex(dp), allocatable :: obs(:,:), s(:,:), matrix(:,:), correction(:,:), &
            delta(:), ic(:), zeros(:), states(:,:), forward(:), backward(:), bf(:), v(:)
        real(dp) :: nanval
        integer :: n, m, nsection, order, k, code
        logical :: truncated, solved

        nanval = ieee_value(0.0_dp,ieee_quiet_nan)
        n = size(x)
        allocate(y(n))
        y = cmplx(nanval,nanval,dp)
        status = gust_status_invalid_input
        if (n < 1 .or. size(sos,1) < 1 .or. size(sos,2) /= 6) return
        if (any(sos(:,4) /= (1.0_dp,0.0_dp))) return
        if (.not. all(ieee_is_finite(real(sos,dp)))) return
        if (.not. all(ieee_is_finite(aimag(sos)))) return
        if (.not. all(ieee_is_finite(real(x,dp)))) return
        if (.not. all(ieee_is_finite(aimag(x)))) return
        if (present(irlen)) then
            if (irlen < 1) return
        end if
        nsection = size(sos,1)
        if (nsection > huge(nsection)/2) return
        order = 2*nsection
        m = n
        if (present(irlen)) then
            if (irlen < (n+1)/2) m = irlen
        end if
        truncated = m < n

        allocate(obs(m,order),s(m,order),zeros(m),states(nsection,2))
        zeros = (0.0_dp,0.0_dp)
        do k = 1,order
            states = (0.0_dp,0.0_dp)
            states(mod(k-1,nsection)+1,(k-1)/nsection+1) = (1.0_dp,0.0_dp)
            obs(:,k) = sosfilt(sos,zeros,zi=states,status=code)
            if (code /= sos_status_success) return
            s(:,k) = sosfilt(sos,obs(m:1:-1,k),status=code)
            if (code /= sos_status_success) return
        end do
        forward = sosfilt(sos,x,status=code)
        if (code /= sos_status_success) return
        y = sosfilt(sos,forward(n:1:-1),status=code)
        if (code /= sos_status_success) return
        y = y(n:1:-1)
        backward = sosfilt(sos,x(n:1:-1),status=code)
        if (code /= sos_status_success) return
        bf = sosfilt(sos,backward(n:1:-1),status=code)
        if (code /= sos_status_success) return
        if (truncated) then
            allocate(matrix(2*m,2*order),correction(2*m,2*order),delta(2*m))
            matrix = (0.0_dp,0.0_dp)
            correction = (0.0_dp,0.0_dp)
            matrix(1:m,1:order) = s(m:1:-1,:) - obs
            matrix(m+1:2*m,order+1:2*order) = obs(m:1:-1,:) - s
            correction(1:m,1:order) = s(m:1:-1,:)
            correction(m+1:2*m,order+1:2*order) = obs(m:1:-1,:)
            delta(1:m) = bf(1:m) - y(1:m)
            delta(m+1:2*m) = bf(n-m+1:n) - y(n-m+1:n)
        else
            allocate(matrix(n,2*order),correction(n,2*order),delta(n))
            matrix(:,1:order) = s(m:1:-1,:) - obs
            matrix(:,order+1:2*order) = obs(m:1:-1,:) - s
            correction(:,1:order) = s(m:1:-1,:)
            correction(:,order+1:2*order) = obs(m:1:-1,:)
            delta = bf - y
        end if
        if (.not. all(ieee_is_finite(real(matrix,dp))) .or. &
            .not. all(ieee_is_finite(aimag(matrix))) .or. &
            .not. all(ieee_is_finite(real(delta,dp))) .or. &
            .not. all(ieee_is_finite(aimag(delta)))) then
            y = cmplx(nanval,nanval,dp)
            status = gust_status_numerical_failure
            return
        end if
        solved = .false.
        if (size(matrix,1) <= size(matrix,2)) then
            call householder_minimum_norm(matrix,delta,ic,solved)
        end if
        if (.not. solved) then
            if (size(matrix,1) <= 1200 .and. size(matrix,2) <= 32) then
                call jacobi_least_squares(matrix,delta,ic,solved)
                if (.not. solved) call pivoted_least_squares(matrix,delta,ic,solved)
            else
                call pivoted_least_squares(matrix,delta,ic,solved)
            end if
        end if
        if (.not. solved) then
            y = cmplx(nanval,nanval,dp)
            status = gust_status_numerical_failure
            return
        end if
        ! Ill-conditioned, short, high-order cascades can make the initial
        ! state system numerically rank-deficient. A pivoted solution with a
        ! large consistency residual must not be reported as a valid filter.
        if (.not. truncated .and. maxval(abs(matmul(matrix,ic)-delta)) > &
            1.0e-9_dp * max(maxval(abs(delta)),tiny(1.0_dp))) then
            y = cmplx(nanval,nanval,dp)
            status = gust_status_numerical_failure
            return
        end if
        v = matmul(correction,ic)
        if (truncated) then
            y(1:m) = y(1:m) + v(1:m)
            y(n-m+1:n) = y(n-m+1:n) + v(m+1:2*m)
        else
            y = y + v
        end if
        if (.not. all(ieee_is_finite(real(y,dp))) .or. &
            .not. all(ieee_is_finite(aimag(y)))) then
            y = cmplx(nanval,nanval,dp)
            status = gust_status_numerical_failure
            return
        end if
        status = gust_status_success
    end function sos_gust_kernel

    ! For an underdetermined system A*x=b with independent rows, compute
    ! the minimum-norm solution by a pivoted Householder QR factorization
    ! of A^H.  In contrast to normal equations this does not square the
    ! condition number, and Householder reflections retain orthogonality
    ! more accurately than modified Gram-Schmidt on ill-scaled SOS states.
    ! Rank-deficient systems fall back to the Jacobi SVD below.
    subroutine householder_minimum_norm(matrix, rhs, solution, ok)
        complex(dp), intent(in) :: matrix(:,:) !! original n-by-p design matrix
        complex(dp), intent(in) :: rhs(:) !! n observations
        complex(dp), allocatable, intent(out) :: solution(:) !! p minimum-norm states
        logical, intent(out) :: ok !! rows independent and finite solution
        complex(dp), allocatable :: hmat(:,:), reflect(:,:), v(:), tempcol(:), &
            z(:), reordered(:)
        real(dp), allocatable :: tau(:)
        integer, allocatable :: permutation(:)
        complex(dp) :: alpha, proj
        real(dp) :: diag_norm, pivot_norm, limit_norm
        integer :: n, p, j, k, pick, oldindex

        n = size(matrix,1)
        p = size(matrix,2)
        allocate(solution(p))
        solution = (0.0_dp,0.0_dp)
        ok = .false.
        if (n < 1 .or. n > p .or. size(rhs) /= n) return
        if (.not. all(ieee_is_finite(real(matrix,dp)))) return
        if (.not. all(ieee_is_finite(aimag(matrix)))) return
        if (.not. all(ieee_is_finite(real(rhs,dp)))) return
        if (.not. all(ieee_is_finite(aimag(rhs)))) return

        allocate(hmat(p,n), reflect(p,n), permutation(n), tau(n), &
            tempcol(p), z(n), reordered(n))
        hmat = conjg(transpose(matrix))
        reflect = (0.0_dp,0.0_dp)
        do j = 1,n
            permutation(j) = j
        end do
        limit_norm = real(max(n,p),dp)*epsilon(1.0_dp)* &
            maxval(abs(hmat))
        do j = 1,n
            ! Column pivoting protects near-dependent rows; swap the full
            ! columns to carry already formed triangular coefficients.
            pick = j
            pivot_norm = -1.0_dp
            do k = j,n
                diag_norm = sqrt(sum(abs(hmat(j:p,k))**2))
                if (diag_norm > pivot_norm) then
                    pick = k
                    pivot_norm = diag_norm
                end if
            end do
            if (.not. ieee_is_finite(pivot_norm)) return
            if (pivot_norm <= max(limit_norm,tiny(1.0_dp))) return
            if (pick /= j) then
                tempcol = hmat(:,j)
                hmat(:,j) = hmat(:,pick)
                hmat(:,pick) = tempcol
                oldindex = permutation(j)
                permutation(j) = permutation(pick)
                permutation(pick) = oldindex
            end if
            allocate(v(p-j+1))
            v = hmat(j:p,j)
            if (abs(v(1)) > 0.0_dp) then
                alpha = -v(1)/abs(v(1))*pivot_norm
            else
                alpha = cmplx(-pivot_norm,0.0_dp,dp)
            end if
            v(1) = v(1) - alpha
            tau(j) = 2.0_dp/real(dot_product(v,v),dp)
            if (.not. ieee_is_finite(tau(j))) return
            reflect(j:p,j) = v
            do k = j,n
                proj = tau(j)*dot_product(v,hmat(j:p,k))
                hmat(j:p,k) = hmat(j:p,k) - v*proj
            end do
            hmat(j,j) = alpha
            if (j < p) hmat(j+1:p,j) = (0.0_dp,0.0_dp)
            deallocate(v)
        end do

        do j = 1,n
            reordered(j) = rhs(permutation(j))
            z(j) = (reordered(j) - &
                sum(conjg(hmat(1:j-1,j))*z(1:j-1))) / conjg(hmat(j,j))
        end do
        solution(1:n) = z
        ! Q = H_1 * ... * H_n; apply its reflectors in reverse order.
        do j = n,1,-1
            proj = tau(j)*dot_product(reflect(j:p,j),solution(j:p))
            solution(j:p) = solution(j:p) - reflect(j:p,j)*proj
        end do
        ok = all(ieee_is_finite(real(solution,dp))) .and. &
            all(ieee_is_finite(aimag(solution)))
    end subroutine householder_minimum_norm

    ! One-sided complex Jacobi SVD for small, potentially rank-deficient
    ! least-squares problems. Orthogonalization is performed on the actual
    ! design columns, avoiding the squared condition of normal equations.
    ! The resulting Moore-Penrose solution is the minimum-norm one.
    subroutine jacobi_least_squares(matrix, rhs, solution, ok)
        complex(dp), intent(in) :: matrix(:,:) !! least-squares design
        complex(dp), intent(in) :: rhs(:) !! right-hand side
        complex(dp), allocatable, intent(out) :: solution(:) !! minimum-norm states
        logical, intent(out) :: ok !! true when converged with finite results
        complex(dp), allocatable :: cols(:,:), vectors(:,:), first(:), second(:), &
            vf(:), vs(:)
        complex(dp) :: gamma, phase
        real(dp) :: alpha, beta, cross, tau, tang, cosine, sine, scale, thresh, sigma2
        integer :: n, p, j, k, sweep
        logical :: changed

        n = size(matrix,1)
        p = size(matrix,2)
        allocate(solution(p),cols(n,p),vectors(p,p),first(n),second(n),vf(p),vs(p))
        solution = (0.0_dp,0.0_dp)
        ok = .false.
        if (n /= size(rhs)) return
        cols = matrix
        vectors = (0.0_dp,0.0_dp)
        do j=1,p
            vectors(j,j) = (1.0_dp,0.0_dp)
        end do
        do sweep=1,80
            changed=.false.
            do j=1,p-1
                alpha=sum(abs(cols(:,j))**2)
                if (alpha <= tiny(1.0_dp)) cycle
                do k=j+1,p
                    beta=sum(abs(cols(:,k))**2)
                    if (beta <= tiny(1.0_dp)) cycle
                    gamma=dot_product(cols(:,j),cols(:,k))
                    cross=abs(gamma)
                    if (cross <= 2.0_dp*epsilon(1.0_dp)*sqrt(alpha)*sqrt(beta)) cycle
                    changed=.true.
                    phase=gamma/cross
                    tau=(beta-alpha)/(2.0_dp*cross)
                    tang=sign(1.0_dp,tau)/(abs(tau)+sqrt(1.0_dp+tau*tau))
                    cosine=1.0_dp/sqrt(1.0_dp+tang*tang)
                    sine=tang*cosine
                    first=cols(:,j)
                    second=cols(:,k)
                    cols(:,j)=cosine*first-sine*conjg(phase)*second
                    cols(:,k)=sine*phase*first+cosine*second
                    vf=vectors(:,j)
                    vs=vectors(:,k)
                    vectors(:,j)=cosine*vf-sine*conjg(phase)*vs
                    vectors(:,k)=sine*phase*vf+cosine*vs
                    alpha=sum(abs(cols(:,j))**2)
                end do
            end do
            if (.not. changed) exit
        end do
        if (changed) return
        scale=0.0_dp
        do j=1,p
            scale=max(scale,sqrt(sum(abs(cols(:,j))**2)))
        end do
        thresh=real(max(n,p),dp)*epsilon(1.0_dp)*scale
        do j=1,p
            sigma2=sum(abs(cols(:,j))**2)
            if (sqrt(sigma2) > thresh .and. sigma2 > 0.0_dp) then
                solution=solution+vectors(:,j)*(dot_product(cols(:,j),rhs)/sigma2)
            end if
        end do
        ok=all(ieee_is_finite(real(solution,dp))) .and. &
            all(ieee_is_finite(aimag(solution)))
    end subroutine jacobi_least_squares

    ! Reorthogonalized column-pivoted QR solve. The filter orders of interest
    ! are small compared with a signal; rank-deficient columns are omitted.
    subroutine pivoted_least_squares(matrix, rhs, solution, ok)
        complex(dp), intent(in) :: matrix(:,:) !! least-squares matrix, nobs by nstates
        complex(dp), intent(in) :: rhs(:) !! discrepancy between filtering orders
        complex(dp), allocatable, intent(out) :: solution(:) !! direct-form-II optimal states
        logical, intent(out) :: ok !! true if a finite solution was computed
        complex(dp), allocatable :: residual(:,:), q(:,:), r(:,:), v(:), factors(:), &
            bmat(:,:), hmat(:,:), rhs_spd(:), tail(:), w(:)
        real(dp), allocatable :: norms(:)
        real(dp) :: norm0, pivot_norm, cutoff, temp_norm
        complex(dp) :: temp
        integer, allocatable :: permutation(:)
        integer :: nr, nc, rank, i, j, k, pick, item, extra
        logical :: solved

        nr = size(matrix,1)
        nc = size(matrix,2)
        allocate(solution(nc), permutation(nc), norms(nc))
        solution = (0.0_dp, 0.0_dp)
        ok = .false.
        if (size(rhs) /= nr) return
        allocate(residual(nr,nc), q(nr,min(nr,nc)), &
            r(min(nr,nc),nc), factors(min(nr,nc)), v(nr))
        residual = matrix
        r = (0.0_dp, 0.0_dp)
        do j = 1, nc
            permutation(j) = j
            norms(j) = sqrt(sum(abs(residual(:,j))**2))
        end do
        norm0 = maxval(norms)
        cutoff = max(10.0_dp*epsilon(1.0_dp)*real(max(nr,nc),dp)*norm0, tiny(1.0_dp))
        rank = 0
        do j = 1, min(nr,nc)
            pick = maxloc(norms(j:nc),dim=1) + j - 1
            pivot_norm = norms(pick)
            if (pivot_norm < cutoff) exit
            if (pick /= j) then
                v = residual(:,j)
                residual(:,j) = residual(:,pick)
                residual(:,pick) = v
                temp_norm = norms(j)
                norms(j) = norms(pick)
                norms(pick) = temp_norm
                item = permutation(j)
                permutation(j) = permutation(pick)
                permutation(pick) = item
                if (j > 1) then
                    do i = 1, j-1
                        temp = r(i,j)
                        r(i,j) = r(i,pick)
                        r(i,pick) = temp
                    end do
                end if
            end if
            rank = j
            r(j,j) = cmplx(sqrt(sum(abs(residual(:,j))**2)), 0.0_dp, dp)
            if (abs(r(j,j)) < cutoff) then
                rank = j-1
                exit
            end if
            q(:,j) = residual(:,j) / r(j,j)
            do k = j+1, nc
                temp = dot_product(q(:,j),residual(:,k))
                r(j,k) = temp
                residual(:,k) = residual(:,k) - q(:,j) * temp
                ! A second orthogonalization pass reduces lost digits.
                temp = dot_product(q(:,j),residual(:,k))
                r(j,k) = r(j,k) + temp
                residual(:,k) = residual(:,k) - q(:,j) * temp
                norms(k) = sqrt(sum(abs(residual(:,k))**2))
            end do
        end do
        factors = (0.0_dp, 0.0_dp)
        if (rank > 0) then
            do i = 1, rank
                factors(i) = dot_product(q(:,i),rhs)
            end do
            do i = rank, 1, -1
                factors(i) = (factors(i) - sum(r(i,i+1:rank) * factors(i+1:rank))) / r(i,i)
            end do
            extra = nc-rank
            if (extra > 0) then
                ! QR yields an upper-trapezoidal full-row-rank system
                ! [R11 R12] z = q^H rhs. Its minimum-norm solution has
                ! z = [w - B*t; t], where B=R11^-1 R12 and
                ! (I+B^H B)*t = B^H*w. This matters for short signals.
                allocate(bmat(rank,extra), w(rank), hmat(extra,extra), &
                    rhs_spd(extra), tail(extra))
                w = factors(1:rank)
                do k = 1, extra
                    bmat(:,k) = r(1:rank,rank+k)
                    do i = rank,1,-1
                        bmat(i,k) = (bmat(i,k) - &
                            sum(r(i,i+1:rank)*bmat(i+1:rank,k))) / r(i,i)
                    end do
                end do
                hmat = matmul(conjg(transpose(bmat)), bmat)
                do k = 1, extra
                    hmat(k,k) = hmat(k,k) + (1.0_dp,0.0_dp)
                end do
                rhs_spd = matmul(conjg(transpose(bmat)),w)
                call solve_hermitian_spd(hmat,rhs_spd,tail,solved)
                if (.not. solved) return
                factors(1:rank) = w - matmul(bmat,tail)
                do i = 1, extra
                    solution(permutation(rank+i)) = tail(i)
                end do
            end if
            do i = 1, rank
                solution(permutation(i)) = factors(i)
            end do
        end if
        ok = all(ieee_is_finite(real(solution, dp))) .and. &
            all(ieee_is_finite(aimag(solution)))
    end subroutine pivoted_least_squares

    ! Cholesky solve for a Hermitian positive definite matrix.
    subroutine solve_hermitian_spd(hmat, rhs, sol, ok)
        complex(dp), intent(in) :: hmat(:,:) !! Hermitian matrix
        complex(dp), intent(in) :: rhs(:) !! right-hand side
        complex(dp), intent(out) :: sol(:) !! solution
        logical, intent(out) :: ok !! false if matrix cannot be factored
        complex(dp), allocatable :: lower(:,:), intermediate(:)
        complex(dp) :: diagonal
        integer :: n, i, j

        n = size(rhs)
        ok = .false.
        allocate(lower(n,n), intermediate(n))
        lower = (0.0_dp, 0.0_dp)
        do j = 1, n
            diagonal = hmat(j,j) - sum(lower(j,1:j-1)*conjg(lower(j,1:j-1)))
            if (.not. ieee_is_finite(real(diagonal, dp))) return
            if (real(diagonal,dp) <= 0.0_dp) return
            lower(j,j) = cmplx(sqrt(real(diagonal,dp)),0.0_dp,dp)
            do i = j+1,n
                lower(i,j) = (hmat(i,j) - &
                    dot_product(lower(j,1:j-1),lower(i,1:j-1))) / lower(j,j)
            end do
        end do
        do i = 1,n
            intermediate(i) = (rhs(i) - &
                sum(lower(i,1:i-1)*intermediate(1:i-1))) / lower(i,i)
        end do
        do i = n,1,-1
            sol(i) = (intermediate(i) - &
                sum(conjg(lower(i+1:n,i))*sol(i+1:n))) / conjg(lower(i,i))
        end do
        ok = all(ieee_is_finite(real(sol,dp))) .and. all(ieee_is_finite(aimag(sol)))
    end subroutine solve_hermitian_spd

end module scifort_signal_filtfilt_gust
