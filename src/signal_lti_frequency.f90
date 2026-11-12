! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! SciPy-style continuous/discrete LTI frequency response and Bode diagrams.
! Transfer-function polynomials use descending powers of s or z.  Discrete
! frequency inputs are angular frequencies in radians per sample; dbode's
! returned frequency axis is in radians per second (scaled by sampling dt).
module scifort_signal_lti_frequency
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_negative_inf, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal_design, only : frequency_response, freqz, freqz_zpk, &
        design_status_success, design_status_invalid_input, design_status_singular
    use scifort_signal_analog, only : freqs, freqs_zpk
    implicit none
    private

    real(dp), parameter :: pi = acos(-1.0_dp)

    type, public :: bode_result
        real(dp), allocatable :: frequencies(:) !! angular frequency, rad/s
        real(dp), allocatable :: magnitude(:) !! magnitude response, dB
        real(dp), allocatable :: phase(:) !! unwrapped phase, degrees
        integer :: status = design_status_invalid_input !! evaluation status
    end type bode_result

    ! MIMO results are indexed as (output, input, frequency).
    ! Separate public names preserve SISO result types and existing calls.
    type, public :: mimo_frequency_response
        real(dp), allocatable :: frequencies(:) !! analog rad/s or digital rad/sample
        complex(dp), allocatable :: response(:,:,:) !! (output, input, frequency)
        integer :: status = design_status_invalid_input !! evaluation status
    end type mimo_frequency_response

    type, public :: mimo_bode_result
        real(dp), allocatable :: frequencies(:) !! angular frequency, rad/s
        real(dp), allocatable :: magnitude(:,:,:) !! decibels (output, input, frequency)
        real(dp), allocatable :: phase(:,:,:) !! unwrapped degrees (output, input, frequency)
        integer :: status = design_status_invalid_input !! evaluation status
    end type mimo_bode_result

    ! A plan retains a unitary Hessenberg or Schur reduction for repeated sweeps.
    ! Matrix components must not be modified after preparation.
    type, public :: mimo_frequency_plan
        complex(dp), allocatable :: h(:,:) !! reduced state matrix
        complex(dp), allocatable :: br(:,:) !! transformed input matrix
        complex(dp), allocatable :: cr(:,:) !! transformed output matrix
        complex(dp), allocatable :: d(:,:) !! direct feedthrough matrix
        real(dp) :: analog_log_range(2) = [-1.0_dp,1.0_dp] !! default analog grid range
        integer :: status = design_status_invalid_input !! preparation status
        character(len=11) :: method = 'hessenberg' !! actual prepared solver
    end type mimo_frequency_plan

    interface prepare_mimo_frequency
        module procedure prepare_mimo_frequency_real, prepare_mimo_frequency_complex
    end interface
    interface freqresp_mimo
        module procedure freqresp_mimo_real, freqresp_mimo_complex, freqresp_mimo_plan
    end interface
    interface dfreqresp_mimo
        module procedure dfreqresp_mimo_real, dfreqresp_mimo_complex, dfreqresp_mimo_plan
    end interface
    interface bode_mimo
        module procedure bode_mimo_real, bode_mimo_complex, bode_mimo_plan
    end interface
    interface dbode_mimo
        module procedure dbode_mimo_real, dbode_mimo_complex, dbode_mimo_plan
    end interface

    interface freqresp
        module procedure freqresp_real_tf, freqresp_complex_tf
        module procedure freqresp_real_ss, freqresp_complex_ss
    end interface
    interface bode
        module procedure bode_real_tf, bode_complex_tf
        module procedure bode_real_ss, bode_complex_ss
    end interface
    interface dfreqresp
        module procedure dfreqresp_real_tf, dfreqresp_complex_tf
        module procedure dfreqresp_real_ss, dfreqresp_complex_ss
    end interface
    interface dbode
        module procedure dbode_real_tf, dbode_complex_tf
        module procedure dbode_real_ss, dbode_complex_ss
    end interface
    interface freqresp_zpk
        module procedure freqresp_zpk_rr, freqresp_zpk_rc
        module procedure freqresp_zpk_cr, freqresp_zpk_cc
    end interface
    interface dfreqresp_zpk
        module procedure dfreqresp_zpk_rr, dfreqresp_zpk_rc
        module procedure dfreqresp_zpk_cr, dfreqresp_zpk_cc
    end interface
    interface bode_zpk
        module procedure bode_zpk_rr, bode_zpk_rc
        module procedure bode_zpk_cr, bode_zpk_cc
    end interface
    interface dbode_zpk
        module procedure dbode_zpk_rr, dbode_zpk_rc
        module procedure dbode_zpk_cr, dbode_zpk_cc
    end interface

    public :: freqresp, bode, dfreqresp, dbode
    public :: freqresp_mimo, bode_mimo, dfreqresp_mimo, dbode_mimo
    public :: prepare_mimo_frequency
    public :: freqresp_zpk, bode_zpk, dfreqresp_zpk, dbode_zpk

contains

    pure logical function valid_dt(dt) result(ok)
        real(dp), intent(in), optional :: dt !! positive sampling interval (seconds)
        ok = .true.
        if (present(dt)) ok = ieee_is_finite(dt) .and. dt > 0.0_dp
    end function valid_dt

    ! Evaluate descending-power z polynomials via the existing ascending-power
    ! delay-coefficient evaluator.  The common degree padding is essential:
    ! H(z) = b(z)/a(z) = bpad(z**-1)/apad(z**-1).
    function discrete_tf(b, a, w, n, whole, dt) result(out)
        complex(dp), intent(in) :: b(:) !! numerator in descending z powers
        complex(dp), intent(in) :: a(:) !! denominator in descending z powers
        real(dp), intent(in), optional :: w(:) !! angular frequencies, radians/sample
        integer, intent(in), optional :: n !! automatic frequency-bin count (default 10000)
        logical, intent(in), optional :: whole !! entire unit circle, default false
        real(dp), intent(in), optional :: dt !! positive sample interval for validation
        type(frequency_response) :: out
        complex(dp), allocatable :: bp(:), ap(:)
        integer :: degree, count

        if (.not. valid_dt(dt)) return
        if (size(b) == 0 .or. size(a) == 0) return
        degree = max(size(b),size(a))
        allocate(bp(degree),ap(degree))
        bp = (0.0_dp,0.0_dp)
        ap = (0.0_dp,0.0_dp)
        bp(degree-size(b)+1:) = b
        ap(degree-size(a)+1:) = a
        count = 10000
        if (present(n)) count = n
        out = freqz(bp,ap,worn=count,whole=whole,w=w)
    end function discrete_tf

    function freqresp_real_tf(b, a, w, n) result(out)
        real(dp), intent(in) :: b(:) !! numerator coefficients, descending s powers
        real(dp), intent(in) :: a(:) !! denominator coefficients, descending s powers
        real(dp), intent(in), optional :: w(:) !! angular frequencies, rad/s
        integer, intent(in), optional :: n !! automatic log-frequency count, default 10000
        type(frequency_response) :: out
        integer :: count
        count = 10000
        if (present(n)) count = n
        out = freqs(b,a,worn=count,w=w)
    end function freqresp_real_tf

    function freqresp_complex_tf(b, a, w, n) result(out)
        complex(dp), intent(in) :: b(:) !! numerator coefficients, descending s powers
        complex(dp), intent(in) :: a(:) !! denominator coefficients, descending s powers
        real(dp), intent(in), optional :: w(:) !! angular frequencies, rad/s
        integer, intent(in), optional :: n !! automatic log-frequency count, default 10000
        type(frequency_response) :: out
        integer :: count
        count = 10000
        if (present(n)) count = n
        out = freqs(b,a,worn=count,w=w)
    end function freqresp_complex_tf

    function dfreqresp_real_tf(b, a, w, n, whole, dt) result(out)
        real(dp), intent(in) :: b(:) !! numerator coefficients, descending z powers
        real(dp), intent(in) :: a(:) !! denominator coefficients, descending z powers
        real(dp), intent(in), optional :: w(:) !! angular frequencies, radians/sample
        integer, intent(in), optional :: n !! automatic frequency count, default 10000
        logical, intent(in), optional :: whole !! full circle rather than [0,pi)
        real(dp), intent(in), optional :: dt !! positive sample interval
        type(frequency_response) :: out
        out = discrete_tf(cmplx(b,0.0_dp,dp),cmplx(a,0.0_dp,dp),w,n,whole,dt)
    end function dfreqresp_real_tf

    function dfreqresp_complex_tf(b, a, w, n, whole, dt) result(out)
        complex(dp), intent(in) :: b(:) !! numerator coefficients, descending z powers
        complex(dp), intent(in) :: a(:) !! denominator coefficients, descending z powers
        real(dp), intent(in), optional :: w(:) !! angular frequencies, radians/sample
        integer, intent(in), optional :: n !! automatic frequency count, default 10000
        logical, intent(in), optional :: whole !! full circle rather than [0,pi)
        real(dp), intent(in), optional :: dt !! positive sample interval
        type(frequency_response) :: out
        out = discrete_tf(b,a,w,n,whole,dt)
    end function dfreqresp_complex_tf

    ! State-space frequency responses are evaluated using shifted linear
    ! systems rather than expanding det(z*I-A) and C*adj(z*I-A)*B.
    ! Polynomial coefficients become extremely ill-conditioned at high order.
    pure logical function finite_complex_matrix(x) result(ok)
        complex(dp), intent(in) :: x(:,:) !! matrix to validate
        ok = all(ieee_is_finite(real(x,dp))) .and. all(ieee_is_finite(aimag(x)))
    end function finite_complex_matrix

    subroutine shifted_solve(a, b, z, x, ok)
        ! Solve (z*I-A)*x=b by Gaussian elimination with partial pivoting.
        ! No polynomial coefficients, determinants, or external LAPACK needed.
        complex(dp), intent(in) :: a(:,:) !! square state matrix
        complex(dp), intent(in) :: b(:) !! right hand side
        complex(dp), intent(in) :: z !! complex Laplace or z-domain argument
        complex(dp), intent(out) :: x(:) !! solution if ok
        logical, intent(out) :: ok !! true when a finite solution exists
        complex(dp), allocatable :: work(:,:), rhs(:)
        complex(dp) :: factor, tmp
        integer :: nstate, k, row, pivot

        ok = .false.
        nstate = size(b)
        allocate(work(nstate,nstate),rhs(nstate))
        work = -a
        rhs = b
        do k = 1, nstate
            work(k,k) = work(k,k)+z
        end do
        do k = 1, nstate
            pivot = k-1+maxloc(abs(work(k:nstate,k)),dim=1)
            if (abs(work(pivot,k)) <= tiny(1.0_dp)) return
            if (pivot /= k) then
                do row = k,nstate
                    tmp = work(k,row)
                    work(k,row) = work(pivot,row)
                    work(pivot,row) = tmp
                end do
                tmp = rhs(k)
                rhs(k) = rhs(pivot)
                rhs(pivot) = tmp
            end if
            do row = k+1,nstate
                factor = work(row,k)/work(k,k)
                work(row,k) = (0.0_dp,0.0_dp)
                work(row,k+1:nstate) = work(row,k+1:nstate)-factor*work(k,k+1:nstate)
                rhs(row) = rhs(row)-factor*rhs(k)
            end do
        end do
        do k = nstate,1,-1
            tmp = rhs(k)
            if (k < nstate) tmp = tmp-sum(work(k,k+1:nstate)*x(k+1:nstate))
            x(k) = tmp/work(k,k)
        end do
        ok = all(ieee_is_finite(real(x,dp))) .and. all(ieee_is_finite(aimag(x)))
    end subroutine shifted_solve


    subroutine shifted_solve_multiple(a, b, z, x, ok)
        ! Factor z*I-A once and solve simultaneously for every input column.
        ! Partial pivoting is shared across right-hand sides; no LAPACK needed.
        complex(dp), intent(in) :: a(:,:) !! state matrix (nstate,nstate)
        complex(dp), intent(in) :: b(:,:) !! input columns (nstate,ninput)
        complex(dp), intent(in) :: z !! s or z evaluation argument
        complex(dp), intent(out) :: x(:,:) !! solutions (nstate,ninput)
        logical, intent(out) :: ok !! true if solve is finite and nonsingular
        complex(dp), allocatable :: work(:,:), rhs(:,:), swap_row(:), swap_rhs(:)
        complex(dp) :: factor
        integer :: ns, ni, k, row, pivot, j

        ok = .false.
        ns = size(a,1)
        ni = size(b,2)
        allocate(work(ns,ns), rhs(ns,ni), swap_row(ns), swap_rhs(ni))
        work = -a
        rhs = b
        do k = 1,ns
            work(k,k) = work(k,k)+z
        end do
        do k = 1,ns
            pivot = k-1+maxloc(abs(work(k:ns,k)),dim=1)
            if (abs(work(pivot,k)) <= tiny(1.0_dp)) return
            if (pivot /= k) then
                swap_row = work(k,:)
                work(k,:) = work(pivot,:)
                work(pivot,:) = swap_row
                swap_rhs = rhs(k,:)
                rhs(k,:) = rhs(pivot,:)
                rhs(pivot,:) = swap_rhs
            end if
            do row = k+1,ns
                factor = work(row,k)/work(k,k)
                work(row,k) = (0.0_dp,0.0_dp)
                work(row,k+1:ns) = work(row,k+1:ns)-factor*work(k,k+1:ns)
                rhs(row,:) = rhs(row,:)-factor*rhs(k,:)
            end do
        end do
        x = (0.0_dp,0.0_dp)
        do k = ns,1,-1
            x(k,:) = rhs(k,:)
            do j = k+1,ns
                x(k,:) = x(k,:)-work(k,j)*x(j,:)
            end do
            x(k,:) = x(k,:)/work(k,k)
        end do
        ok = finite_complex_matrix(x)
    end subroutine shifted_solve_multiple

    subroutine reduce_to_hessenberg(a, b, c, h, br, cr, ok)
        ! Unitary Householder similarity: h=Q**H*A*Q, br=Q**H*B, cr=C*Q.
        ! A is reduced once for a complete frequency grid.  In contrast to
        ! an eigenvector basis, unitary Q does not amplify conditioning.
        complex(dp), intent(in) :: a(:,:) !! original state matrix
        complex(dp), intent(in) :: b(:,:) !! original input matrix
        complex(dp), intent(in) :: c(:,:) !! original output matrix
        complex(dp), allocatable, intent(out) :: h(:,:) !! upper Hessenberg matrix
        complex(dp), allocatable, intent(out) :: br(:,:) !! transformed inputs
        complex(dp), allocatable, intent(out) :: cr(:,:) !! transformed outputs
        logical, intent(out) :: ok !! true if all reduced matrices are finite
        complex(dp), allocatable :: v(:), column(:)
        complex(dp) :: phase, projection
        real(dp) :: scale, norm_v, tau
        integer :: ns, ni, no, k, i, j

        h = a
        br = b
        cr = c
        ns = size(a,1)
        ni = size(b,2)
        no = size(c,1)
        do k = 1, ns-2
            scale = maxval(abs(h(k+1:ns,k)))
            if (scale == 0.0_dp) cycle
            column = h(k+1:ns,k)/scale
            norm_v = sqrt(sum(abs(column)**2))
            phase = (1.0_dp,0.0_dp)
            if (abs(column(1)) > 0.0_dp) phase = column(1)/abs(column(1))
            v = column
            v(1) = v(1)+phase*norm_v
            tau = 2.0_dp/real(dot_product(v,v),dp)
            ! Left application to A and B; right application to A and C.
            do j = k,ns
                projection = sum(conjg(v)*h(k+1:ns,j))
                h(k+1:ns,j) = h(k+1:ns,j)-tau*v*projection
            end do
            do i = 1,ns
                projection = sum(h(i,k+1:ns)*v)
                h(i,k+1:ns) = h(i,k+1:ns)-tau*projection*conjg(v)
            end do
            do j = 1,ni
                projection = sum(conjg(v)*br(k+1:ns,j))
                br(k+1:ns,j) = br(k+1:ns,j)-tau*v*projection
            end do
            do i = 1,no
                projection = sum(cr(i,k+1:ns)*v)
                cr(i,k+1:ns) = cr(i,k+1:ns)-tau*projection*conjg(v)
            end do
            h(k+2:ns,k) = (0.0_dp,0.0_dp)
        end do
        ok = finite_complex_matrix(h) .and. finite_complex_matrix(br) .and. &
            finite_complex_matrix(cr)
    end subroutine reduce_to_hessenberg

    subroutine hessenberg_solve_multiple(h, b, z, x, ok)
        ! Pivoted solve of (z*I-H)*X=B for upper Hessenberg H.
        ! Only adjacent rows require pivoting; O(n**2) work per frequency
        ! plus O(n**2*nin) for multiple right-hand sides.
        complex(dp), intent(in) :: h(:,:) !! Hessenberg state matrix
        complex(dp), intent(in) :: b(:,:) !! input columns
        complex(dp), intent(in) :: z !! complex frequency argument
        complex(dp), intent(out) :: x(:,:) !! state responses
        logical, intent(out) :: ok !! solution available and finite
        complex(dp), allocatable :: work(:,:), rhs(:,:), saved_row(:), saved_rhs(:)
        complex(dp) :: factor
        integer :: ns, ni, k, j

        ok = .false.
        ns = size(h,1)
        ni = size(b,2)
        work = -h
        rhs = b
        allocate(saved_row(ns),saved_rhs(ni))
        do k = 1,ns
            work(k,k) = work(k,k)+z
        end do
        do k = 1,ns-1
            if (abs(work(k+1,k)) > abs(work(k,k))) then
                saved_row(k:ns) = work(k,k:ns)
                work(k,k:ns) = work(k+1,k:ns)
                work(k+1,k:ns) = saved_row(k:ns)
                saved_rhs = rhs(k,:)
                rhs(k,:) = rhs(k+1,:)
                rhs(k+1,:) = saved_rhs
            end if
            if (abs(work(k,k)) <= tiny(1.0_dp)) return
            factor = work(k+1,k)/work(k,k)
            work(k+1,k) = (0.0_dp,0.0_dp)
            work(k+1,k+1:ns) = work(k+1,k+1:ns)-factor*work(k,k+1:ns)
            rhs(k+1,:) = rhs(k+1,:)-factor*rhs(k,:)
        end do
        if (ns > 0) then
            if (abs(work(ns,ns)) <= tiny(1.0_dp)) return
        end if
        x = (0.0_dp,0.0_dp)
        do k = ns,1,-1
            x(k,:) = rhs(k,:)
            do j = k+1,ns
                x(k,:) = x(k,:)-work(k,j)*x(j,:)
            end do
            x(k,:) = x(k,:)/work(k,k)
        end do
        ok = finite_complex_matrix(x)
    end subroutine hessenberg_solve_multiple

    subroutine state_space_analog_grid(a, count, grid, log_range)
        ! A pole-scale estimate avoids characteristic-polynomial formation.
        ! The grid is logarithmic, covering nonzero row and diagonal scales.
        complex(dp), intent(in) :: a(:,:) !! state matrix
        integer, intent(in) :: count !! requested frequency count, at least 1
        real(dp), allocatable, intent(out) :: grid(:) !! frequencies in rad/s
        real(dp), intent(out), optional :: log_range(2) !! grid decade bounds
        real(dp) :: low, high, scale, magnitude
        integer :: i, j, nstate

        nstate = size(a,1)
        low = huge(1.0_dp)
        high = 0.0_dp
        do i = 1,nstate
            scale = 0.0_dp
            do j = 1,nstate
                magnitude = abs(a(i,j))
                scale = max(scale,magnitude)
                if (magnitude > 0.0_dp) low = min(low,magnitude)
            end do
            high = max(high,scale)
        end do
        if (high == 0.0_dp) then
            low = 1.0_dp
            high = 1.0_dp
        end if
        ! Clamp decades to finite powers of ten (including extreme scales).
        low = max(-300.0_dp, min(300.0_dp, floor(log10(low))-1.0_dp))
        high = max(low, max(-300.0_dp, min(300.0_dp, ceiling(log10(high))+1.0_dp)))
        if (present(log_range)) log_range = [low,high]
        call fill_mimo_log_grid([low,high],count,grid)
    end subroutine state_space_analog_grid

    subroutine fill_mimo_log_grid(log_range,count,grid)
        real(dp), intent(in) :: log_range(2) !! start and end decades
        integer, intent(in) :: count !! at least one frequency
        real(dp), allocatable, intent(out) :: grid(:) !! generated frequencies
        real(dp) :: width
        integer :: i
        width = max(0.0_dp,log_range(2)-log_range(1))
        allocate(grid(count))
        do i = 1,count
            if (count == 1) then
                grid(i) = 10.0_dp**log_range(1)
            else
                grid(i) = 10.0_dp**(log_range(1)+width*real(i-1,dp)/real(count-1,dp))
            end if
        end do
    end subroutine fill_mimo_log_grid

    function state_space_response(a,b,c,d,w,n,discrete,whole,dt) result(out)
        complex(dp), intent(in) :: a(:,:) !! state matrix
        complex(dp), intent(in) :: b(:,:) !! input matrix, one column
        complex(dp), intent(in) :: c(:,:) !! output matrix, one row
        complex(dp), intent(in) :: d(:,:) !! direct feedthrough, 1 by 1
        real(dp), intent(in), optional :: w(:) !! explicit angular frequencies
        integer, intent(in), optional :: n !! automatic frequency count
        logical, intent(in) :: discrete !! false: analog, true: digital
        logical, intent(in), optional :: whole !! full unit circle for digital grids
        real(dp), intent(in), optional :: dt !! sampling interval for validation
        type(frequency_response) :: out
        complex(dp), allocatable :: x(:)
        complex(dp) :: z
        real(dp) :: nan
        integer :: i, count, ns
        logical :: ok, full_circle

        ns = size(a,1)
        if (size(a,2) /= ns .or. size(b,1) /= ns .or. size(b,2) /= 1) return
        if (size(c,1) /= 1 .or. size(c,2) /= ns) return
        if (size(d,1) /= 1 .or. size(d,2) /= 1) return
        if (.not. finite_complex_matrix(a) .or. .not. finite_complex_matrix(b)) return
        if (.not. finite_complex_matrix(c) .or. .not. finite_complex_matrix(d)) return
        if (discrete .and. .not. valid_dt(dt)) return
        if (present(w)) then
            if (.not. all(ieee_is_finite(w))) return
            out%frequencies = w
        else
            count = 10000
            if (present(n)) count = n
            if (count < 1) return
            if (discrete) then
                full_circle = .false.
                if (present(whole)) full_circle = whole
                allocate(out%frequencies(count))
                do i = 1,count
                    if (full_circle) then
                        out%frequencies(i) = 2.0_dp*pi*real(i-1,dp)/real(count,dp)
                    else
                        out%frequencies(i) = pi*real(i-1,dp)/real(count,dp)
                    end if
                end do
            else
                call state_space_analog_grid(a,count,out%frequencies)
            end if
        end if
        allocate(out%response(size(out%frequencies)),x(ns))
        out%status = design_status_success
        nan = ieee_value(0.0_dp,ieee_quiet_nan)
        do i = 1,size(out%frequencies)
            z = cmplx(0.0_dp,out%frequencies(i),dp)
            if (discrete) z = exp(z)
            if (ns == 0) then
                out%response(i) = d(1,1)
            else
                call shifted_solve(a,b(:,1),z,x,ok)
                if (ok) then
                    out%response(i) = d(1,1)+sum(c(1,:)*x)
                else
                    out%response(i) = cmplx(nan,nan,dp)
                    out%status = design_status_singular
                end if
            end if
            if (.not. ieee_is_finite(real(out%response(i),dp)) .or. &
                .not. ieee_is_finite(aimag(out%response(i)))) then
                out%status = design_status_singular
                out%response(i) = cmplx(nan,nan,dp)
            end if
        end do
    end function state_space_response

    function freqresp_real_ss(a, b, c, d, w, n) result(out)
        real(dp), intent(in) :: a(:,:) !! state matrix
        real(dp), intent(in) :: b(:,:) !! input matrix
        real(dp), intent(in) :: c(:,:) !! output matrix
        real(dp), intent(in) :: d(:,:) !! direct feedthrough
        real(dp), intent(in), optional :: w(:) !! angular frequencies, rad/s
        integer, intent(in), optional :: n !! automatic frequency count
        type(frequency_response) :: out
        out = state_space_response(cmplx(a,0.0_dp,dp),cmplx(b,0.0_dp,dp), &
            cmplx(c,0.0_dp,dp),cmplx(d,0.0_dp,dp),w,n,.false.)
    end function freqresp_real_ss

    function freqresp_complex_ss(a, b, c, d, w, n) result(out)
        complex(dp), intent(in) :: a(:,:) !! state matrix
        complex(dp), intent(in) :: b(:,:) !! input matrix
        complex(dp), intent(in) :: c(:,:) !! output matrix
        complex(dp), intent(in) :: d(:,:) !! direct feedthrough
        real(dp), intent(in), optional :: w(:) !! angular frequencies, rad/s
        integer, intent(in), optional :: n !! automatic frequency count
        type(frequency_response) :: out
        out = state_space_response(a,b,c,d,w,n,.false.)
    end function freqresp_complex_ss

    function dfreqresp_real_ss(a, b, c, d, w, n, whole, dt) result(out)
        real(dp), intent(in) :: a(:,:) !! discrete state matrix
        real(dp), intent(in) :: b(:,:) !! input matrix
        real(dp), intent(in) :: c(:,:) !! output matrix
        real(dp), intent(in) :: d(:,:) !! direct feedthrough
        real(dp), intent(in), optional :: w(:) !! angular frequencies, radians/sample
        integer, intent(in), optional :: n !! automatic frequency count
        logical, intent(in), optional :: whole !! full unit-circle response
        real(dp), intent(in), optional :: dt !! positive sampling interval
        type(frequency_response) :: out
        out = state_space_response(cmplx(a,0.0_dp,dp),cmplx(b,0.0_dp,dp), &
            cmplx(c,0.0_dp,dp),cmplx(d,0.0_dp,dp),w,n,.true.,whole,dt)
    end function dfreqresp_real_ss

    function dfreqresp_complex_ss(a, b, c, d, w, n, whole, dt) result(out)
        complex(dp), intent(in) :: a(:,:) !! discrete state matrix
        complex(dp), intent(in) :: b(:,:) !! input matrix
        complex(dp), intent(in) :: c(:,:) !! output matrix
        complex(dp), intent(in) :: d(:,:) !! direct feedthrough
        real(dp), intent(in), optional :: w(:) !! angular frequencies, radians/sample
        integer, intent(in), optional :: n !! automatic frequency count
        logical, intent(in), optional :: whole !! full unit-circle response
        real(dp), intent(in), optional :: dt !! positive sampling interval
        type(frequency_response) :: out
        out = state_space_response(a,b,c,d,w,n,.true.,whole,dt)
    end function dfreqresp_complex_ss

    function bode_from_response(response, dt) result(out)
        type(frequency_response), intent(in) :: response !! complex frequency response
        real(dp), intent(in), optional :: dt !! scale output axis to rad/s
        type(bode_result) :: out
        real(dp) :: theta, previous, delta, corrected, current, scale
        integer :: i

        out%status = response%status
        if (response%status /= design_status_success) return
        if (.not. allocated(response%frequencies)) return
        if (.not. allocated(response%response)) return
        if (.not. valid_dt(dt)) then
            out%status = design_status_invalid_input
            return
        end if
        scale = 1.0_dp
        if (present(dt)) scale = dt
        out%frequencies = response%frequencies / scale
        allocate(out%magnitude(size(response%response)),out%phase(size(response%response)))
        previous = 0.0_dp
        theta = 0.0_dp
        do i = 1,size(response%response)
            current = atan2(aimag(response%response(i)),real(response%response(i),dp))
            if (i == 1) then
                theta = current
            else
                delta = current-previous
                corrected = modulo(delta+pi,2.0_dp*pi)-pi
                if (corrected == -pi .and. delta > 0.0_dp) corrected = pi
                theta = theta+corrected
            end if
            out%phase(i) = theta*180.0_dp/pi
            if (abs(response%response(i)) == 0.0_dp) then
                out%magnitude(i) = ieee_value(0.0_dp,ieee_negative_inf)
            else
                out%magnitude(i) = 20.0_dp*log10(abs(response%response(i)))
            end if
            previous = current
        end do
    end function bode_from_response
    function bode_real_tf(b, a, w, n) result(out)
        real(dp), intent(in) :: b(:) !! transfer coefficients
        real(dp), intent(in) :: a(:) !! transfer coefficients
        real(dp), intent(in), optional :: w(:) !! evaluation frequency grid
        integer, intent(in), optional :: n !! automatic grid size, default 100
        type(bode_result) :: out
        type(frequency_response) :: response
        integer :: count
        count = 100
        if (present(n)) count = n
        response = freqresp(b,a,w=w,n=count)
        out = bode_from_response(response)
    end function bode_real_tf

    function dbode_real_tf(b, a, w, n, dt) result(out)
        real(dp), intent(in) :: b(:) !! transfer coefficients
        real(dp), intent(in) :: a(:) !! transfer coefficients
        real(dp), intent(in), optional :: w(:) !! evaluation frequency grid
        integer, intent(in), optional :: n !! automatic grid size, default 100
        real(dp), intent(in), optional :: dt !! positive sample interval, default 1
        type(bode_result) :: out
        type(frequency_response) :: response
        integer :: count
        count = 100
        if (present(n)) count = n
        response = dfreqresp(b,a,w=w,n=count,dt=dt)
        out = bode_from_response(response,dt)
    end function dbode_real_tf

    function bode_real_ss(a, b, c, d, w, n) result(out)
        real(dp), intent(in) :: a(:,:) !! state-space matrix
        real(dp), intent(in) :: b(:,:) !! state-space matrix
        real(dp), intent(in) :: c(:,:) !! state-space matrix
        real(dp), intent(in) :: d(:,:) !! state-space matrix
        real(dp), intent(in), optional :: w(:) !! evaluation frequency grid
        integer, intent(in), optional :: n !! automatic grid size, default 100
        type(bode_result) :: out
        type(frequency_response) :: response
        integer :: count
        count = 100
        if (present(n)) count = n
        response = freqresp(a,b,c,d,w=w,n=count)
        out = bode_from_response(response)
    end function bode_real_ss

    function dbode_real_ss(a, b, c, d, w, n, dt) result(out)
        real(dp), intent(in) :: a(:,:) !! state-space matrix
        real(dp), intent(in) :: b(:,:) !! state-space matrix
        real(dp), intent(in) :: c(:,:) !! state-space matrix
        real(dp), intent(in) :: d(:,:) !! state-space matrix
        real(dp), intent(in), optional :: w(:) !! evaluation frequency grid
        integer, intent(in), optional :: n !! automatic grid size, default 100
        real(dp), intent(in), optional :: dt !! positive sample interval, default 1
        type(bode_result) :: out
        type(frequency_response) :: response
        integer :: count
        count = 100
        if (present(n)) count = n
        response = dfreqresp(a,b,c,d,w=w,n=count,dt=dt)
        out = bode_from_response(response,dt)
    end function dbode_real_ss

    function bode_complex_tf(b, a, w, n) result(out)
        complex(dp), intent(in) :: b(:) !! transfer coefficients
        complex(dp), intent(in) :: a(:) !! transfer coefficients
        real(dp), intent(in), optional :: w(:) !! evaluation frequency grid
        integer, intent(in), optional :: n !! automatic grid size, default 100
        type(bode_result) :: out
        type(frequency_response) :: response
        integer :: count
        count = 100
        if (present(n)) count = n
        response = freqresp(b,a,w=w,n=count)
        out = bode_from_response(response)
    end function bode_complex_tf

    function dbode_complex_tf(b, a, w, n, dt) result(out)
        complex(dp), intent(in) :: b(:) !! transfer coefficients
        complex(dp), intent(in) :: a(:) !! transfer coefficients
        real(dp), intent(in), optional :: w(:) !! evaluation frequency grid
        integer, intent(in), optional :: n !! automatic grid size, default 100
        real(dp), intent(in), optional :: dt !! positive sample interval, default 1
        type(bode_result) :: out
        type(frequency_response) :: response
        integer :: count
        count = 100
        if (present(n)) count = n
        response = dfreqresp(b,a,w=w,n=count,dt=dt)
        out = bode_from_response(response,dt)
    end function dbode_complex_tf

    function bode_complex_ss(a, b, c, d, w, n) result(out)
        complex(dp), intent(in) :: a(:,:) !! state-space matrix
        complex(dp), intent(in) :: b(:,:) !! state-space matrix
        complex(dp), intent(in) :: c(:,:) !! state-space matrix
        complex(dp), intent(in) :: d(:,:) !! state-space matrix
        real(dp), intent(in), optional :: w(:) !! evaluation frequency grid
        integer, intent(in), optional :: n !! automatic grid size, default 100
        type(bode_result) :: out
        type(frequency_response) :: response
        integer :: count
        count = 100
        if (present(n)) count = n
        response = freqresp(a,b,c,d,w=w,n=count)
        out = bode_from_response(response)
    end function bode_complex_ss

    function dbode_complex_ss(a, b, c, d, w, n, dt) result(out)
        complex(dp), intent(in) :: a(:,:) !! state-space matrix
        complex(dp), intent(in) :: b(:,:) !! state-space matrix
        complex(dp), intent(in) :: c(:,:) !! state-space matrix
        complex(dp), intent(in) :: d(:,:) !! state-space matrix
        real(dp), intent(in), optional :: w(:) !! evaluation frequency grid
        integer, intent(in), optional :: n !! automatic grid size, default 100
        real(dp), intent(in), optional :: dt !! positive sample interval, default 1
        type(bode_result) :: out
        type(frequency_response) :: response
        integer :: count
        count = 100
        if (present(n)) count = n
        response = dfreqresp(a,b,c,d,w=w,n=count,dt=dt)
        out = bode_from_response(response,dt)
    end function dbode_complex_ss

    function freqresp_zpk_rr(z,p,k,w,n) result(out)
        real(dp), intent(in) :: z(:) !! system zeros
        real(dp), intent(in) :: p(:) !! system poles
        real(dp), intent(in) :: k !! scalar system gain
        real(dp), intent(in), optional :: w(:) !! angular frequency grid
        integer, intent(in), optional :: n !! frequency point count, default 10000
        type(frequency_response) :: out
        integer :: count
        count = 10000
        if (present(n)) count = n
        out = freqs_zpk(z,p,k,worn=count,w=w)
    end function freqresp_zpk_rr

    function bode_zpk_rr(z,p,k,w,n) result(out)
        real(dp), intent(in) :: z(:) !! system zeros
        real(dp), intent(in) :: p(:) !! system poles
        real(dp), intent(in) :: k !! scalar system gain
        real(dp), intent(in), optional :: w(:) !! angular frequency grid
        integer, intent(in), optional :: n !! frequency point count, default 100
        type(bode_result) :: out
        type(frequency_response) :: response
        integer :: count
        count = 100
        if (present(n)) count = n
        response = freqresp_zpk(z,p,k,w=w,n=count)
        out = bode_from_response(response)
    end function bode_zpk_rr

    function dfreqresp_zpk_rr(z,p,k,w,n,whole,dt) result(out)
        real(dp), intent(in) :: z(:) !! system zeros
        real(dp), intent(in) :: p(:) !! system poles
        real(dp), intent(in) :: k !! scalar system gain
        real(dp), intent(in), optional :: w(:) !! angular frequency grid
        integer, intent(in), optional :: n !! frequency point count, default 10000
        logical, intent(in), optional :: whole !! evaluate complete unit circle
        real(dp), intent(in), optional :: dt !! positive sampling interval
        type(frequency_response) :: out
        integer :: count
        count = 10000
        if (present(n)) count = n
        if (.not. valid_dt(dt)) return
        out = freqz_zpk(z,p,k,worn=count,whole=whole,fs=2.0_dp*pi,w=w)
    end function dfreqresp_zpk_rr

    function dbode_zpk_rr(z,p,k,w,n,dt) result(out)
        real(dp), intent(in) :: z(:) !! system zeros
        real(dp), intent(in) :: p(:) !! system poles
        real(dp), intent(in) :: k !! scalar system gain
        real(dp), intent(in), optional :: w(:) !! angular frequency grid
        integer, intent(in), optional :: n !! frequency point count, default 100
        real(dp), intent(in), optional :: dt !! sampling interval, default 1
        type(bode_result) :: out
        type(frequency_response) :: response
        integer :: count
        count = 100
        if (present(n)) count = n
        response = dfreqresp_zpk(z,p,k,w=w,n=count,dt=dt)
        out = bode_from_response(response,dt)
    end function dbode_zpk_rr

    function freqresp_zpk_rc(z,p,k,w,n) result(out)
        real(dp), intent(in) :: z(:) !! system zeros
        real(dp), intent(in) :: p(:) !! system poles
        complex(dp), intent(in) :: k !! scalar system gain
        real(dp), intent(in), optional :: w(:) !! angular frequency grid
        integer, intent(in), optional :: n !! frequency point count, default 10000
        type(frequency_response) :: out
        integer :: count
        count = 10000
        if (present(n)) count = n
        out = freqs_zpk(z,p,k,worn=count,w=w)
    end function freqresp_zpk_rc

    function bode_zpk_rc(z,p,k,w,n) result(out)
        real(dp), intent(in) :: z(:) !! system zeros
        real(dp), intent(in) :: p(:) !! system poles
        complex(dp), intent(in) :: k !! scalar system gain
        real(dp), intent(in), optional :: w(:) !! angular frequency grid
        integer, intent(in), optional :: n !! frequency point count, default 100
        type(bode_result) :: out
        type(frequency_response) :: response
        integer :: count
        count = 100
        if (present(n)) count = n
        response = freqresp_zpk(z,p,k,w=w,n=count)
        out = bode_from_response(response)
    end function bode_zpk_rc

    function dfreqresp_zpk_rc(z,p,k,w,n,whole,dt) result(out)
        real(dp), intent(in) :: z(:) !! system zeros
        real(dp), intent(in) :: p(:) !! system poles
        complex(dp), intent(in) :: k !! scalar system gain
        real(dp), intent(in), optional :: w(:) !! angular frequency grid
        integer, intent(in), optional :: n !! frequency point count, default 10000
        logical, intent(in), optional :: whole !! evaluate complete unit circle
        real(dp), intent(in), optional :: dt !! positive sampling interval
        type(frequency_response) :: out
        integer :: count
        count = 10000
        if (present(n)) count = n
        if (.not. valid_dt(dt)) return
        out = freqz_zpk(z,p,k,worn=count,whole=whole,fs=2.0_dp*pi,w=w)
    end function dfreqresp_zpk_rc

    function dbode_zpk_rc(z,p,k,w,n,dt) result(out)
        real(dp), intent(in) :: z(:) !! system zeros
        real(dp), intent(in) :: p(:) !! system poles
        complex(dp), intent(in) :: k !! scalar system gain
        real(dp), intent(in), optional :: w(:) !! angular frequency grid
        integer, intent(in), optional :: n !! frequency point count, default 100
        real(dp), intent(in), optional :: dt !! sampling interval, default 1
        type(bode_result) :: out
        type(frequency_response) :: response
        integer :: count
        count = 100
        if (present(n)) count = n
        response = dfreqresp_zpk(z,p,k,w=w,n=count,dt=dt)
        out = bode_from_response(response,dt)
    end function dbode_zpk_rc

    function freqresp_zpk_cr(z,p,k,w,n) result(out)
        complex(dp), intent(in) :: z(:) !! system zeros
        complex(dp), intent(in) :: p(:) !! system poles
        real(dp), intent(in) :: k !! scalar system gain
        real(dp), intent(in), optional :: w(:) !! angular frequency grid
        integer, intent(in), optional :: n !! frequency point count, default 10000
        type(frequency_response) :: out
        integer :: count
        count = 10000
        if (present(n)) count = n
        out = freqs_zpk(z,p,k,worn=count,w=w)
    end function freqresp_zpk_cr

    function bode_zpk_cr(z,p,k,w,n) result(out)
        complex(dp), intent(in) :: z(:) !! system zeros
        complex(dp), intent(in) :: p(:) !! system poles
        real(dp), intent(in) :: k !! scalar system gain
        real(dp), intent(in), optional :: w(:) !! angular frequency grid
        integer, intent(in), optional :: n !! frequency point count, default 100
        type(bode_result) :: out
        type(frequency_response) :: response
        integer :: count
        count = 100
        if (present(n)) count = n
        response = freqresp_zpk(z,p,k,w=w,n=count)
        out = bode_from_response(response)
    end function bode_zpk_cr

    function dfreqresp_zpk_cr(z,p,k,w,n,whole,dt) result(out)
        complex(dp), intent(in) :: z(:) !! system zeros
        complex(dp), intent(in) :: p(:) !! system poles
        real(dp), intent(in) :: k !! scalar system gain
        real(dp), intent(in), optional :: w(:) !! angular frequency grid
        integer, intent(in), optional :: n !! frequency point count, default 10000
        logical, intent(in), optional :: whole !! evaluate complete unit circle
        real(dp), intent(in), optional :: dt !! positive sampling interval
        type(frequency_response) :: out
        integer :: count
        count = 10000
        if (present(n)) count = n
        if (.not. valid_dt(dt)) return
        out = freqz_zpk(z,p,k,worn=count,whole=whole,fs=2.0_dp*pi,w=w)
    end function dfreqresp_zpk_cr

    function dbode_zpk_cr(z,p,k,w,n,dt) result(out)
        complex(dp), intent(in) :: z(:) !! system zeros
        complex(dp), intent(in) :: p(:) !! system poles
        real(dp), intent(in) :: k !! scalar system gain
        real(dp), intent(in), optional :: w(:) !! angular frequency grid
        integer, intent(in), optional :: n !! frequency point count, default 100
        real(dp), intent(in), optional :: dt !! sampling interval, default 1
        type(bode_result) :: out
        type(frequency_response) :: response
        integer :: count
        count = 100
        if (present(n)) count = n
        response = dfreqresp_zpk(z,p,k,w=w,n=count,dt=dt)
        out = bode_from_response(response,dt)
    end function dbode_zpk_cr

    function freqresp_zpk_cc(z,p,k,w,n) result(out)
        complex(dp), intent(in) :: z(:) !! system zeros
        complex(dp), intent(in) :: p(:) !! system poles
        complex(dp), intent(in) :: k !! scalar system gain
        real(dp), intent(in), optional :: w(:) !! angular frequency grid
        integer, intent(in), optional :: n !! frequency point count, default 10000
        type(frequency_response) :: out
        integer :: count
        count = 10000
        if (present(n)) count = n
        out = freqs_zpk(z,p,k,worn=count,w=w)
    end function freqresp_zpk_cc

    function bode_zpk_cc(z,p,k,w,n) result(out)
        complex(dp), intent(in) :: z(:) !! system zeros
        complex(dp), intent(in) :: p(:) !! system poles
        complex(dp), intent(in) :: k !! scalar system gain
        real(dp), intent(in), optional :: w(:) !! angular frequency grid
        integer, intent(in), optional :: n !! frequency point count, default 100
        type(bode_result) :: out
        type(frequency_response) :: response
        integer :: count
        count = 100
        if (present(n)) count = n
        response = freqresp_zpk(z,p,k,w=w,n=count)
        out = bode_from_response(response)
    end function bode_zpk_cc

    function dfreqresp_zpk_cc(z,p,k,w,n,whole,dt) result(out)
        complex(dp), intent(in) :: z(:) !! system zeros
        complex(dp), intent(in) :: p(:) !! system poles
        complex(dp), intent(in) :: k !! scalar system gain
        real(dp), intent(in), optional :: w(:) !! angular frequency grid
        integer, intent(in), optional :: n !! frequency point count, default 10000
        logical, intent(in), optional :: whole !! evaluate complete unit circle
        real(dp), intent(in), optional :: dt !! positive sampling interval
        type(frequency_response) :: out
        integer :: count
        count = 10000
        if (present(n)) count = n
        if (.not. valid_dt(dt)) return
        out = freqz_zpk(z,p,k,worn=count,whole=whole,fs=2.0_dp*pi,w=w)
    end function dfreqresp_zpk_cc

    function dbode_zpk_cc(z,p,k,w,n,dt) result(out)
        complex(dp), intent(in) :: z(:) !! system zeros
        complex(dp), intent(in) :: p(:) !! system poles
        complex(dp), intent(in) :: k !! scalar system gain
        real(dp), intent(in), optional :: w(:) !! angular frequency grid
        integer, intent(in), optional :: n !! frequency point count, default 100
        real(dp), intent(in), optional :: dt !! sampling interval, default 1
        type(bode_result) :: out
        type(frequency_response) :: response
        integer :: count
        count = 100
        if (present(n)) count = n
        response = dfreqresp_zpk(z,p,k,w=w,n=count,dt=dt)
        out = bode_from_response(response,dt)
    end function dbode_zpk_cc


    subroutine complex_givens(f,g,cs,sn)
        ! Unitary rotation [cs,sn;-conjg(sn),cs] annihilates g.
        complex(dp), intent(in) :: f !! leading value
        complex(dp), intent(in) :: g !! value to eliminate
        real(dp), intent(out) :: cs !! real cosine
        complex(dp), intent(out) :: sn !! complex sine
        real(dp) :: scale, norm
        complex(dp) :: phase

        if (abs(g) == 0.0_dp) then
            cs = 1.0_dp
            sn = (0.0_dp,0.0_dp)
        else if (abs(f) == 0.0_dp) then
            cs = 0.0_dp
            sn = conjg(g)/abs(g)
        else
            scale = max(abs(f),abs(g))
            norm = scale*sqrt((abs(f)/scale)**2+(abs(g)/scale)**2)
            phase = f/abs(f)
            cs = abs(f)/norm
            sn = phase*conjg(g)/norm
        end if
    end subroutine complex_givens

    subroutine schur_from_hessenberg(h,b,c,t,bt,ct,ok)
        ! Shifted implicit single-shift QR sweeps on upper Hessenberg H.
        ! Each rotation is applied to T, Q**H*B and C*Q simultaneously;
        ! accumulated Q is never stored. Failure leaves the caller's H intact.
        complex(dp), intent(in) :: h(:,:) !! upper Hessenberg state matrix
        complex(dp), intent(in) :: b(:,:) !! transformed inputs
        complex(dp), intent(in) :: c(:,:) !! transformed outputs
        complex(dp), allocatable, intent(out) :: t(:,:) !! upper triangular Schur form
        complex(dp), allocatable, intent(out) :: bt(:,:) !! Schur-basis input matrix
        complex(dp), allocatable, intent(out) :: ct(:,:) !! Schur-basis output matrix
        logical, intent(out) :: ok !! true if QR iterations converged
        complex(dp) :: f, g, sn, mu, ev1, ev2, disc
        complex(dp), allocatable :: temp(:)
        real(dp) :: cs, threshold
        integer :: ns, lo, hi, j, k, steps, max_steps

        ok = .false.
        t = h
        bt = b
        ct = c
        ns = size(t,1)
        if (ns <= 1) then
            ok = .true.
            return
        end if
        max_steps = max(200,120*ns)
        hi = ns
        steps = 0
        do while (hi > 1)
            do j = 2,hi
                threshold = 32.0_dp*epsilon(1.0_dp)* &
                    (abs(t(j,j))+abs(t(j-1,j-1)))
                if (abs(t(j,j-1)) <= threshold) t(j,j-1) = (0.0_dp,0.0_dp)
            end do
            if (t(hi,hi-1) == (0.0_dp,0.0_dp)) then
                hi = hi-1
                cycle
            end if
            if (steps >= max_steps) return
            lo = hi-1
            do while (lo > 1)
                if (t(lo,lo-1) == (0.0_dp,0.0_dp)) exit
                lo = lo-1
            end do
            disc = sqrt((t(hi-1,hi-1)-t(hi,hi))**2+ &
                4.0_dp*t(hi-1,hi)*t(hi,hi-1))
            ev1 = 0.5_dp*(t(hi-1,hi-1)+t(hi,hi)+disc)
            ev2 = 0.5_dp*(t(hi-1,hi-1)+t(hi,hi)-disc)
            mu = ev1
            if (abs(ev2-t(hi,hi)) < abs(ev1-t(hi,hi))) mu = ev2
            if (mod(steps+1,40) == 0) mu = t(hi,hi)+0.7_dp*t(hi,hi-1)
            f = t(lo,lo)-mu
            g = t(lo+1,lo)
            do k = lo,hi-1
                if (k > lo) then
                    f = t(k,k-1)
                    g = t(k+1,k-1)
                end if
                call complex_givens(f,g,cs,sn)
                ! Similarity T <- G*T*G**H, B <- G*B, C <- C*G**H.
                temp = t(k,:)
                t(k,:) = cs*temp+sn*t(k+1,:)
                t(k+1,:) = -conjg(sn)*temp+cs*t(k+1,:)
                temp = t(:,k)
                t(:,k) = cs*temp+conjg(sn)*t(:,k+1)
                t(:,k+1) = -sn*temp+cs*t(:,k+1)
                temp = bt(k,:)
                bt(k,:) = cs*temp+sn*bt(k+1,:)
                bt(k+1,:) = -conjg(sn)*temp+cs*bt(k+1,:)
                temp = ct(:,k)
                ct(:,k) = cs*temp+conjg(sn)*ct(:,k+1)
                ct(:,k+1) = -sn*temp+cs*ct(:,k+1)
                if (k > lo) t(k+1,k-1) = (0.0_dp,0.0_dp)
            end do
            steps = steps+1
            if (.not. finite_complex_matrix(t)) return
        end do
        ok = finite_complex_matrix(bt) .and. finite_complex_matrix(ct)
    end subroutine schur_from_hessenberg

    subroutine triangular_solve_multiple(t,b,z,x,ok)
        ! Back substitutions in z*I-T for already prepared complex Schur T.
        complex(dp), intent(in) :: t(:,:) !! upper triangular state matrix
        complex(dp), intent(in) :: b(:,:) !! right hand sides
        complex(dp), intent(in) :: z !! Laplace/z-domain evaluation point
        complex(dp), intent(out) :: x(:,:) !! response states
        logical, intent(out) :: ok !! true if all pivots and results are finite
        complex(dp) :: pivot
        integer :: k, j, ns

        ok = .false.
        ns = size(t,1)
        x = (0.0_dp,0.0_dp)
        do k = ns,1,-1
            pivot = z-t(k,k)
            if (abs(pivot) <= tiny(1.0_dp)) return
            x(k,:) = b(k,:)
            do j = k+1,ns
                x(k,:) = x(k,:)+t(k,j)*x(j,:)
            end do
            x(k,:) = x(k,:)/pivot
        end do
        ok = finite_complex_matrix(x)
    end subroutine triangular_solve_multiple

    subroutine select_mimo_channels(a,b,c,d,inputs,outputs,bs,cs,ds,ok)
        ! Channels are one-based and may be reordered or repeated.
        complex(dp), intent(in) :: a(:,:) !! state matrix
        complex(dp), intent(in) :: b(:,:) !! input matrix
        complex(dp), intent(in) :: c(:,:) !! output matrix
        complex(dp), intent(in) :: d(:,:) !! feedthrough matrix
        integer, intent(in), optional :: inputs(:) !! selected input indices
        integer, intent(in), optional :: outputs(:) !! selected output indices
        complex(dp), allocatable, intent(out) :: bs(:,:) !! chosen input columns
        complex(dp), allocatable, intent(out) :: cs(:,:) !! chosen output rows
        complex(dp), allocatable, intent(out) :: ds(:,:) !! chosen direct terms
        logical, intent(out) :: ok !! true if dimensions and indices are valid
        integer, allocatable :: ii(:), oo(:)
        integer :: i, ns, ni, no

        ok = .false.
        ns = size(a,1)
        ni = size(b,2)
        no = size(c,1)
        if (ni < 1 .or. no < 1) return
        if (size(a,2) /= ns .or. size(b,1) /= ns) return
        if (size(c,2) /= ns .or. size(d,1) /= no .or. size(d,2) /= ni) return
        if (.not. finite_complex_matrix(a) .or. .not. finite_complex_matrix(b)) return
        if (.not. finite_complex_matrix(c) .or. .not. finite_complex_matrix(d)) return
        ii = [(i,i=1,ni)]
        oo = [(i,i=1,no)]
        if (present(inputs)) ii = inputs
        if (present(outputs)) oo = outputs
        if (size(ii) == 0 .or. size(oo) == 0) return
        if (any(ii < 1) .or. any(ii > ni)) return
        if (any(oo < 1) .or. any(oo > no)) return
        bs = b(:,ii)
        cs = c(oo,:)
        ds = d(oo,ii)
        ok = .true.
    end subroutine select_mimo_channels

    function prepare_mimo_frequency_complex(a,b,c,d,solver) result(plan)
        ! Prepare Q**H*A*Q once, reusable for analog/digital and Bode grids.
        complex(dp), intent(in) :: a(:,:) !! state matrix
        complex(dp), intent(in) :: b(:,:) !! input matrix
        complex(dp), intent(in) :: c(:,:) !! output matrix
        complex(dp), intent(in) :: d(:,:) !! direct feedthrough matrix
        character(len=*), intent(in), optional :: solver !! hessenberg (default) or schur
        type(mimo_frequency_plan) :: plan
        complex(dp), allocatable :: bs(:,:), cs(:,:), ds(:,:)
        complex(dp), allocatable :: t(:,:), bt(:,:), ct(:,:)
        real(dp), allocatable :: grid(:)
        logical :: ok

        if (present(solver)) then
            if (solver /= 'hessenberg' .and. solver /= 'schur') return
        end if
        call select_mimo_channels(a,b,c,d,bs=bs,cs=cs,ds=ds,ok=ok)
        if (.not. ok) return
        call reduce_to_hessenberg(a,bs,cs,plan%h,plan%br,plan%cr,ok)
        if (.not. ok) return
        if (present(solver)) then
            if (solver == 'schur') then
                call schur_from_hessenberg(plan%h,plan%br,plan%cr,t,bt,ct,ok)
                if (ok) then
                    call move_alloc(t,plan%h)
                    call move_alloc(bt,plan%br)
                    call move_alloc(ct,plan%cr)
                    plan%method = 'schur'
                end if
            end if
        end if
        plan%d = ds
        call state_space_analog_grid(a,1,grid,plan%analog_log_range)
        plan%status = design_status_success
    end function prepare_mimo_frequency_complex

    function prepare_mimo_frequency_real(a,b,c,d,solver) result(plan)
        real(dp), intent(in) :: a(:,:) !! state matrix
        real(dp), intent(in) :: b(:,:) !! input matrix
        real(dp), intent(in) :: c(:,:) !! output matrix
        real(dp), intent(in) :: d(:,:) !! feedthrough matrix
        character(len=*), intent(in), optional :: solver !! hessenberg (default) or schur
        type(mimo_frequency_plan) :: plan
        plan = prepare_mimo_frequency_complex(cmplx(a,0.0_dp,dp),cmplx(b,0.0_dp,dp), &
            cmplx(c,0.0_dp,dp),cmplx(d,0.0_dp,dp),solver)
    end function prepare_mimo_frequency_real

    function state_space_response_mimo(a,b,c,d,w,n,discrete,whole,dt,solver,inputs,outputs,pre_reduced) result(out)
        ! Evaluate C*(z*I-A)**(-1)*B+D at each frequency (all I/O pairs).
        complex(dp), intent(in) :: a(:,:) !! state matrix (nstate,nstate)
        complex(dp), intent(in) :: b(:,:) !! input matrix (nstate,ninput)
        complex(dp), intent(in) :: c(:,:) !! output matrix (noutput,nstate)
        complex(dp), intent(in) :: d(:,:) !! direct terms (noutput,ninput)
        real(dp), intent(in), optional :: w(:) !! explicit rad/s or rad/sample
        integer, intent(in), optional :: n !! auto frequency count (default 10000)
        logical, intent(in) :: discrete !! digital if true
        logical, intent(in), optional :: whole !! digital full unit circle
        real(dp), intent(in), optional :: dt !! positive sample interval
        character(len=*), intent(in), optional :: solver !! auto, direct, hessenberg, schur
        integer, intent(in), optional :: inputs(:) !! one-based input channels
        integer, intent(in), optional :: outputs(:) !! one-based output channels
        logical, intent(in), optional :: pre_reduced !! state matrix is already reduced
        type(mimo_frequency_response) :: out
        complex(dp), allocatable :: x(:,:), h(:,:), br(:,:), cr(:,:)
        complex(dp), allocatable :: bs(:,:), cs(:,:), ds(:,:)
        complex(dp) :: z
        real(dp) :: nan
        integer :: ns, ni, no, count, i
        logical :: ok, full_circle, use_reduction, reduced_in
        character(len=:), allocatable :: strategy

        ns = size(a,1)
        call select_mimo_channels(a,b,c,d,inputs,outputs,bs,cs,ds,ok)
        if (.not. ok) return
        ni = size(bs,2)
        no = size(cs,1)
        if (discrete .and. .not. valid_dt(dt)) return
        reduced_in = .false.
        if (present(pre_reduced)) reduced_in = pre_reduced
        strategy = 'auto'
        if (present(solver)) strategy = solver
        select case(strategy)
        case('auto','direct','hessenberg','schur')
            continue
        case default
            return
        end select
        if (present(w)) then
            if (.not. all(ieee_is_finite(w))) return
            out%frequencies = w
        else
            count = 10000
            if (present(n)) count = n
            if (count < 1) return
            if (discrete) then
                full_circle = .false.
                if (present(whole)) full_circle = whole
                allocate(out%frequencies(count))
                do i = 1,count
                    if (full_circle) then
                        out%frequencies(i) = 2.0_dp*pi*real(i-1,dp)/real(count,dp)
                    else
                        out%frequencies(i) = pi*real(i-1,dp)/real(count,dp)
                    end if
                end do
            else
                call state_space_analog_grid(a,count,out%frequencies)
            end if
        end if
        select case(strategy)
        case('auto')
            use_reduction = ns >= 8 .and. size(out%frequencies) >= 8
        case('direct')
            use_reduction = .false.
        case('hessenberg','schur')
            use_reduction = .true.
        case default
            out%status = design_status_invalid_input
            return
        end select
        if (reduced_in) use_reduction = .true.
        if (use_reduction .and. ns > 0 .and. .not. reduced_in) then
            call reduce_to_hessenberg(a,bs,cs,h,br,cr,ok)
            if (.not. ok) then
                out%status = design_status_invalid_input
                return
            end if
            if (strategy == 'schur') then
                block
                    complex(dp), allocatable :: t(:,:), bt(:,:), ct(:,:)
                    call schur_from_hessenberg(h,br,cr,t,bt,ct,ok)
                    if (ok) then
                        call move_alloc(t,h)
                        call move_alloc(bt,br)
                        call move_alloc(ct,cr)
                    else
                        ! Retain the original Hessenberg reduction on nonconvergence.
                        strategy = 'hessenberg'
                    end if
                end block
            end if
        end if
        allocate(out%response(no,ni,size(out%frequencies)),x(ns,ni))
        out%status = design_status_success
        nan = ieee_value(0.0_dp,ieee_quiet_nan)
        do i = 1,size(out%frequencies)
            z = cmplx(0.0_dp,out%frequencies(i),dp)
            if (discrete) z = exp(z)
            if (ns == 0) then
                out%response(:,:,i) = ds
            else
                if (reduced_in) then
                    if (strategy == 'schur') then
                        call triangular_solve_multiple(a,bs,z,x,ok)
                    else
                        call hessenberg_solve_multiple(a,bs,z,x,ok)
                    end if
                else if (use_reduction) then
                    if (strategy == 'schur') then
                        call triangular_solve_multiple(h,br,z,x,ok)
                    else
                        call hessenberg_solve_multiple(h,br,z,x,ok)
                    end if
                else
                    call shifted_solve_multiple(a,bs,z,x,ok)
                end if
                if (ok) then
                    if (reduced_in) then
                        out%response(:,:,i) = matmul(cs,x)+ds
                    else if (use_reduction) then
                        out%response(:,:,i) = matmul(cr,x)+ds
                    else
                        out%response(:,:,i) = matmul(cs,x)+ds
                    end if
                else
                    out%response(:,:,i) = cmplx(nan,nan,dp)
                    out%status = design_status_singular
                    cycle
                end if
            end if
            if (.not. finite_complex_matrix(out%response(:,:,i))) then
                out%response(:,:,i) = cmplx(nan,nan,dp)
                out%status = design_status_singular
            end if
        end do
    end function state_space_response_mimo

    function freqresp_mimo_real(a,b,c,d,w,n,solver,inputs,outputs) result(out)
        real(dp), intent(in) :: a(:,:) !! state matrix
        real(dp), intent(in) :: b(:,:) !! input matrix
        real(dp), intent(in) :: c(:,:) !! output matrix
        real(dp), intent(in) :: d(:,:) !! feedthrough matrix
        real(dp), intent(in), optional :: w(:) !! angular frequency, rad/s
        integer, intent(in), optional :: n !! frequency count
        character(len=*), intent(in), optional :: solver !! auto, direct, hessenberg, schur
        integer, intent(in), optional :: inputs(:) !! one-based input channel indices
        integer, intent(in), optional :: outputs(:) !! one-based output channel indices
        type(mimo_frequency_response) :: out
        out = state_space_response_mimo(cmplx(a,0.0_dp,dp),cmplx(b,0.0_dp,dp), &
            cmplx(c,0.0_dp,dp),cmplx(d,0.0_dp,dp),w,n,.false., &
            solver=solver,inputs=inputs,outputs=outputs)
    end function freqresp_mimo_real

    function freqresp_mimo_complex(a,b,c,d,w,n,solver,inputs,outputs) result(out)
        complex(dp), intent(in) :: a(:,:) !! state matrix
        complex(dp), intent(in) :: b(:,:) !! input matrix
        complex(dp), intent(in) :: c(:,:) !! output matrix
        complex(dp), intent(in) :: d(:,:) !! feedthrough matrix
        real(dp), intent(in), optional :: w(:) !! angular frequency, rad/s
        integer, intent(in), optional :: n !! frequency count
        character(len=*), intent(in), optional :: solver !! auto, direct, hessenberg, schur
        integer, intent(in), optional :: inputs(:) !! one-based input channel indices
        integer, intent(in), optional :: outputs(:) !! one-based output channel indices
        type(mimo_frequency_response) :: out
        out = state_space_response_mimo(a,b,c,d,w,n,.false.,solver=solver,inputs=inputs,outputs=outputs)
    end function freqresp_mimo_complex

    function dfreqresp_mimo_real(a,b,c,d,w,n,whole,dt,solver,inputs,outputs) result(out)
        real(dp), intent(in) :: a(:,:) !! discrete state matrix
        real(dp), intent(in) :: b(:,:) !! input matrix
        real(dp), intent(in) :: c(:,:) !! output matrix
        real(dp), intent(in) :: d(:,:) !! feedthrough matrix
        real(dp), intent(in), optional :: w(:) !! radians per sample
        integer, intent(in), optional :: n !! frequency count
        logical, intent(in), optional :: whole !! full-circle grid
        real(dp), intent(in), optional :: dt !! positive sample interval
        character(len=*), intent(in), optional :: solver !! auto, direct, hessenberg, schur
        integer, intent(in), optional :: inputs(:) !! one-based input channel indices
        integer, intent(in), optional :: outputs(:) !! one-based output channel indices
        type(mimo_frequency_response) :: out
        out = state_space_response_mimo(cmplx(a,0.0_dp,dp),cmplx(b,0.0_dp,dp), &
            cmplx(c,0.0_dp,dp),cmplx(d,0.0_dp,dp),w,n,.true.,whole,dt,solver,inputs,outputs)
    end function dfreqresp_mimo_real

    function dfreqresp_mimo_complex(a,b,c,d,w,n,whole,dt,solver,inputs,outputs) result(out)
        complex(dp), intent(in) :: a(:,:) !! discrete state matrix
        complex(dp), intent(in) :: b(:,:) !! input matrix
        complex(dp), intent(in) :: c(:,:) !! output matrix
        complex(dp), intent(in) :: d(:,:) !! feedthrough matrix
        real(dp), intent(in), optional :: w(:) !! radians per sample
        integer, intent(in), optional :: n !! frequency count
        logical, intent(in), optional :: whole !! full-circle grid
        real(dp), intent(in), optional :: dt !! positive sample interval
        character(len=*), intent(in), optional :: solver !! auto, direct, hessenberg, schur
        integer, intent(in), optional :: inputs(:) !! one-based input channel indices
        integer, intent(in), optional :: outputs(:) !! one-based output channel indices
        type(mimo_frequency_response) :: out
        out = state_space_response_mimo(a,b,c,d,w,n,.true.,whole,dt,solver,inputs,outputs)
    end function dfreqresp_mimo_complex

    function freqresp_mimo_plan(plan,w,n,inputs,outputs) result(out)
        type(mimo_frequency_plan), intent(in) :: plan !! precomputed reduction
        real(dp), intent(in), optional :: w(:) !! angular frequencies, rad/s
        integer, intent(in), optional :: n !! automatic frequency count
        integer, intent(in), optional :: inputs(:) !! one-based input channels
        integer, intent(in), optional :: outputs(:) !! one-based output channels
        type(mimo_frequency_response) :: out
        real(dp), allocatable :: auto_grid(:)
        integer :: count
        if (plan%status /= design_status_success) return
        if (present(w)) then
            out = state_space_response_mimo(plan%h,plan%br,plan%cr,plan%d,w,n,.false., &
                solver=plan%method,inputs=inputs,outputs=outputs,pre_reduced=.true.)
        else
            count = 10000
            if (present(n)) count = n
            if (count < 1) return
            call fill_mimo_log_grid(plan%analog_log_range,count,auto_grid)
            out = state_space_response_mimo(plan%h,plan%br,plan%cr,plan%d,w=auto_grid, &
                discrete=.false.,solver=plan%method,inputs=inputs,outputs=outputs, &
                pre_reduced=.true.)
        end if
    end function freqresp_mimo_plan

    function dfreqresp_mimo_plan(plan,w,n,whole,dt,inputs,outputs) result(out)
        type(mimo_frequency_plan), intent(in) :: plan !! precomputed reduction
        real(dp), intent(in), optional :: w(:) !! angular frequencies, rad/sample
        integer, intent(in), optional :: n !! automatic frequency count
        logical, intent(in), optional :: whole !! full-circle frequency grid
        real(dp), intent(in), optional :: dt !! positive sample interval
        integer, intent(in), optional :: inputs(:) !! one-based input channels
        integer, intent(in), optional :: outputs(:) !! one-based output channels
        type(mimo_frequency_response) :: out
        if (plan%status /= design_status_success) return
        out = state_space_response_mimo(plan%h,plan%br,plan%cr,plan%d,w,n,.true.,whole,dt, &
            plan%method,inputs,outputs,pre_reduced=.true.)
    end function dfreqresp_mimo_plan

    function bode_from_mimo(response,dt) result(out)
        ! Unwrap phase independently for each input/output channel.
        type(mimo_frequency_response), intent(in) :: response !! complex responses
        real(dp), intent(in), optional :: dt !! scale digital axis to rad/s
        type(mimo_bode_result) :: out
        real(dp) :: nan, scale, current, previous, theta, delta, corrected, amplitude
        integer :: no, ni, i, j, k
        logical :: have_previous

        out%status = response%status
        if (response%status == design_status_invalid_input) return
        if (.not. allocated(response%frequencies)) return
        if (.not. allocated(response%response)) return
        if (.not. valid_dt(dt)) then
            out%status = design_status_invalid_input
            return
        end if
        scale = 1.0_dp
        if (present(dt)) scale = dt
        out%frequencies = response%frequencies/scale
        no = size(response%response,1)
        ni = size(response%response,2)
        allocate(out%magnitude(no,ni,size(out%frequencies)))
        allocate(out%phase(no,ni,size(out%frequencies)))
        nan = ieee_value(0.0_dp,ieee_quiet_nan)
        do j = 1,ni
            do i = 1,no
                theta = 0.0_dp
                previous = 0.0_dp
                have_previous = .false.
                do k = 1,size(out%frequencies)
                    if (.not. ieee_is_finite(real(response%response(i,j,k),dp)) .or. &
                        .not. ieee_is_finite(aimag(response%response(i,j,k)))) then
                        out%magnitude(i,j,k) = nan
                        out%phase(i,j,k) = nan
                        have_previous = .false.
                        cycle
                    end if
                    current = atan2(aimag(response%response(i,j,k)), &
                        real(response%response(i,j,k),dp))
                    if (.not. have_previous) then
                        theta = current
                    else
                        delta = current-previous
                        corrected = modulo(delta+pi,2.0_dp*pi)-pi
                        if (corrected == -pi .and. delta > 0.0_dp) corrected = pi
                        theta = theta+corrected
                    end if
                    out%phase(i,j,k) = theta*180.0_dp/pi
                    amplitude = abs(response%response(i,j,k))
                    if (amplitude == 0.0_dp) then
                        out%magnitude(i,j,k) = ieee_value(0.0_dp,ieee_negative_inf)
                    else
                        out%magnitude(i,j,k) = 20.0_dp*log10(amplitude)
                    end if
                    previous = current
                    have_previous = .true.
                end do
            end do
        end do
    end function bode_from_mimo

    function bode_mimo_real(a,b,c,d,w,n,solver,inputs,outputs) result(out)
        real(dp), intent(in) :: a(:,:) !! continuous state matrix
        real(dp), intent(in) :: b(:,:) !! input matrix
        real(dp), intent(in) :: c(:,:) !! output matrix
        real(dp), intent(in) :: d(:,:) !! feedthrough matrix
        real(dp), intent(in), optional :: w(:) !! angular frequencies, rad/s
        integer, intent(in), optional :: n !! sample count, default 100
        character(len=*), intent(in), optional :: solver !! auto, direct, hessenberg, schur
        integer, intent(in), optional :: inputs(:) !! one-based input channel indices
        integer, intent(in), optional :: outputs(:) !! one-based output channel indices
        type(mimo_bode_result) :: out
        type(mimo_frequency_response) :: response
        integer :: count
        count = 100
        if (present(n)) count = n
        response = freqresp_mimo(a,b,c,d,w=w,n=count,solver=solver,inputs=inputs,outputs=outputs)
        out = bode_from_mimo(response)
    end function bode_mimo_real

    function bode_mimo_complex(a,b,c,d,w,n,solver,inputs,outputs) result(out)
        complex(dp), intent(in) :: a(:,:) !! continuous state matrix
        complex(dp), intent(in) :: b(:,:) !! input matrix
        complex(dp), intent(in) :: c(:,:) !! output matrix
        complex(dp), intent(in) :: d(:,:) !! feedthrough matrix
        real(dp), intent(in), optional :: w(:) !! angular frequencies, rad/s
        integer, intent(in), optional :: n !! sample count, default 100
        character(len=*), intent(in), optional :: solver !! auto, direct, hessenberg, schur
        integer, intent(in), optional :: inputs(:) !! one-based input channel indices
        integer, intent(in), optional :: outputs(:) !! one-based output channel indices
        type(mimo_bode_result) :: out
        type(mimo_frequency_response) :: response
        integer :: count
        count = 100
        if (present(n)) count = n
        response = freqresp_mimo(a,b,c,d,w=w,n=count,solver=solver,inputs=inputs,outputs=outputs)
        out = bode_from_mimo(response)
    end function bode_mimo_complex

    function dbode_mimo_real(a,b,c,d,w,n,dt,solver,inputs,outputs) result(out)
        real(dp), intent(in) :: a(:,:) !! discrete state matrix
        real(dp), intent(in) :: b(:,:) !! input matrix
        real(dp), intent(in) :: c(:,:) !! output matrix
        real(dp), intent(in) :: d(:,:) !! feedthrough matrix
        real(dp), intent(in), optional :: w(:) !! radians per sample
        integer, intent(in), optional :: n !! sample count, default 100
        real(dp), intent(in), optional :: dt !! sample interval, default one
        character(len=*), intent(in), optional :: solver !! auto, direct, hessenberg, schur
        integer, intent(in), optional :: inputs(:) !! one-based input channel indices
        integer, intent(in), optional :: outputs(:) !! one-based output channel indices
        type(mimo_bode_result) :: out
        type(mimo_frequency_response) :: response
        integer :: count
        count = 100
        if (present(n)) count = n
        response = dfreqresp_mimo(a,b,c,d,w=w,n=count,dt=dt,solver=solver,inputs=inputs,outputs=outputs)
        out = bode_from_mimo(response,dt)
    end function dbode_mimo_real

    function dbode_mimo_complex(a,b,c,d,w,n,dt,solver,inputs,outputs) result(out)
        complex(dp), intent(in) :: a(:,:) !! discrete state matrix
        complex(dp), intent(in) :: b(:,:) !! input matrix
        complex(dp), intent(in) :: c(:,:) !! output matrix
        complex(dp), intent(in) :: d(:,:) !! feedthrough matrix
        real(dp), intent(in), optional :: w(:) !! radians per sample
        integer, intent(in), optional :: n !! sample count, default 100
        real(dp), intent(in), optional :: dt !! sample interval, default one
        character(len=*), intent(in), optional :: solver !! auto, direct, hessenberg, schur
        integer, intent(in), optional :: inputs(:) !! one-based input channel indices
        integer, intent(in), optional :: outputs(:) !! one-based output channel indices
        type(mimo_bode_result) :: out
        type(mimo_frequency_response) :: response
        integer :: count
        count = 100
        if (present(n)) count = n
        response = dfreqresp_mimo(a,b,c,d,w=w,n=count,dt=dt,solver=solver,inputs=inputs,outputs=outputs)
        out = bode_from_mimo(response,dt)
    end function dbode_mimo_complex

    function bode_mimo_plan(plan,w,n,inputs,outputs) result(out)
        type(mimo_frequency_plan), intent(in) :: plan !! precomputed reduction
        real(dp), intent(in), optional :: w(:) !! angular frequency grid
        integer, intent(in), optional :: n !! number of frequency samples
        integer, intent(in), optional :: inputs(:) !! one-based input channels
        integer, intent(in), optional :: outputs(:) !! one-based output channels
        type(mimo_bode_result) :: out
        type(mimo_frequency_response) :: response
        integer :: count
        count = 100
        if (present(n)) count = n
        response = freqresp_mimo(plan,w=w,n=count,inputs=inputs,outputs=outputs)
        out = bode_from_mimo(response)
    end function bode_mimo_plan

    function dbode_mimo_plan(plan,w,n,dt,inputs,outputs) result(out)
        type(mimo_frequency_plan), intent(in) :: plan !! precomputed reduction
        real(dp), intent(in), optional :: w(:) !! angular frequency grid
        integer, intent(in), optional :: n !! number of frequency samples
        real(dp), intent(in), optional :: dt !! sample interval
        integer, intent(in), optional :: inputs(:) !! one-based input channels
        integer, intent(in), optional :: outputs(:) !! one-based output channels
        type(mimo_bode_result) :: out
        type(mimo_frequency_response) :: response
        integer :: count
        count = 100
        if (present(n)) count = n
        response = dfreqresp_mimo(plan,w=w,n=count,dt=dt,inputs=inputs,outputs=outputs)
        out = bode_from_mimo(response,dt)
    end function dbode_mimo_plan

end module scifort_signal_lti_frequency
