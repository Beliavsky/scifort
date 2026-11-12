! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Actions of the exponential of complex sparse matrices, exp(t*A)*B.
! Scaling and bounded Taylor expansion, with an optional trace shift.
! The operator interface requires a certified/conservative infinity-norm
! bound for A. It never attempts to estimate that bound from opaque callbacks.
module scifort_sparse_linalg_complex_expm
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_sparse_complex, only : csr_complex_matrix, csc_complex_matrix, &
        valid_csr_complex, valid_csc_complex, csc_complex_to_csr
    implicit none
    private

    integer, parameter, public :: complex_expm_success = 0
    integer, parameter, public :: complex_expm_no_convergence = 1
    integer, parameter, public :: complex_expm_invalid_input = -1
    integer, parameter, public :: complex_expm_breakdown = -2
    integer, parameter, public :: complex_expm_callback_failure = -3

    public :: expm_multiply_complex, expm_multiply_complex_operator
    public :: expm_multiply_complex_interval, expm_multiply_complex_operator_interval
    public :: complex_expm_matvec

    abstract interface
        subroutine complex_expm_matvec(x, y, info)
            import dp
            complex(dp), intent(in) :: x(:) !! Input vector, length n.
            complex(dp), intent(out) :: y(:) !! A*x, length n.
            integer, intent(out) :: info !! Zero on success; nonzero otherwise.
        end subroutine complex_expm_matvec
    end interface

    interface expm_multiply_complex
        module procedure expm_csr_vector
        module procedure expm_csr_block
        module procedure expm_csc_vector
        module procedure expm_csc_block
    end interface expm_multiply_complex

    interface expm_multiply_complex_operator
        module procedure expm_operator_vector
        module procedure expm_operator_block
    end interface expm_multiply_complex_operator

    interface expm_multiply_complex_interval
        module procedure interval_csr_vector
        module procedure interval_csr_block
        module procedure interval_csc_vector
        module procedure interval_csc_block
    end interface expm_multiply_complex_interval

    interface expm_multiply_complex_operator_interval
        module procedure interval_operator_vector
        module procedure interval_operator_block
    end interface expm_multiply_complex_operator_interval

    abstract interface
        subroutine complex_expm_step(time, b, y, info, products)
            import dp
            real(dp), intent(in) :: time !! Requested time increment.
            complex(dp), intent(in) :: b(:,:) !! Dense input block.
            complex(dp), allocatable, intent(out) :: y(:,:) !! Propagated block.
            integer, intent(out) :: info !! Zero on success.
            integer, intent(out) :: products !! Number of operator actions.
        end subroutine complex_expm_step
    end interface

contains

    pure elemental logical function finite_complex(z) result(ok)
        complex(dp), intent(in) :: z !! Candidate number.
        ok = ieee_is_finite(real(z,dp)) .and. ieee_is_finite(aimag(z))
    end function finite_complex

    subroutine shifted_norm(a, tracea, mu, bound, info)
        type(csr_complex_matrix), intent(in) :: a !! Valid square canonical CSR matrix.
        complex(dp), intent(in), optional :: tracea !! Optional trace-shift hint.
        complex(dp), intent(out) :: mu !! Chosen diagonal shift.
        real(dp), intent(out) :: bound !! Exact infinity norm of A-mu*I.
        integer, intent(out) :: info !! Zero on success.
        complex(dp) :: actual_mu, suggested_mu
        real(dp) :: actual_bound, suggested_bound
        integer :: i, k

        info = complex_expm_success
        mu = cmplx(0.0_dp,0.0_dp,dp)
        bound = 0.0_dp
        if (a%n_rows == 0) return
        actual_mu = cmplx(0.0_dp,0.0_dp,dp)
        do i = 1,a%n_rows
            do k = a%indptr(i),a%indptr(i+1)-1
                if (a%indices(k) == i) actual_mu = actual_mu + a%data(k)/real(a%n_rows,dp)
            end do
        end do
        if (.not. finite_complex(actual_mu)) then
            info = complex_expm_breakdown
            return
        end if
        call norm_with_shift(a,actual_mu,actual_bound,info)
        if (info /= complex_expm_success) return
        mu = actual_mu
        bound = actual_bound
        if (present(tracea)) then
            suggested_mu = tracea/real(a%n_rows,dp)
            if (.not. finite_complex(suggested_mu)) return
            call norm_with_shift(a,suggested_mu,suggested_bound,info)
            if (info /= complex_expm_success) return
            if (suggested_bound < bound) then
                mu = suggested_mu
                bound = suggested_bound
            end if
        end if
    end subroutine shifted_norm

    subroutine norm_with_shift(a, mu, bound, info)
        type(csr_complex_matrix), intent(in) :: a !! Valid square CSR input.
        complex(dp), intent(in) :: mu !! Diagonal shift of A.
        real(dp), intent(out) :: bound !! Maximum absolute row sum.
        integer, intent(out) :: info !! Zero if finite.
        real(dp) :: row_sum
        logical :: has_diagonal
        integer :: i, k, j

        bound = 0.0_dp
        info = complex_expm_success
        do i = 1,a%n_rows
            row_sum = 0.0_dp
            has_diagonal = .false.
            do k = a%indptr(i),a%indptr(i+1)-1
                j = a%indices(k)
                if (j == i) then
                    row_sum = row_sum + abs(a%data(k)-mu)
                    has_diagonal = .true.
                else
                    row_sum = row_sum + abs(a%data(k))
                end if
            end do
            if (.not. has_diagonal) row_sum = row_sum + abs(mu)
            if (.not. ieee_is_finite(row_sum)) then
                info = complex_expm_breakdown
                return
            end if
            bound = max(bound,row_sum)
        end do
    end subroutine norm_with_shift

    subroutine expm_csr_block(a, b, y, info, t, tracea, tol, scaling_steps, products)
        type(csr_complex_matrix), intent(in) :: a !! Square canonical CSR matrix.
        complex(dp), intent(in) :: b(:,:) !! Finite RHS block (n,nrhs).
        complex(dp), allocatable, intent(out) :: y(:,:) !! exp(t*A)*b, same shape as b.
        integer, intent(out) :: info !! 0 success, 1 Taylor exhaustion, negative error.
        real(dp), intent(in), optional :: t !! Real time, default 1.
        complex(dp), intent(in), optional :: tracea !! Optional trace hint.
        real(dp), intent(in), optional :: tol !! Relative Taylor remainder target.
        integer, intent(out), optional :: scaling_steps !! Number of time substeps.
        integer, intent(out), optional :: products !! Number of matrix-vector actions.
        complex(dp) :: mu
        real(dp) :: bound
        integer :: status

        allocate(y(0,0))
        info = complex_expm_invalid_input
        if (present(scaling_steps)) scaling_steps = 0
        if (present(products)) products = 0
        if (.not. valid_csr_complex(a)) return
        if (a%n_rows /= a%n_cols .or. size(b,1) /= a%n_rows) return
        if (.not. all(finite_complex(b))) return
        if (present(tracea)) then
            if (.not. finite_complex(tracea)) return
        end if
        call shifted_norm(a,tracea,mu,bound,status)
        if (status /= complex_expm_success) then
            info = status
            return
        end if
        call expm_core(a%n_rows,b,csr_action,mu,bound,y,info,t,tol,scaling_steps,products)

    contains
        subroutine csr_action(x, result, ierr)
            complex(dp), intent(in) :: x(:) !! Input vector.
            complex(dp), intent(out) :: result(:) !! CSR times x.
            integer, intent(out) :: ierr !! Zero on success.
            integer :: i, k

            ierr = 0
            result = cmplx(0.0_dp,0.0_dp,dp)
            do i = 1,a%n_rows
                do k = a%indptr(i),a%indptr(i+1)-1
                    result(i) = result(i) + a%data(k)*x(a%indices(k))
                end do
            end do
        end subroutine csr_action
    end subroutine expm_csr_block

    subroutine expm_csr_vector(a, b, y, info, t, tracea, tol, scaling_steps, products)
        type(csr_complex_matrix), intent(in) :: a !! Square canonical CSR matrix.
        complex(dp), intent(in) :: b(:) !! Finite vector of length n.
        complex(dp), allocatable, intent(out) :: y(:) !! exp(t*A)*b.
        integer, intent(out) :: info !! Status code.
        real(dp), intent(in), optional :: t !! Real time, default 1.
        complex(dp), intent(in), optional :: tracea !! Optional trace hint.
        real(dp), intent(in), optional :: tol !! Relative Taylor remainder target.
        integer, intent(out), optional :: scaling_steps !! Number of time substeps.
        integer, intent(out), optional :: products !! Number of matrix-vector actions.
        complex(dp), allocatable :: block(:,:), result(:,:)

        allocate(block(size(b),1))
        block(:,1) = b
        call expm_csr_block(a,block,result,info,t,tracea,tol,scaling_steps,products)
        if (size(result,2) == 1) then
            y = result(:,1)
        else
            allocate(y(0))
        end if
    end subroutine expm_csr_vector

    subroutine expm_csc_block(a, b, y, info, t, tracea, tol, scaling_steps, products)
        type(csc_complex_matrix), intent(in) :: a !! Square canonical CSC matrix.
        complex(dp), intent(in) :: b(:,:) !! Finite RHS block (n,nrhs).
        complex(dp), allocatable, intent(out) :: y(:,:) !! exp(t*A)*b.
        integer, intent(out) :: info !! Status code.
        real(dp), intent(in), optional :: t !! Real time, default 1.
        complex(dp), intent(in), optional :: tracea !! Optional trace hint.
        real(dp), intent(in), optional :: tol !! Relative Taylor remainder target.
        integer, intent(out), optional :: scaling_steps !! Number of substeps.
        integer, intent(out), optional :: products !! Matrix-vector action count.
        type(csr_complex_matrix) :: csr
        integer :: status

        allocate(y(0,0))
        info = complex_expm_invalid_input
        if (present(scaling_steps)) scaling_steps = 0
        if (present(products)) products = 0
        if (.not. valid_csc_complex(a)) return
        if (a%n_rows /= a%n_cols .or. size(b,1) /= a%n_rows) return
        csr = csc_complex_to_csr(a,status)
        if (status /= 0) return
        deallocate(y)
        call expm_csr_block(csr,b,y,info,t,tracea,tol,scaling_steps,products)
    end subroutine expm_csc_block

    subroutine expm_csc_vector(a, b, y, info, t, tracea, tol, scaling_steps, products)
        type(csc_complex_matrix), intent(in) :: a !! Square canonical CSC matrix.
        complex(dp), intent(in) :: b(:) !! Finite vector of length n.
        complex(dp), allocatable, intent(out) :: y(:) !! exp(t*A)*b.
        integer, intent(out) :: info !! Status code.
        real(dp), intent(in), optional :: t !! Real time, default 1.
        complex(dp), intent(in), optional :: tracea !! Optional trace hint.
        real(dp), intent(in), optional :: tol !! Relative Taylor remainder target.
        integer, intent(out), optional :: scaling_steps !! Number of substeps.
        integer, intent(out), optional :: products !! Matrix-vector action count.
        complex(dp), allocatable :: block(:,:), result(:,:)

        allocate(block(size(b),1))
        block(:,1) = b
        call expm_csc_block(a,block,result,info,t,tracea,tol,scaling_steps,products)
        if (size(result,2) == 1) then
            y = result(:,1)
        else
            allocate(y(0))
        end if
    end subroutine expm_csc_vector

    subroutine expm_operator_block(n, b, matvec, norm_bound, y, info, t, tracea, tol, &
                                    scaling_steps, products)
        integer, intent(in) :: n !! Operator dimension, nonnegative.
        complex(dp), intent(in) :: b(:,:) !! Dense input block (n,nrhs).
        procedure(complex_expm_matvec) :: matvec !! Callback computing A*x.
        real(dp), intent(in) :: norm_bound !! Conservative bound >= infinity norm of A.
        complex(dp), allocatable, intent(out) :: y(:,:) !! exp(t*A)*b.
        integer, intent(out) :: info !! Status code.
        real(dp), intent(in), optional :: t !! Real time, default 1.
        complex(dp), intent(in), optional :: tracea !! Optional supplied trace.
        real(dp), intent(in), optional :: tol !! Relative Taylor remainder target.
        integer, intent(out), optional :: scaling_steps !! Number of substeps.
        integer, intent(out), optional :: products !! Number of callback invocations.
        complex(dp) :: mu
        real(dp) :: bound

        allocate(y(0,0))
        info = complex_expm_invalid_input
        if (present(scaling_steps)) scaling_steps = 0
        if (present(products)) products = 0
        if (n < 0 .or. n == huge(n) .or. size(b,1) /= n) return
        if (.not. all(finite_complex(b))) return
        if (.not. ieee_is_finite(norm_bound) .or. norm_bound < 0.0_dp) return
        mu = cmplx(0.0_dp,0.0_dp,dp)
        if (present(tracea)) then
            if (.not. finite_complex(tracea)) return
            if (n > 0) mu = tracea/real(n,dp)
        end if
        if (norm_bound > huge(1.0_dp)-abs(mu)) return
        bound = norm_bound+abs(mu)
        deallocate(y)
        call expm_core(n,b,matvec,mu,bound,y,info,t,tol,scaling_steps,products)
    end subroutine expm_operator_block

    subroutine expm_operator_vector(n, b, matvec, norm_bound, y, info, t, tracea, tol, &
                                     scaling_steps, products)
        integer, intent(in) :: n !! Operator dimension.
        complex(dp), intent(in) :: b(:) !! Input vector.
        procedure(complex_expm_matvec) :: matvec !! Callback A*x.
        real(dp), intent(in) :: norm_bound !! Conservative bound on infinity norm of A.
        complex(dp), allocatable, intent(out) :: y(:) !! exp(t*A)*b.
        integer, intent(out) :: info !! Status code.
        real(dp), intent(in), optional :: t !! Real time, default 1.
        complex(dp), intent(in), optional :: tracea !! Optional trace of A.
        real(dp), intent(in), optional :: tol !! Relative Taylor remainder target.
        integer, intent(out), optional :: scaling_steps !! Number of substeps.
        integer, intent(out), optional :: products !! Callback invocation count.
        complex(dp), allocatable :: block(:,:), result(:,:)

        allocate(block(size(b),1))
        block(:,1) = b
        call expm_operator_block(n,block,matvec,norm_bound,result,info,t,tracea,tol, &
                                 scaling_steps,products)
        if (size(result,2) == 1) then
            y = result(:,1)
        else
            allocate(y(0))
        end if
    end subroutine expm_operator_vector

    subroutine expm_core(n, b, matvec, mu, norm_bound, y, info, t, tol, scaling_steps, products)
        integer, intent(in) :: n !! Operator dimension.
        complex(dp), intent(in) :: b(:,:) !! Dense RHS block (n,nrhs).
        procedure(complex_expm_matvec) :: matvec !! Callback A*x, status zero if successful.
        complex(dp), intent(in) :: mu !! Scalar diagonal shift.
        real(dp), intent(in) :: norm_bound !! Upper bound on norm_infinity(A-mu*I).
        complex(dp), allocatable, intent(out) :: y(:,:) !! Result block (n,nrhs).
        integer, intent(out) :: info !! Success, exhaustion, invalid, breakdown, callback error.
        real(dp), intent(in), optional :: t !! Scalar real time.
        real(dp), intent(in), optional :: tol !! Relative Taylor remainder threshold.
        integer, intent(out), optional :: scaling_steps !! Number of substeps.
        integer, intent(out), optional :: products !! Callback invocation count.
        complex(dp), allocatable :: term(:,:), accumulator(:,:), product(:,:)
        complex(dp) :: shift_factor
        real(dp) :: time, threshold, scaled_norm, h, beta, ratio, denominator, tail
        integer :: s, step, j, col, status, nprod
        logical :: converged
        integer, parameter :: max_taylor_terms = 80
        real(dp), parameter :: max_step_norm = 2.0_dp
        integer, parameter :: max_scaling_steps = 200000

        info = complex_expm_invalid_input
        if (present(scaling_steps)) scaling_steps = 0
        if (present(products)) products = 0
        allocate(y(0,0))
        if (n < 0 .or. size(b,1) /= n .or. norm_bound < 0.0_dp) return
        if (.not. ieee_is_finite(norm_bound) .or. .not. finite_complex(mu)) return
        if (.not. all(finite_complex(b))) return
        time = 1.0_dp
        if (present(t)) time = t
        threshold = 4.0e-15_dp
        if (present(tol)) threshold = tol
        if (.not. ieee_is_finite(time) .or. .not. ieee_is_finite(threshold)) return
        if (threshold < 0.0_dp .or. threshold >= 1.0_dp) return
        threshold = max(threshold,4.0_dp*epsilon(1.0_dp))
        deallocate(y)
        y = b
        info = complex_expm_success
        if (n == 0 .or. size(b,2) == 0 .or. time == 0.0_dp) return
        if (all(b == cmplx(0.0_dp,0.0_dp,dp))) return
        if (norm_bound == 0.0_dp) then
            shift_factor = exp(time*mu)
            if (.not. finite_complex(shift_factor)) then
                info = complex_expm_breakdown
                return
            end if
            y = shift_factor*b
            if (.not. all(finite_complex(y))) info = complex_expm_breakdown
            return
        end if
        if (abs(time) > huge(1.0_dp)/norm_bound) then
            info = complex_expm_breakdown
            return
        end if
        scaled_norm = abs(time)*norm_bound
        if (scaled_norm > real(max_scaling_steps,dp)*max_step_norm) then
            info = complex_expm_invalid_input
            return
        end if
        s = max(1,ceiling(scaled_norm/max_step_norm))
        h = time/real(s,dp)
        beta = abs(h)*norm_bound
        shift_factor = exp(h*mu)
        if (.not. finite_complex(shift_factor)) then
            info = complex_expm_breakdown
            return
        end if
        if (present(scaling_steps)) scaling_steps = s
        allocate(term(n,size(b,2)),accumulator(n,size(b,2)),product(n,size(b,2)))
        nprod = 0
        do step = 1,s
            term = y
            accumulator = y
            converged = .false.
            do j = 1,max_taylor_terms
                do col = 1,size(b,2)
                    call matvec(term(:,col),product(:,col),status)
                    nprod = nprod+1
                    if (status /= 0) then
                        info = complex_expm_callback_failure
                        exit
                    end if
                end do
                if (info /= complex_expm_success) exit
                if (.not. all(finite_complex(product))) then
                    info = complex_expm_breakdown
                    exit
                end if
                term = (h/real(j,dp))*(product-mu*term)
                accumulator = accumulator+term
                if (.not. all(finite_complex(accumulator)) .or. .not. all(finite_complex(term))) then
                    info = complex_expm_breakdown
                    exit
                end if
                ratio = beta/real(j+1,dp)
                denominator = 1.0_dp-beta/real(j+2,dp)
                if (denominator > 0.0_dp .and. ratio < 1.0_dp) then
                    converged = .true.
                    do col = 1,size(b,2)
                        tail = maxval(abs(term(:,col)))*ratio/denominator
                        if (tail > threshold*maxval(abs(accumulator(:,col)))) then
                            converged = .false.
                            exit
                        end if
                    end do
                    if (converged) exit
                end if
            end do
            if (info /= complex_expm_success) exit
            if (.not. converged) then
                info = complex_expm_no_convergence
                exit
            end if
            y = shift_factor*accumulator
            if (.not. all(finite_complex(y))) then
                info = complex_expm_breakdown
                exit
            end if
        end do
        if (present(products)) products = nprod
    end subroutine expm_core


    ! Propagate from the requested time closest to zero. This avoids losing
    ! information by propagating forward from a strongly decayed first point.
    subroutine interval_core(n, b, start, stop, num, step, y, info, endpoint, products)
        integer, intent(in) :: n !! Square operator dimension.
        complex(dp), intent(in) :: b(:,:) !! Input RHS block (n,nrhs).
        real(dp), intent(in) :: start !! First requested time.
        real(dp), intent(in) :: stop !! Last time bound.
        integer, intent(in) :: num !! Number of samples, at least one.
        procedure(complex_expm_step) :: step !! Apply exponential at one scalar time.
        complex(dp), allocatable, intent(out) :: y(:,:,:) !! (n,nrhs,num) result.
        integer, intent(out) :: info !! Status code.
        logical, intent(in), optional :: endpoint !! Include stop, default true.
        integer, intent(out), optional :: products !! Total matrix-vector actions.
        complex(dp), allocatable :: current(:,:), next(:,:)
        real(dp) :: delta, time_here, anchor_time
        integer :: i, anchor, denominator, nprod, total_products
        logical :: include_endpoint

        info = complex_expm_invalid_input
        if (present(products)) products = 0
        allocate(y(0,0,0))
        if (n < 0 .or. size(b,1) /= n) return
        if (.not. all(finite_complex(b))) return
        if (num < 1 .or. .not. ieee_is_finite(start) .or. .not. ieee_is_finite(stop)) return
        include_endpoint = .true.
        if (present(endpoint)) include_endpoint = endpoint
        delta = 0.0_dp
        if (num > 1) then
            denominator = num
            if (include_endpoint) denominator = num-1
            delta = (stop-start)/real(denominator,dp)
            if (.not. ieee_is_finite(delta)) return
        end if
        anchor = 1
        anchor_time = start
        do i = 2,num
            time_here = start+real(i-1,dp)*delta
            if (.not. ieee_is_finite(time_here)) return
            if (abs(time_here) < abs(anchor_time)) then
                anchor = i
                anchor_time = time_here
            end if
        end do
        call step(anchor_time,b,current,info,nprod)
        if (info /= complex_expm_success) return
        total_products = nprod
        deallocate(y)
        allocate(y(n,size(b,2),num))
        y = cmplx(0.0_dp,0.0_dp,dp)
        y(:,:,anchor) = current
        do i = anchor+1,num
            call step(delta,current,next,info,nprod)
            total_products = total_products+nprod
            if (info /= complex_expm_success) exit
            y(:,:,i) = next
            call move_alloc(next,current)
        end do
        if (info == complex_expm_success) then
            current = y(:,:,anchor)
            do i = anchor-1,1,-1
                call step(-delta,current,next,info,nprod)
                total_products = total_products+nprod
                if (info /= complex_expm_success) exit
                y(:,:,i) = next
                call move_alloc(next,current)
            end do
        end if
        if (present(products)) products = total_products
    end subroutine interval_core

    subroutine interval_csr_block(a, b, start, stop, num, y, info, endpoint, tracea, tol, products)
        type(csr_complex_matrix), intent(in) :: a !! Square canonical CSR input.
        complex(dp), intent(in) :: b(:,:) !! Finite input block.
        real(dp), intent(in) :: start !! First time.
        real(dp), intent(in) :: stop !! Last time bound.
        integer, intent(in) :: num !! Number of samples.
        complex(dp), allocatable, intent(out) :: y(:,:,:) !! Output (n,nrhs,num).
        integer, intent(out) :: info !! Status code.
        logical, intent(in), optional :: endpoint !! Include stop.
        complex(dp), intent(in), optional :: tracea !! Optional trace hint.
        real(dp), intent(in), optional :: tol !! Relative Taylor remainder tolerance.
        integer, intent(out), optional :: products !! Matrix-vector action count.

        allocate(y(0,0,0))
        info = complex_expm_invalid_input
        if (present(products)) products = 0
        if (.not. valid_csr_complex(a)) return
        if (a%n_cols /= a%n_rows .or. size(b,1) /= a%n_rows) return
        if (present(tracea)) then
            if (.not. finite_complex(tracea)) return
        end if
        deallocate(y)
        call interval_core(a%n_rows,b,start,stop,num,step,y,info,endpoint,products)
    contains
        subroutine step(time,x,z,status,nprod)
            real(dp), intent(in) :: time !! Requested time.
            complex(dp), intent(in) :: x(:,:) !! Input block.
            complex(dp), allocatable, intent(out) :: z(:,:) !! Output block.
            integer, intent(out) :: status !! Computation status.
            integer, intent(out) :: nprod !! Matrix-vector action count.
            call expm_csr_block(a,x,z,status,t=time,tracea=tracea,tol=tol,products=nprod)
        end subroutine step
    end subroutine interval_csr_block

    subroutine interval_csr_vector(a, b, start, stop, num, y, info, endpoint, tracea, tol, products)
        type(csr_complex_matrix), intent(in) :: a !! Square canonical CSR input.
        complex(dp), intent(in) :: b(:) !! Finite input vector.
        real(dp), intent(in) :: start !! First time.
        real(dp), intent(in) :: stop !! Last time bound.
        integer, intent(in) :: num !! Number of samples.
        complex(dp), allocatable, intent(out) :: y(:,:) !! Output (n,num).
        integer, intent(out) :: info !! Status code.
        logical, intent(in), optional :: endpoint !! Include stop.
        complex(dp), intent(in), optional :: tracea !! Optional trace hint.
        real(dp), intent(in), optional :: tol !! Relative Taylor remainder tolerance.
        integer, intent(out), optional :: products !! Matrix-vector action count.
        complex(dp), allocatable :: block(:,:), result(:,:,:)

        allocate(block(size(b),1))
        block(:,1) = b
        call interval_csr_block(a,block,start,stop,num,result,info,endpoint,tracea,tol,products)
        if (size(result,2) == 1) then
            y = result(:,1,:)
        else
            allocate(y(0,0))
        end if
    end subroutine interval_csr_vector

    subroutine interval_csc_block(a, b, start, stop, num, y, info, endpoint, tracea, tol, products)
        type(csc_complex_matrix), intent(in) :: a !! Square canonical CSC input.
        complex(dp), intent(in) :: b(:,:) !! Finite input block.
        real(dp), intent(in) :: start !! First time.
        real(dp), intent(in) :: stop !! Last time bound.
        integer, intent(in) :: num !! Number of samples.
        complex(dp), allocatable, intent(out) :: y(:,:,:) !! Output (n,nrhs,num).
        integer, intent(out) :: info !! Status code.
        logical, intent(in), optional :: endpoint !! Include stop.
        complex(dp), intent(in), optional :: tracea !! Optional trace hint.
        real(dp), intent(in), optional :: tol !! Relative Taylor remainder tolerance.
        integer, intent(out), optional :: products !! Matrix-vector action count.
        type(csr_complex_matrix) :: csr
        integer :: status

        info = complex_expm_invalid_input
        if (present(products)) products = 0
        allocate(y(0,0,0))
        if (.not. valid_csc_complex(a)) return
        if (a%n_rows /= a%n_cols .or. size(b,1) /= a%n_rows) return
        csr = csc_complex_to_csr(a,status)
        if (status /= 0) return
        deallocate(y)
        call interval_csr_block(csr,b,start,stop,num,y,info,endpoint,tracea,tol,products)
    end subroutine interval_csc_block

    subroutine interval_csc_vector(a, b, start, stop, num, y, info, endpoint, tracea, tol, products)
        type(csc_complex_matrix), intent(in) :: a !! Square canonical CSC input.
        complex(dp), intent(in) :: b(:) !! Finite input vector.
        real(dp), intent(in) :: start !! First time.
        real(dp), intent(in) :: stop !! Last time bound.
        integer, intent(in) :: num !! Number of samples.
        complex(dp), allocatable, intent(out) :: y(:,:) !! Output (n,num).
        integer, intent(out) :: info !! Status code.
        logical, intent(in), optional :: endpoint !! Include stop.
        complex(dp), intent(in), optional :: tracea !! Optional trace hint.
        real(dp), intent(in), optional :: tol !! Relative Taylor remainder tolerance.
        integer, intent(out), optional :: products !! Matrix-vector action count.
        complex(dp), allocatable :: block(:,:), result(:,:,:)

        allocate(block(size(b),1))
        block(:,1) = b
        call interval_csc_block(a,block,start,stop,num,result,info,endpoint,tracea,tol,products)
        if (size(result,2) == 1) then
            y = result(:,1,:)
        else
            allocate(y(0,0))
        end if
    end subroutine interval_csc_vector

    subroutine interval_operator_block(n, b, matvec, norm_bound, start, stop, num, y, info, &
                                        endpoint, tracea, tol, products)
        integer, intent(in) :: n !! Operator dimension.
        complex(dp), intent(in) :: b(:,:) !! Finite input block (n,nrhs).
        procedure(complex_expm_matvec) :: matvec !! Callback A*x.
        real(dp), intent(in) :: norm_bound !! Upper bound on infinity norm of A.
        real(dp), intent(in) :: start !! First time.
        real(dp), intent(in) :: stop !! Last time bound.
        integer, intent(in) :: num !! Number of samples.
        complex(dp), allocatable, intent(out) :: y(:,:,:) !! Output (n,nrhs,num).
        integer, intent(out) :: info !! Status code.
        logical, intent(in), optional :: endpoint !! Include stop.
        complex(dp), intent(in), optional :: tracea !! Optional trace of A.
        real(dp), intent(in), optional :: tol !! Relative Taylor remainder tolerance.
        integer, intent(out), optional :: products !! Callback count.

        allocate(y(0,0,0))
        info = complex_expm_invalid_input
        if (present(products)) products = 0
        if (n < 0 .or. size(b,1) /= n) return
        if (.not. ieee_is_finite(norm_bound) .or. norm_bound < 0.0_dp) return
        if (present(tracea)) then
            if (.not. finite_complex(tracea)) return
        end if
        deallocate(y)
        call interval_core(n,b,start,stop,num,step,y,info,endpoint,products)
    contains
        subroutine step(time,x,z,status,nprod)
            real(dp), intent(in) :: time !! Requested time.
            complex(dp), intent(in) :: x(:,:) !! Input block.
            complex(dp), allocatable, intent(out) :: z(:,:) !! Output block.
            integer, intent(out) :: status !! Computation status.
            integer, intent(out) :: nprod !! Callback count.
            call expm_operator_block(n,x,matvec,norm_bound,z,status,t=time, &
                                     tracea=tracea,tol=tol,products=nprod)
        end subroutine step
    end subroutine interval_operator_block

    subroutine interval_operator_vector(n, b, matvec, norm_bound, start, stop, num, y, info, &
                                         endpoint, tracea, tol, products)
        integer, intent(in) :: n !! Operator dimension.
        complex(dp), intent(in) :: b(:) !! Finite input vector.
        procedure(complex_expm_matvec) :: matvec !! Callback A*x.
        real(dp), intent(in) :: norm_bound !! Upper bound on infinity norm of A.
        real(dp), intent(in) :: start !! First time.
        real(dp), intent(in) :: stop !! Last time bound.
        integer, intent(in) :: num !! Number of samples.
        complex(dp), allocatable, intent(out) :: y(:,:) !! Output (n,num).
        integer, intent(out) :: info !! Status code.
        logical, intent(in), optional :: endpoint !! Include stop.
        complex(dp), intent(in), optional :: tracea !! Optional trace of A.
        real(dp), intent(in), optional :: tol !! Relative Taylor remainder tolerance.
        integer, intent(out), optional :: products !! Callback count.
        complex(dp), allocatable :: block(:,:), result(:,:,:)

        allocate(block(size(b),1))
        block(:,1) = b
        call interval_operator_block(n,block,matvec,norm_bound,start,stop,num,result,info, &
                                     endpoint,tracea,tol,products)
        if (size(result,2) == 1) then
            y = result(:,1,:)
        else
            allocate(y(0,0))
        end if
    end subroutine interval_operator_vector

end module scifort_sparse_linalg_complex_expm
