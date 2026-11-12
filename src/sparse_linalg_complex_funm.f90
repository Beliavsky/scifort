! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Actions of general analytic functions of complex sparse matrices.
! Reorthogonalized Arnoldi (or Hermitian Lanczos) with block-Hessenberg
! restarted matrix-function projections, following the public SciFort real
! funm_multiply_krylov interface. The dense f(H) evaluator is supplied by
! the caller. Update-based stopping is heuristic, NOT a rigorous error bound.
module scifort_sparse_linalg_complex_funm
    use, intrinsic :: iso_fortran_env, only : int64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_sparse_complex, only : csr_complex_matrix, csc_complex_matrix, &
        valid_csr_complex, valid_csc_complex, csc_complex_to_csr
    implicit none
    private

    integer, parameter, public :: complex_funm_success = 0
    integer, parameter, public :: complex_funm_no_convergence = 1
    integer, parameter, public :: complex_funm_invalid_input = -1
    integer, parameter, public :: complex_funm_breakdown = -2
    integer, parameter, public :: complex_funm_callback_failure = -3

    public :: funm_multiply_krylov_complex, funm_multiply_krylov_complex_operator
    public :: complex_funm_dense_function, complex_funm_matvec

    abstract interface
        subroutine complex_funm_dense_function(matrix, value, info)
            import dp
            complex(dp), intent(in) :: matrix(:,:) !! Small square projected matrix.
            complex(dp), intent(out) :: value(:,:) !! f(matrix), same square shape.
            integer, intent(out) :: info !! Zero on success; nonzero otherwise.
        end subroutine complex_funm_dense_function

        subroutine complex_funm_matvec(x, y, info)
            import dp
            complex(dp), intent(in) :: x(:) !! Input vector, length n.
            complex(dp), intent(out) :: y(:) !! A*x, length n.
            integer, intent(out) :: info !! Zero on success; nonzero otherwise.
        end subroutine complex_funm_matvec
    end interface

    interface funm_multiply_krylov_complex
        module procedure funm_csr_vector, funm_csc_vector
        module procedure funm_csr_block, funm_csc_block
    end interface funm_multiply_krylov_complex

    interface funm_multiply_krylov_complex_operator
        module procedure funm_operator_vector, funm_operator_block
    end interface funm_multiply_krylov_complex_operator

contains

    pure elemental logical function finite_complex(z) result(ok)
        complex(dp), intent(in) :: z !! Candidate complex value.
        ok = ieee_is_finite(real(z,dp)) .and. ieee_is_finite(aimag(z))
    end function finite_complex

    ! Scaled Euclidean norm, avoiding overflow in abs(z)**2.
    pure real(dp) function complex_norm(x) result(value)
        complex(dp), intent(in) :: x(:) !! Input complex vector.
        real(dp) :: scale, ssq, component
        integer :: k, part

        scale = 0.0_dp
        ssq = 1.0_dp
        do k = 1, size(x)
            do part = 1, 2
                if (part == 1) then
                    component = abs(real(x(k),dp))
                else
                    component = abs(aimag(x(k)))
                end if
                if (component == 0.0_dp) cycle
                if (scale < component) then
                    ssq = 1.0_dp + ssq*(scale/component)**2
                    scale = component
                else
                    ssq = ssq + (component/scale)**2
                end if
            end do
        end do
        value = scale*sqrt(ssq)
    end function complex_norm

    ! Validate options shared by vector and block entry points, including
    ! zero-column blocks, where no vector solver would otherwise be called.
    pure logical function valid_options(t,rtol,atol,restart_every_m,max_restarts) result(ok)
        complex(dp), intent(in), optional :: t !! Complex time scale.
        real(dp), intent(in), optional :: rtol !! Relative stopping tolerance.
        real(dp), intent(in), optional :: atol !! Absolute stopping tolerance.
        integer, intent(in), optional :: restart_every_m !! Basis size, positive.
        integer, intent(in), optional :: max_restarts !! Cycle budget, positive.
        ok = .false.
        if (present(t)) then
            if (.not. finite_complex(t)) return
        end if
        if (present(rtol)) then
            if (.not. ieee_is_finite(rtol)) return
            if (rtol < 0.0_dp) return
        end if
        if (present(atol)) then
            if (.not. ieee_is_finite(atol)) return
            if (atol < 0.0_dp) return
        end if
        if (present(restart_every_m)) then
            if (restart_every_m <= 0) return
        end if
        if (present(max_restarts)) then
            if (max_restarts <= 0) return
        end if
        ok = .true.
    end function valid_options

    subroutine funm_csr_vector(a, b, f, y, info, t, rtol, atol, &
                               restart_every_m, max_restarts, assume_hermitian, &
                               cycles, products, update_norm)
        type(csr_complex_matrix), intent(in) :: a !! Canonical complex CSR matrix.
        complex(dp), intent(in) :: b(:) !! Right hand vector, length a%n_rows.
        procedure(complex_funm_dense_function) :: f !! Projected matrix function.
        complex(dp), allocatable, intent(out) :: y(:) !! Approximation to f(t*A)*b.
        integer, intent(out) :: info !! Zero success, positive unconverged, negative error.
        complex(dp), intent(in), optional :: t !! Complex scalar multiplier, default 1.
        real(dp), intent(in), optional :: rtol !! Relative update tolerance, default 1e-6.
        real(dp), intent(in), optional :: atol !! Absolute update tolerance, default 0.
        integer, intent(in), optional :: restart_every_m !! Krylov vectors/cycle; default min(n,20).
        integer, intent(in), optional :: max_restarts !! Maximum cycles, default 20.
        logical, intent(in), optional :: assume_hermitian !! A is Hermitian; use Lanczos.
        integer, intent(out), optional :: cycles !! Number of completed Arnoldi cycles.
        integer, intent(out), optional :: products !! Number of A*x operations.
        real(dp), intent(out), optional :: update_norm !! Norm of last correction.
        if (.not. valid_csr_complex(a)) then
            allocate(y(0))
            info = complex_funm_invalid_input
            if (present(cycles)) cycles = 0
            if (present(products)) products = 0
            if (present(update_norm)) update_norm = 0.0_dp
            return
        end if
        if (a%n_rows /= a%n_cols .or. size(b) /= a%n_rows) then
            allocate(y(0))
            info = complex_funm_invalid_input
            if (present(cycles)) cycles = 0
            if (present(products)) products = 0
            if (present(update_norm)) update_norm = 0.0_dp
            return
        end if
        call funm_operator_vector(a%n_rows,b,csr_action,f,y,info,t,rtol,atol, &
            restart_every_m,max_restarts,assume_hermitian,cycles,products,update_norm)
    contains
        subroutine csr_action(x, ax, status)
            complex(dp), intent(in) :: x(:) !! Vector to multiply.
            complex(dp), intent(out) :: ax(:) !! A*x.
            integer, intent(out) :: status !! Zero for success.
            integer :: row, k
            status = 0
            ax = cmplx(0.0_dp,0.0_dp,dp)
            do row = 1, a%n_rows
                do k = a%indptr(row), a%indptr(row+1)-1
                    ax(row) = ax(row) + a%data(k)*x(a%indices(k))
                end do
            end do
        end subroutine csr_action
    end subroutine funm_csr_vector

    subroutine funm_csc_vector(a, b, f, y, info, t, rtol, atol, &
                               restart_every_m, max_restarts, assume_hermitian, &
                               cycles, products, update_norm)
        type(csc_complex_matrix), intent(in) :: a !! Canonical complex CSC matrix.
        complex(dp), intent(in) :: b(:) !! Right hand vector.
        procedure(complex_funm_dense_function) :: f !! Matrix-function callback.
        complex(dp), allocatable, intent(out) :: y(:) !! f(t*A)*b approximation.
        integer, intent(out) :: info !! Solver status.
        complex(dp), intent(in), optional :: t !! Complex multiplier, default 1.
        real(dp), intent(in), optional :: rtol !! Relative update tolerance.
        real(dp), intent(in), optional :: atol !! Absolute update tolerance.
        integer, intent(in), optional :: restart_every_m !! Krylov vectors per cycle.
        integer, intent(in), optional :: max_restarts !! Maximum cycles.
        logical, intent(in), optional :: assume_hermitian !! Hermitian A assumed.
        integer, intent(out), optional :: cycles !! Number of completed cycles.
        integer, intent(out), optional :: products !! A*x count.
        real(dp), intent(out), optional :: update_norm !! Last update magnitude.
        type(csr_complex_matrix) :: csr
        integer :: status

        if (.not. valid_csc_complex(a)) then
            allocate(y(0))
            info = complex_funm_invalid_input
            if (present(cycles)) cycles = 0
            if (present(products)) products = 0
            if (present(update_norm)) update_norm = 0.0_dp
            return
        end if
        csr = csc_complex_to_csr(a,status)
        if (status /= 0) then
            allocate(y(0))
            info = complex_funm_invalid_input
            if (present(cycles)) cycles = 0
            if (present(products)) products = 0
            if (present(update_norm)) update_norm = 0.0_dp
            return
        end if
        call funm_csr_vector(csr,b,f,y,info,t,rtol,atol, &
            restart_every_m,max_restarts,assume_hermitian,cycles,products,update_norm)
    end subroutine funm_csc_vector

    subroutine funm_csr_block(a, b, f, y, info, t, rtol, atol, &
                              restart_every_m, max_restarts, assume_hermitian, &
                              cycles, products, update_norm)
        type(csr_complex_matrix), intent(in) :: a !! Square CSR input.
        complex(dp), intent(in) :: b(:,:) !! Dense RHS block, (n,nrhs).
        procedure(complex_funm_dense_function) :: f !! Dense function evaluator.
        complex(dp), allocatable, intent(out) :: y(:,:) !! f(t*A)*b.
        integer, intent(out) :: info !! Status of whole block solve.
        complex(dp), intent(in), optional :: t !! Complex multiplier.
        real(dp), intent(in), optional :: rtol !! Relative tolerance.
        real(dp), intent(in), optional :: atol !! Absolute tolerance.
        integer, intent(in), optional :: restart_every_m !! Local Krylov size.
        integer, intent(in), optional :: max_restarts !! Maximum cycles per RHS.
        logical, intent(in), optional :: assume_hermitian !! Hermitian A assumed.
        integer, intent(out), optional :: cycles !! Sum of completed cycles across RHS.
        integer, intent(out), optional :: products !! Total A*x calls across RHS.
        real(dp), intent(out), optional :: update_norm !! Largest last correction among RHS.
        complex(dp), allocatable :: column(:)
        integer :: j, done, applied, cumulative_cycles, cumulative_products, status
        real(dp) :: last_update, largest_update

        allocate(y(0,0))
        info = complex_funm_invalid_input
        if (present(cycles)) cycles = 0
        if (present(products)) products = 0
        if (present(update_norm)) update_norm = 0.0_dp
        if (.not. valid_csr_complex(a)) return
        if (a%n_rows /= a%n_cols .or. size(b,1) /= a%n_rows) return
        if (.not. all(finite_complex(b))) return
        if (.not. valid_options(t,rtol,atol,restart_every_m,max_restarts)) return
        deallocate(y)
        allocate(y(size(b,1),size(b,2)))
        y = cmplx(0.0_dp,0.0_dp,dp)
        cumulative_cycles = 0
        cumulative_products = 0
        largest_update = 0.0_dp
        info = complex_funm_success
        do j = 1, size(b,2)
            call funm_csr_vector(a,b(:,j),f,column,status,t,rtol,atol, &
                restart_every_m,max_restarts,assume_hermitian,done,applied,last_update)
            cumulative_cycles = cumulative_cycles + done
            cumulative_products = cumulative_products + applied
            largest_update = max(largest_update,last_update)
            if (size(column) == size(b,1)) y(:,j) = column
            if (status /= complex_funm_success) then
                info = status
                exit
            end if
        end do
        if (present(cycles)) cycles = cumulative_cycles
        if (present(products)) products = cumulative_products
        if (present(update_norm)) update_norm = largest_update
    end subroutine funm_csr_block

    subroutine funm_csc_block(a, b, f, y, info, t, rtol, atol, &
                              restart_every_m, max_restarts, assume_hermitian, &
                              cycles, products, update_norm)
        type(csc_complex_matrix), intent(in) :: a !! Square CSC input.
        complex(dp), intent(in) :: b(:,:) !! RHS block.
        procedure(complex_funm_dense_function) :: f !! Dense function evaluator.
        complex(dp), allocatable, intent(out) :: y(:,:) !! Result block.
        integer, intent(out) :: info !! Solver status.
        complex(dp), intent(in), optional :: t !! Complex multiplier.
        real(dp), intent(in), optional :: rtol !! Relative tolerance.
        real(dp), intent(in), optional :: atol !! Absolute tolerance.
        integer, intent(in), optional :: restart_every_m !! Local Krylov size.
        integer, intent(in), optional :: max_restarts !! Maximum cycles per RHS.
        logical, intent(in), optional :: assume_hermitian !! Hermitian A assumed.
        integer, intent(out), optional :: cycles !! Sum of cycles.
        integer, intent(out), optional :: products !! Total A*x calls.
        real(dp), intent(out), optional :: update_norm !! Largest last update.
        type(csr_complex_matrix) :: csr
        integer :: status

        allocate(y(0,0))
        info = complex_funm_invalid_input
        if (present(cycles)) cycles = 0
        if (present(products)) products = 0
        if (present(update_norm)) update_norm = 0.0_dp
        if (.not. valid_csc_complex(a)) return
        csr = csc_complex_to_csr(a,status)
        if (status /= 0) return
        call funm_csr_block(csr,b,f,y,info,t,rtol,atol, &
            restart_every_m,max_restarts,assume_hermitian,cycles,products,update_norm)
    end subroutine funm_csc_block

    subroutine funm_operator_block(n, b, matvec, f, y, info, t, rtol, atol, &
                                   restart_every_m, max_restarts, assume_hermitian, &
                                   cycles, products, update_norm)
        integer, intent(in) :: n !! Matrix-free operator order.
        complex(dp), intent(in) :: b(:,:) !! Input RHS block, shape (n,nrhs).
        procedure(complex_funm_matvec) :: matvec !! Callback for A*x.
        procedure(complex_funm_dense_function) :: f !! Dense projected matrix function.
        complex(dp), allocatable, intent(out) :: y(:,:) !! Computed f(t*A)*b.
        integer, intent(out) :: info !! Solver status.
        complex(dp), intent(in), optional :: t !! Complex argument scale.
        real(dp), intent(in), optional :: rtol !! Relative tolerance.
        real(dp), intent(in), optional :: atol !! Absolute tolerance.
        integer, intent(in), optional :: restart_every_m !! Local Krylov size.
        integer, intent(in), optional :: max_restarts !! Max cycles per RHS.
        logical, intent(in), optional :: assume_hermitian !! Hermitian operator assumed.
        integer, intent(out), optional :: cycles !! Sum of cycles across RHS.
        integer, intent(out), optional :: products !! Total matvec calls.
        real(dp), intent(out), optional :: update_norm !! Largest last update.
        complex(dp), allocatable :: column(:)
        integer :: j, done, applied, status, total_cycles, total_products
        real(dp) :: last_update, largest_update

        allocate(y(0,0))
        info = complex_funm_invalid_input
        if (present(cycles)) cycles = 0
        if (present(products)) products = 0
        if (present(update_norm)) update_norm = 0.0_dp
        if (n < 0 .or. size(b,1) /= n) return
        if (.not. all(finite_complex(b))) return
        if (.not. valid_options(t,rtol,atol,restart_every_m,max_restarts)) return
        deallocate(y)
        allocate(y(n,size(b,2)))
        y = cmplx(0.0_dp,0.0_dp,dp)
        total_cycles = 0
        total_products = 0
        largest_update = 0.0_dp
        info = complex_funm_success
        do j = 1, size(b,2)
            call funm_operator_vector(n,b(:,j),matvec,f,column,status,t,rtol,atol, &
                restart_every_m,max_restarts,assume_hermitian,done,applied,last_update)
            total_cycles = total_cycles + done
            total_products = total_products + applied
            largest_update = max(largest_update,last_update)
            if (size(column) == n) y(:,j) = column
            if (status /= complex_funm_success) then
                info = status
                exit
            end if
        end do
        if (present(cycles)) cycles = total_cycles
        if (present(products)) products = total_products
        if (present(update_norm)) update_norm = largest_update
    end subroutine funm_operator_block

    subroutine funm_operator_vector(n, b, matvec, f, y, info, t, rtol, atol, &
                                    restart_every_m, max_restarts, assume_hermitian, &
                                    cycles, products, update_norm)
        integer, intent(in) :: n !! Nonnegative matrix order.
        complex(dp), intent(in) :: b(:) !! Finite vector with n entries.
        procedure(complex_funm_matvec) :: matvec !! A*x; status=0 on success.
        procedure(complex_funm_dense_function) :: f !! f(H) evaluator for small H.
        complex(dp), allocatable, intent(out) :: y(:) !! f(t*A)*b approximation.
        integer, intent(out) :: info !! 0 success, 1 exhausted, <0 invalid or callback error.
        complex(dp), intent(in), optional :: t !! Complex scale, default (1,0).
        real(dp), intent(in), optional :: rtol !! Relative update tolerance, default 1e-6.
        real(dp), intent(in), optional :: atol !! Absolute update tolerance, default 0.
        integer, intent(in), optional :: restart_every_m !! Maximum basis columns/cycle.
        integer, intent(in), optional :: max_restarts !! Maximum cycles, default 20.
        logical, intent(in), optional :: assume_hermitian !! A Hermitian (not t*A); default false.
        integer, intent(out), optional :: cycles !! Finished Arnoldi cycles.
        integer, intent(out), optional :: products !! Matvec calls, including failed calls.
        real(dp), intent(out), optional :: update_norm !! Norm of last correction.
        complex(dp) :: time_scale
        complex(dp), allocatable :: v(:,:), h(:,:), fh(:,:), seed(:), candidate(:), grown_h(:,:)
        real(dp) :: reltol, abstol, bnorm, threshold, change
        integer :: m, cap, capacity, next_capacity, ic, base, j, global_j, &
                   nruns, nproducts, status
        logical :: hermitian, happy

        allocate(y(0))
        info = complex_funm_invalid_input
        nruns = 0
        nproducts = 0
        change = 0.0_dp
        if (present(cycles)) cycles = 0
        if (present(products)) products = 0
        if (present(update_norm)) update_norm = 0.0_dp
        if (n < 0 .or. size(b) /= n) return
        if (.not. all(finite_complex(b))) return
        time_scale = cmplx(1.0_dp,0.0_dp,dp)
        reltol = 1.0e-6_dp
        abstol = 0.0_dp
        cap = 20
        hermitian = .false.
        if (present(t)) time_scale = t
        if (present(rtol)) reltol = rtol
        if (present(atol)) abstol = atol
        if (present(max_restarts)) cap = max_restarts
        if (present(assume_hermitian)) hermitian = assume_hermitian
        if (.not. finite_complex(time_scale)) return
        if (.not. ieee_is_finite(reltol) .or. .not. ieee_is_finite(abstol)) return
        if (reltol < 0.0_dp .or. abstol < 0.0_dp .or. cap <= 0) return
        if (present(restart_every_m)) then
            if (restart_every_m <= 0) return
        end if
        if (n == 0) then
            info = complex_funm_success
            return
        end if
        m = min(20,n)
        if (present(restart_every_m)) m = min(n,restart_every_m)
        if (m >= n) cap = 1
        bnorm = complex_norm(b)
        if (.not. ieee_is_finite(bnorm)) then
            info = complex_funm_breakdown
            return
        end if
        deallocate(y)
        allocate(y(n))
        y = cmplx(0.0_dp,0.0_dp,dp)
        if (bnorm == 0.0_dp) then
            info = complex_funm_success
            return
        end if
        threshold = max(abstol,reltol*bnorm)
        if (.not. ieee_is_finite(threshold)) then
            info = complex_funm_invalid_input
            return
        end if
        capacity = min(cap,2)
        if (m > (huge(m)-1)/capacity) return
        allocate(v(n,m+1), h(m*capacity+1,m*capacity), &
                 fh(m*capacity,m*capacity), seed(n), candidate(n))
        h = cmplx(0.0_dp,0.0_dp,dp)
        seed = b/bnorm
        info = complex_funm_no_convergence
        do ic = 1, cap
            if (ic > capacity) then
                next_capacity = int(min(int(cap,int64),2_int64*int(capacity,int64)))
                if (m > (huge(m)-1)/next_capacity) then
                    info = complex_funm_invalid_input
                    exit
                end if
                allocate(grown_h(m*next_capacity+1,m*next_capacity))
                grown_h = cmplx(0.0_dp,0.0_dp,dp)
                grown_h(1:size(h,1),1:size(h,2)) = h
                call move_alloc(grown_h,h)
                deallocate(fh)
                allocate(fh(m*next_capacity,m*next_capacity))
                capacity = next_capacity
            end if
            base = (ic-1)*m
            call funm_krylov_cycle(matvec,seed,v,h,base,m,hermitian, &
                                  j,happy,nproducts,status)
            if (status /= complex_funm_success) then
                info = status
                exit
            end if
            global_j = base+j
            if (.not. all(finite_complex(time_scale*h(1:global_j,1:global_j)))) then
                info = complex_funm_breakdown
                exit
            end if
            call f(time_scale*h(1:global_j,1:global_j), &
                   fh(1:global_j,1:global_j),status)
            if (status /= 0) then
                info = complex_funm_callback_failure
                exit
            end if
            if (.not. all(finite_complex(fh(1:global_j,1:global_j)))) then
                info = complex_funm_breakdown
                exit
            end if
            candidate = bnorm*matmul(v(:,1:j),fh(base+1:global_j,1))
            if (.not. all(finite_complex(candidate))) then
                info = complex_funm_breakdown
                exit
            end if
            if (.not. all(finite_complex(y+candidate))) then
                info = complex_funm_breakdown
                exit
            end if
            y = y+candidate
            change = complex_norm(candidate)
            nruns = ic
            if (happy .or. (ic > 1 .and. change <= threshold)) then
                info = complex_funm_success
                exit
            end if
            if (ic == cap) exit
            seed = v(:,m+1)
        end do
        if (present(cycles)) cycles = nruns
        if (present(products)) products = nproducts
        if (present(update_norm)) update_norm = change
    end subroutine funm_operator_vector

    subroutine funm_krylov_cycle(matvec,seed,v,h,base,m,hermitian,used,happy,products,info)
        procedure(complex_funm_matvec) :: matvec !! Matrix-vector callback.
        complex(dp), intent(in) :: seed(:) !! Nonzero normalized start direction.
        complex(dp), intent(out) :: v(:,:) !! Local (n,m+1) basis.
        complex(dp), intent(inout) :: h(:,:) !! Extended Hessenberg matrix.
        integer, intent(in) :: base !! Columns in previously built blocks.
        integer, intent(in) :: m !! Maximum Krylov dimension this cycle.
        logical, intent(in) :: hermitian !! Use Lanczos three-term recurrence.
        integer, intent(out) :: used !! Number of basis vectors consumed.
        logical, intent(out) :: happy !! Invariant subspace detected.
        integer, intent(inout) :: products !! Total operator applications.
        integer, intent(out) :: info !! Status code.
        complex(dp) :: w(size(seed)), projection
        real(dp) :: beta, alpha, before
        integer :: k, j, pass, status

        v = cmplx(0.0_dp,0.0_dp,dp)
        v(:,1) = seed/complex_norm(seed)
        used = 0
        happy = .false.
        info = complex_funm_success
        do k = 1, m
            call matvec(v(:,k),w,status)
            products = products+1
            if (status /= 0) then
                info = complex_funm_callback_failure
                return
            end if
            if (.not. all(finite_complex(w))) then
                info = complex_funm_breakdown
                return
            end if
            before = complex_norm(w)
            if (hermitian) then
                if (k > 1) w = w-h(base+k,base+k-1)*v(:,k-1)
                alpha = real(dot_product(v(:,k),w),dp)
                h(base+k,base+k) = cmplx(alpha,0.0_dp,dp)
                w = w-alpha*v(:,k)
            else
                do pass = 1,2
                    do j = 1,k
                        projection = dot_product(v(:,j),w)
                        h(base+j,base+k) = h(base+j,base+k)+projection
                        w = w-projection*v(:,j)
                    end do
                end do
            end if
            beta = complex_norm(w)
            if (.not. ieee_is_finite(beta)) then
                info = complex_funm_breakdown
                return
            end if
            used = k
            h(base+k+1,base+k) = cmplx(beta,0.0_dp,dp)
            happy = beta <= 32.0_dp*epsilon(1.0_dp)*max(before,tiny(1.0_dp))
            if (k == size(seed)) happy = .true.
            if (happy) return
            v(:,k+1) = w/beta
            if (hermitian .and. k < m) h(base+k,base+k+1) = cmplx(beta,0.0_dp,dp)
        end do
    end subroutine funm_krylov_cycle

end module scifort_sparse_linalg_complex_funm
