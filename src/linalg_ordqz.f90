! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Ordered complex generalized Schur form. Adjacent triangular pencil blocks
! are exchanged by separate two-sided unitary transformations. No explicit
! inverse and no LAPACK dependency are needed for the reordering step.
module scifort_linalg_ordqz
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    use scifort_linalg_matrix_functions, only : matrix_function_max_order
    use scifort_linalg_qz, only : qz
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    implicit none
    private

    public :: ordqz, ordqz_factors, ordqz_selector

    abstract interface
        logical function ordqz_selector(alpha, beta)
            import dp
            complex(dp), intent(in) :: alpha !! homogeneous eigenvalue numerator
            complex(dp), intent(in) :: beta !! homogeneous eigenvalue denominator
        end function ordqz_selector
    end interface

    interface ordqz
        module procedure ordqz_rr
        module procedure ordqz_rc
        module procedure ordqz_cr
        module procedure ordqz_cc
    end interface ordqz

contains

    subroutine ordqz_rr(a, b, aa, bb, q, z, status, alpha, beta, sort, sdim, selector, method)
        real(dp), intent(in) :: a(:, :) !! first square matrix, order 0..40
        real(dp), intent(in) :: b(:, :) !! second square matrix, same order
        complex(dp), intent(out) :: aa(:, :) !! reordered upper-triangular factor of A
        complex(dp), intent(out) :: bb(:, :) !! reordered upper-triangular factor of B
        complex(dp), intent(out) :: q(:, :) !! left unitary Schur vectors
        complex(dp), intent(out) :: z(:, :) !! right unitary Schur vectors
        integer, intent(out) :: status !! linalg_status_* completion code
        complex(dp), intent(out), optional :: alpha(:) !! homogeneous numerators, size n
        complex(dp), intent(out), optional :: beta(:) !! homogeneous denominators, size n
        character(len=*), intent(in), optional :: sort !! lhp(default), rhp, iuc, ouc, none
        integer, intent(out), optional :: sdim !! number of selected eigenvalues
        procedure(ordqz_selector), optional :: selector !! custom selector of (alpha,beta)
        character(len=*), intent(in), optional :: method !! QZ method: auto, iterative, inverse

        call ordqz_core(cmplx(a,kind=dp),cmplx(b,kind=dp),aa,bb,q,z,status, &
            alpha,beta,sort,sdim,selector,method)
    end subroutine ordqz_rr

    subroutine ordqz_rc(a, b, aa, bb, q, z, status, alpha, beta, sort, sdim, selector, method)
        real(dp), intent(in) :: a(:, :) !! first real square matrix, order 0..40
        complex(dp), intent(in) :: b(:, :) !! second complex square matrix, same order
        complex(dp), intent(out) :: aa(:, :) !! reordered upper-triangular factor of A
        complex(dp), intent(out) :: bb(:, :) !! reordered upper-triangular factor of B
        complex(dp), intent(out) :: q(:, :) !! left unitary Schur vectors
        complex(dp), intent(out) :: z(:, :) !! right unitary Schur vectors
        integer, intent(out) :: status !! linalg_status_* completion code
        complex(dp), intent(out), optional :: alpha(:) !! homogeneous numerators, size n
        complex(dp), intent(out), optional :: beta(:) !! homogeneous denominators, size n
        character(len=*), intent(in), optional :: sort !! lhp(default), rhp, iuc, ouc, none
        integer, intent(out), optional :: sdim !! number of selected eigenvalues
        procedure(ordqz_selector), optional :: selector !! custom selector of (alpha,beta)
        character(len=*), intent(in), optional :: method !! QZ method: auto, iterative, inverse

        call ordqz_core(cmplx(a,kind=dp),b,aa,bb,q,z,status, &
            alpha,beta,sort,sdim,selector,method)
    end subroutine ordqz_rc

    subroutine ordqz_cr(a, b, aa, bb, q, z, status, alpha, beta, sort, sdim, selector, method)
        complex(dp), intent(in) :: a(:, :) !! first complex square matrix, order 0..40
        real(dp), intent(in) :: b(:, :) !! second real square matrix, same order
        complex(dp), intent(out) :: aa(:, :) !! reordered upper-triangular factor of A
        complex(dp), intent(out) :: bb(:, :) !! reordered upper-triangular factor of B
        complex(dp), intent(out) :: q(:, :) !! left unitary Schur vectors
        complex(dp), intent(out) :: z(:, :) !! right unitary Schur vectors
        integer, intent(out) :: status !! linalg_status_* completion code
        complex(dp), intent(out), optional :: alpha(:) !! homogeneous numerators, size n
        complex(dp), intent(out), optional :: beta(:) !! homogeneous denominators, size n
        character(len=*), intent(in), optional :: sort !! lhp(default), rhp, iuc, ouc, none
        integer, intent(out), optional :: sdim !! number of selected eigenvalues
        procedure(ordqz_selector), optional :: selector !! custom selector of (alpha,beta)
        character(len=*), intent(in), optional :: method !! QZ method: auto, iterative, inverse

        call ordqz_core(a,cmplx(b,kind=dp),aa,bb,q,z,status, &
            alpha,beta,sort,sdim,selector,method)
    end subroutine ordqz_cr

    subroutine ordqz_cc(a, b, aa, bb, q, z, status, alpha, beta, sort, sdim, selector, method)
        complex(dp), intent(in) :: a(:, :) !! first complex square matrix, order 0..40
        complex(dp), intent(in) :: b(:, :) !! second complex square matrix, same order
        complex(dp), intent(out) :: aa(:, :) !! reordered upper-triangular factor of A
        complex(dp), intent(out) :: bb(:, :) !! reordered upper-triangular factor of B
        complex(dp), intent(out) :: q(:, :) !! left unitary Schur vectors
        complex(dp), intent(out) :: z(:, :) !! right unitary Schur vectors
        integer, intent(out) :: status !! linalg_status_* completion code
        complex(dp), intent(out), optional :: alpha(:) !! homogeneous numerators, size n
        complex(dp), intent(out), optional :: beta(:) !! homogeneous denominators, size n
        character(len=*), intent(in), optional :: sort !! lhp(default), rhp, iuc, ouc, none
        integer, intent(out), optional :: sdim !! number of selected eigenvalues
        procedure(ordqz_selector), optional :: selector !! custom selector of (alpha,beta)
        character(len=*), intent(in), optional :: method !! QZ method: auto, iterative, inverse

        call ordqz_core(a,b,aa,bb,q,z,status,alpha,beta,sort,sdim,selector,method)
    end subroutine ordqz_cc

    subroutine ordqz_core(a, b, aa, bb, q, z, status, alpha, beta, sort, sdim, selector, method)
        complex(dp), intent(in) :: a(:, :) !! square coefficient A
        complex(dp), intent(in) :: b(:, :) !! square coefficient B
        complex(dp), intent(out) :: aa(:, :) !! ordered triangular A factor
        complex(dp), intent(out) :: bb(:, :) !! ordered triangular B factor
        complex(dp), intent(out) :: q(:, :) !! left unitary matrix
        complex(dp), intent(out) :: z(:, :) !! right unitary matrix
        integer, intent(out) :: status !! status of factorization and reordering
        complex(dp), intent(out), optional :: alpha(:) !! ordered numerator diagonal
        complex(dp), intent(out), optional :: beta(:) !! ordered denominator diagonal
        character(len=*), intent(in), optional :: sort !! selection region, default lhp
        integer, intent(out), optional :: sdim !! number of selected eigenvalues
        procedure(ordqz_selector), optional :: selector !! custom homogeneous selector
        character(len=*), intent(in), optional :: method !! auto, iterative, inverse
        integer :: n, j

        aa = cmplx(0.0_dp,0.0_dp,dp)
        bb = cmplx(0.0_dp,0.0_dp,dp)
        q = cmplx(0.0_dp,0.0_dp,dp)
        z = cmplx(0.0_dp,0.0_dp,dp)
        if (present(alpha)) alpha = cmplx(0.0_dp,0.0_dp,dp)
        if (present(beta)) beta = cmplx(0.0_dp,0.0_dp,dp)
        if (present(sdim)) sdim = 0
        status = linalg_status_invalid_shape
        n = size(a,1)
        if (n > matrix_function_max_order) return
        if (size(a,2) /= n .or. size(b,1) /= n .or. size(b,2) /= n) return
        if (size(aa,1) /= n .or. size(aa,2) /= n) return
        if (size(bb,1) /= n .or. size(bb,2) /= n) return
        if (size(q,1) /= n .or. size(q,2) /= n) return
        if (size(z,1) /= n .or. size(z,2) /= n) return
        if (present(alpha)) then
            if (size(alpha) /= n) return
        end if
        if (present(beta)) then
            if (size(beta) /= n) return
        end if
        if (.not. valid_sort(sort)) return

        call qz(a,b,aa,bb,q,z,status,method=method)
        if (status /= linalg_status_success) return
        call ordqz_factors(aa,bb,q,z,status,sort,sdim,selector)
        if (status /= linalg_status_success) return
        do j=1,n
            if (present(alpha)) alpha(j) = aa(j,j)
            if (present(beta)) beta(j) = bb(j,j)
        end do
    end subroutine ordqz_core

    pure logical function valid_sort(sort) result(ok)
        character(len=*), intent(in), optional :: sort !! sort code, absent means lhp
        character(len=:), allocatable :: choice
        choice = 'lhp'
        if (present(sort)) choice = trim(adjustl(sort))
        ok = choice == 'lhp' .or. choice == 'rhp' .or. choice == 'iuc' .or. &
            choice == 'ouc' .or. choice == 'none'
    end function valid_sort

    pure logical function finite_matrix(a) result(ok)
        complex(dp), intent(in) :: a(:, :) !! matrix whose components are inspected
        ok = all(ieee_is_finite(real(a,dp))) .and. all(ieee_is_finite(aimag(a)))
    end function finite_matrix

    ! Avoid alpha/beta: beta=0 represents an infinite eigenvalue. The
    ! half-plane test uses the sign of real(alpha*conjg(beta)), normalized to
    ! avoid overflow/underflow when coefficients differ by many orders.
    logical function chosen(alpha, beta, sort, beta_floor, selector) result(yes)
        complex(dp), intent(in) :: alpha !! homogeneous eigenvalue numerator
        complex(dp), intent(in) :: beta !! homogeneous eigenvalue denominator
        character(len=*), intent(in) :: sort !! one of the documented region codes
        real(dp), intent(in) :: beta_floor !! threshold for numerical infinite eigenvalues
        procedure(ordqz_selector), optional :: selector !! optional user predicate
        real(dp) :: ma, mb, projection
        if (present(selector)) then
            yes = selector(alpha,beta)
            return
        end if
        yes = .false.
        ma = abs(alpha)
        mb = abs(beta)
        if (ma <= 0.0_dp .and. mb <= beta_floor) return
        select case(sort)
        case('lhp','rhp')
            if (ma <= 0.0_dp .or. mb <= beta_floor) return
            projection = real((alpha/ma)*conjg(beta/mb),dp)
            if (sort == 'lhp') yes = projection < 0.0_dp
            if (sort == 'rhp') yes = projection > 0.0_dp
        case('iuc')
            yes = mb > beta_floor .and. ma < mb
        case('ouc')
            yes = (mb <= beta_floor .and. ma > 0.0_dp) .or. ma > mb
        case('none')
            yes = .false.
        end select
    end function chosen

    ! Swap two consecutive 1x1 complex Schur blocks in a regular triangular
    ! pencil. The first column of the right unitary rotation is the local
    ! generalized eigenvector of the second block. Left unitary alignment
    ! then makes both transformed (2,1) entries zero.
    subroutine adjacent_swap(a, b, q, z, k, status)
        complex(dp), intent(inout) :: a(:, :) !! normalized first triangular Schur factor
        complex(dp), intent(inout) :: b(:, :) !! normalized second triangular Schur factor
        complex(dp), intent(inout) :: q(:, :) !! accumulated left unitary matrix
        complex(dp), intent(inout) :: z(:, :) !! accumulated right unitary matrix
        integer, intent(in) :: k !! first block index; exchange k and k+1
        integer, intent(out) :: status !! success or numerical failure
        complex(dp) :: x(2), y(2), ya(2), yb(2), r(2,2), l(2,2)
        complex(dp) :: a11, a12, a22, b11, b12, b22
        real(dp) :: sa, sb, nx, ny, na, nb, tolerance

        status = linalg_status_no_convergence
        sa = maxval(abs(a(k:k+1,k:k+1)))
        sb = maxval(abs(b(k:k+1,k:k+1)))
        if (sa <= 0.0_dp) sa = 1.0_dp
        if (sb <= 0.0_dp) sb = 1.0_dp
        a11 = a(k,k)/sa
        a12 = a(k,k+1)/sa
        a22 = a(k+1,k+1)/sa
        b11 = b(k,k)/sb
        b12 = b(k,k+1)/sb
        b22 = b(k+1,k+1)/sb
        x(1) = a22*b12 - b22*a12
        x(2) = b22*a11 - a22*b11
        nx = hypot(abs(x(1)),abs(x(2)))
        if (.not. ieee_is_finite(nx) .or. nx <= 0.0_dp) return
        x = x / nx
        r(1,1) = x(1)
        r(2,1) = x(2)
        r(1,2) = -conjg(x(2))
        r(2,2) = conjg(x(1))
        ya = matmul(a(k:k+1,k:k+1)/sa,x)
        yb = matmul(b(k:k+1,k:k+1)/sb,x)
        na = hypot(abs(ya(1)),abs(ya(2)))
        nb = hypot(abs(yb(1)),abs(yb(2)))
        if (na >= nb) then
            y = ya
            ny = na
        else
            y = yb
            ny = nb
        end if
        if (.not. ieee_is_finite(ny) .or. ny <= 0.0_dp) return
        y = y / ny
        l(1,1) = y(1)
        l(2,1) = y(2)
        l(1,2) = -conjg(y(2))
        l(2,2) = conjg(y(1))

        a(:,k:k+1) = matmul(a(:,k:k+1),r)
        b(:,k:k+1) = matmul(b(:,k:k+1),r)
        a(k:k+1,:) = matmul(transpose(conjg(l)),a(k:k+1,:))
        b(k:k+1,:) = matmul(transpose(conjg(l)),b(k:k+1,:))
        q(:,k:k+1) = matmul(q(:,k:k+1),l)
        z(:,k:k+1) = matmul(z(:,k:k+1),r)
        tolerance = 1024.0_dp*epsilon(1.0_dp)
        if (max(abs(a(k+1,k)),abs(b(k+1,k))) > tolerance) return
        a(k+1,k) = cmplx(0.0_dp,0.0_dp,dp)
        b(k+1,k) = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_success
    end subroutine adjacent_swap

    ! Rearrange already computed complex Schur factors, including factors of
    ! order >40 supplied by another implementation. On failure, inputs are
    ! unchanged. Reordering a full eigensystem never forms B**(-1)A.
    subroutine ordqz_factors(aa, bb, q, z, status, sort, sdim, selector)
        complex(dp), intent(inout) :: aa(:, :) !! first complex triangular Schur factor
        complex(dp), intent(inout) :: bb(:, :) !! second complex triangular Schur factor
        complex(dp), intent(inout) :: q(:, :) !! left unitary Schur matrix
        complex(dp), intent(inout) :: z(:, :) !! right unitary Schur matrix
        integer, intent(out) :: status !! success or error; original data unchanged on error
        character(len=*), intent(in), optional :: sort !! lhp(default), rhp, iuc, ouc, none
        integer, intent(out), optional :: sdim !! number of selected eigenvalues
        procedure(ordqz_selector), optional :: selector !! homogeneous eigenvalue predicate
        complex(dp), allocatable :: wa(:,:), wb(:,:), wq(:,:), wz(:,:), original_a(:,:), original_b(:,:)
        complex(dp), allocatable :: final_a(:,:), final_b(:,:)
        real(dp) :: sa, sb, error_a, error_b, tolerance, beta_floor
        logical, allocatable :: selected(:)
        logical :: swap_flag
        integer :: n, i, j, k, count, local_status
        character(len=:), allocatable :: choice

        if (present(sdim)) sdim = 0
        status = linalg_status_invalid_shape
        n = size(aa,1)
        if (size(aa,2) /= n .or. size(bb,1) /= n .or. size(bb,2) /= n) return
        if (size(q,1) /= n .or. size(q,2) /= n) return
        if (size(z,1) /= n .or. size(z,2) /= n) return
        if (.not. valid_sort(sort)) return
        choice = 'lhp'
        if (present(sort)) choice = trim(adjustl(sort))
        status = linalg_status_no_convergence
        if (.not. finite_matrix(aa) .or. .not. finite_matrix(bb)) return
        if (.not. finite_matrix(q) .or. .not. finite_matrix(z)) return
        if (n == 0) then
            status = linalg_status_success
            return
        end if
        sa = maxval(abs(aa))
        sb = maxval(abs(bb))
        if (sa <= 0.0_dp) sa = 1.0_dp
        if (sb <= 0.0_dp) sb = 1.0_dp
        allocate(wa(n,n),wb(n,n),wq(n,n),wz(n,n),selected(n))
        allocate(original_a(n,n),original_b(n,n))
        wa = aa/sa
        wb = bb/sb
        tolerance = 1024.0_dp*epsilon(1.0_dp)
        do j=1,n-1
            if (maxval(abs(wa(j+1:n,j))) > tolerance) return
            if (maxval(abs(wb(j+1:n,j))) > tolerance) return
            wa(j+1:n,j) = cmplx(0.0_dp,0.0_dp,dp)
            wb(j+1:n,j) = cmplx(0.0_dp,0.0_dp,dp)
        end do
        wq = q
        wz = z
        original_a = matmul(wq,matmul(wa,transpose(conjg(wz))))
        original_b = matmul(wq,matmul(wb,transpose(conjg(wz))))
        beta_floor = 64.0_dp*epsilon(1.0_dp)*maxval(abs(bb))
        do i=1,n
            selected(i) = chosen(aa(i,i),bb(i,i),choice,beta_floor,selector)
        end do
        count = 0
        do j=1,n
            if (.not. selected(j)) cycle
            do k=j-1,count+1,-1
                call adjacent_swap(wa,wb,wq,wz,k,local_status)
                if (local_status /= linalg_status_success) return
                swap_flag = selected(k)
                selected(k) = selected(k+1)
                selected(k+1) = swap_flag
            end do
            count = count+1
        end do
        error_a = maxval(abs(matmul(wq,matmul(wa,transpose(conjg(wz))))-original_a))
        error_b = maxval(abs(matmul(wq,matmul(wb,transpose(conjg(wz))))-original_b))
        if (.not. ieee_is_finite(error_a) .or. .not. ieee_is_finite(error_b)) return
        if (max(error_a,error_b) > 1.0e-10_dp) return
        if (.not. finite_matrix(wa) .or. .not. finite_matrix(wb)) return
        final_a = wa*sa
        final_b = wb*sb
        if (.not. finite_matrix(final_a) .or. .not. finite_matrix(final_b)) return
        aa = final_a
        bb = final_b
        q = wq
        z = wz
        if (present(sdim)) sdim = count
        status = linalg_status_success
    end subroutine ordqz_factors

end module scifort_linalg_ordqz
