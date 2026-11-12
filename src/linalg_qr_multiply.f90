! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Apply the economic Q of a Householder QR factorization without forming Q.
! qr_apply reuses compact factors; qr_multiply computes them for convenience.
module scifort_linalg_qr_multiply
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    use scifort_linalg_qr, only : qr_factor
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    implicit none
    private
    public :: qr_multiply, qr_apply

    interface qr_multiply
        module procedure multiply_real_matrix
        module procedure multiply_complex_matrix
        module procedure multiply_real_vector
        module procedure multiply_complex_vector
        module procedure multiply_real_complex_matrix
        module procedure multiply_complex_real_matrix
        module procedure multiply_real_complex_vector
        module procedure multiply_complex_real_vector
    end interface
    interface qr_apply
        module procedure apply_real_matrix
        module procedure apply_complex_matrix
        module procedure apply_real_vector
        module procedure apply_complex_vector
    end interface

contains
    pure logical function finite_real(a) result(ok)
        real(dp), intent(in) :: a(:, :) !! real matrix with any dimensions
        ok = all(ieee_is_finite(a))
    end function finite_real

    pure logical function finite_complex(a) result(ok)
        complex(dp), intent(in) :: a(:, :) !! complex matrix with any dimensions
        ok = all(ieee_is_finite(real(a,dp))) .and. all(ieee_is_finite(aimag(a)))
    end function finite_complex

    subroutine apply_real_matrix(af, tau, c, product, status, mode, conjugate)
        real(dp), intent(in) :: af(:, :) !! m-by-n compact Householder QR matrix
        real(dp), intent(in) :: tau(:) !! min(m,n) Householder multipliers
        real(dp), intent(in) :: c(:, :) !! k-by-p (left) or p-by-m (right) matrix
        real(dp), intent(out) :: product(:, :) !! m-by-p (left) or p-by-k (right) result
        integer, intent(out) :: status !! linalg_status_* outcome
        character(*), intent(in), optional :: mode !! 'right' (default) or 'left'
        logical, intent(in), optional :: conjugate !! no effect on real Q; default false
        real(dp), allocatable :: work(:, :)
        real(dp) :: step
        integer :: m, n, k, i, j
        character(:), allocatable :: side

        product = 0.0_dp
        status = linalg_status_invalid_shape
        m = size(af,1)
        n = size(af,2)
        k = min(m,n)
        if (m < 1 .or. n < 1 .or. size(tau) /= k) return
        side = 'right'
        if (present(mode)) side = trim(adjustl(mode))
        if (side == 'left') then
            if (size(c,1) /= k) return
            if (size(product,1) /= m .or. size(product,2) /= size(c,2)) return
        else if (side == 'right') then
            if (size(c,2) /= m) return
            if (size(product,1) /= size(c,1) .or. size(product,2) /= k) return
        else
            return
        end if
        status = linalg_status_no_convergence
        if (.not. finite_real(af) .or. .not. all(ieee_is_finite(tau))) return
        if (.not. finite_real(c)) return
        if (side == 'left') then
            allocate(work(m,size(c,2)))
            work = 0.0_dp
            work(1:k,:) = c
            ! H(1)...H(k) times [C;0]: apply H(k),...,H(1).
            do i = k, 1, -1
                if (tau(i) == 0.0_dp) cycle
                do j = 1, size(work,2)
                    step = tau(i)*(work(i,j) + dot_product(af(i+1:m,i),work(i+1:m,j)))
                    work(i,j) = work(i,j) - step
                    work(i+1:m,j) = work(i+1:m,j) - af(i+1:m,i)*step
                end do
            end do
            if (.not. finite_real(work)) return
            product = work
        else
            work = c
            ! C H(1)...H(k), then retain the first k columns.
            do i = 1, k
                if (tau(i) == 0.0_dp) cycle
                do j = 1, size(work,1)
                    step = tau(i)*(work(j,i) + dot_product(work(j,i+1:m),af(i+1:m,i)))
                    work(j,i) = work(j,i) - step
                    work(j,i+1:m) = work(j,i+1:m) - step*af(i+1:m,i)
                end do
            end do
            if (.not. finite_real(work)) return
            product = work(:,1:k)
        end if
        status = linalg_status_success
    end subroutine apply_real_matrix

    subroutine apply_complex_matrix(af, tau, c, product, status, mode, conjugate)
        complex(dp), intent(in) :: af(:, :) !! m-by-n compact Householder QR matrix
        complex(dp), intent(in) :: tau(:) !! min(m,n) Householder multipliers
        complex(dp), intent(in) :: c(:, :) !! k-by-p (left) or p-by-m (right) matrix
        complex(dp), intent(out) :: product(:, :) !! m-by-p (left) or p-by-k (right) result
        integer, intent(out) :: status !! linalg_status_* outcome
        character(*), intent(in), optional :: mode !! 'right' (default) or 'left'
        logical, intent(in), optional :: conjugate !! use elementwise conjugate Q; default false
        complex(dp), allocatable :: work(:, :)
        complex(dp) :: step
        logical :: conjugated
        integer :: m, n, k, i, j
        character(:), allocatable :: side

        product = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_invalid_shape
        m = size(af,1)
        n = size(af,2)
        k = min(m,n)
        if (m < 1 .or. n < 1 .or. size(tau) /= k) return
        side = 'right'
        if (present(mode)) side = trim(adjustl(mode))
        if (side == 'left') then
            if (size(c,1) /= k) return
            if (size(product,1) /= m .or. size(product,2) /= size(c,2)) return
        else if (side == 'right') then
            if (size(c,2) /= m) return
            if (size(product,1) /= size(c,1) .or. size(product,2) /= k) return
        else
            return
        end if
        status = linalg_status_no_convergence
        if (.not. finite_complex(af)) return
        if (.not. all(ieee_is_finite(real(tau,dp)))) return
        if (.not. all(ieee_is_finite(aimag(tau)))) return
        if (.not. finite_complex(c)) return
        conjugated = .false.
        if (present(conjugate)) conjugated = conjugate
        if (side == 'left') then
            allocate(work(m,size(c,2)))
            work = cmplx(0.0_dp,0.0_dp,dp)
            work(1:k,:) = c
            do i = k, 1, -1
                if (tau(i) == cmplx(0.0_dp,0.0_dp,dp)) cycle
                do j = 1, size(work,2)
                    if (conjugated) then
                        step = conjg(tau(i))*(work(i,j) + &
                            dot_product(conjg(af(i+1:m,i)),work(i+1:m,j)))
                        work(i+1:m,j) = work(i+1:m,j) - conjg(af(i+1:m,i))*step
                    else
                        step = tau(i)*(work(i,j) + dot_product(af(i+1:m,i),work(i+1:m,j)))
                        work(i+1:m,j) = work(i+1:m,j) - af(i+1:m,i)*step
                    end if
                    work(i,j) = work(i,j) - step
                end do
            end do
            if (.not. finite_complex(work)) return
            product = work
        else
            work = c
            do i = 1, k
                if (tau(i) == cmplx(0.0_dp,0.0_dp,dp)) cycle
                do j = 1, size(work,1)
                    if (conjugated) then
                        step = conjg(tau(i))*(work(j,i) + &
                            sum(work(j,i+1:m)*conjg(af(i+1:m,i))))
                        work(j,i+1:m) = work(j,i+1:m) - step*af(i+1:m,i)
                    else
                        step = tau(i)*(work(j,i) + sum(work(j,i+1:m)*af(i+1:m,i)))
                        work(j,i+1:m) = work(j,i+1:m) - step*conjg(af(i+1:m,i))
                    end if
                    work(j,i) = work(j,i) - step
                end do
            end do
            if (.not. finite_complex(work)) return
            product = work(:,1:k)
        end if
        status = linalg_status_success
    end subroutine apply_complex_matrix

    subroutine apply_real_vector(af, tau, c, product, status, mode, conjugate)
        real(dp), intent(in) :: af(:, :) !! compact real QR factors
        real(dp), intent(in) :: tau(:) !! Householder multipliers
        real(dp), intent(in) :: c(:) !! k (left) or m (right) elements
        real(dp), intent(out) :: product(:) !! m (left) or k (right) elements
        integer, intent(out) :: status !! linalg_status_* outcome
        character(*), intent(in), optional :: mode !! 'right' (default) or 'left'
        logical, intent(in), optional :: conjugate !! ignored for real Q
        real(dp), allocatable :: result(:, :)
        character(:), allocatable :: side

        product = 0.0_dp
        side = 'right'
        if (present(mode)) side = trim(adjustl(mode))
        if (side == 'left') then
            allocate(result(size(product),1))
            call apply_real_matrix(af,tau,reshape(c,[size(c),1]),result,status,side,conjugate)
            if (status == linalg_status_success) product = result(:,1)
        else
            allocate(result(1,size(product)))
            call apply_real_matrix(af,tau,reshape(c,[1,size(c)]),result,status,side,conjugate)
            if (status == linalg_status_success) product = result(1,:)
        end if
    end subroutine apply_real_vector

    subroutine apply_complex_vector(af, tau, c, product, status, mode, conjugate)
        complex(dp), intent(in) :: af(:, :) !! compact complex QR factors
        complex(dp), intent(in) :: tau(:) !! Householder multipliers
        complex(dp), intent(in) :: c(:) !! k (left) or m (right) elements
        complex(dp), intent(out) :: product(:) !! m (left) or k (right) elements
        integer, intent(out) :: status !! linalg_status_* outcome
        character(*), intent(in), optional :: mode !! 'right' (default) or 'left'
        logical, intent(in), optional :: conjugate !! use elementwise conjugate Q
        complex(dp), allocatable :: result(:, :)
        character(:), allocatable :: side

        product = cmplx(0.0_dp,0.0_dp,dp)
        side = 'right'
        if (present(mode)) side = trim(adjustl(mode))
        if (side == 'left') then
            allocate(result(size(product),1))
            call apply_complex_matrix(af,tau,reshape(c,[size(c),1]),result,status,side,conjugate)
            if (status == linalg_status_success) product = result(:,1)
        else
            allocate(result(1,size(product)))
            call apply_complex_matrix(af,tau,reshape(c,[1,size(c)]),result,status,side,conjugate)
            if (status == linalg_status_success) product = result(1,:)
        end if
    end subroutine apply_complex_vector

    subroutine multiply_real_matrix(a, c, product, r, status, mode, pivots, conjugate)
        real(dp), intent(in) :: a(:, :) !! m-by-n matrix to QR-factor
        real(dp), intent(in) :: c(:, :) !! k-by-p (left) or p-by-m (right)
        real(dp), intent(out) :: product(:, :) !! m-by-p (left) or p-by-k (right)
        real(dp), intent(out) :: r(:, :) !! k-by-n economic upper trapezoidal factor
        integer, intent(out) :: status !! linalg_status_* outcome
        character(*), intent(in), optional :: mode !! 'right' (default) or 'left'
        integer, intent(out), optional :: pivots(:) !! one-based permutation; activates pivoting
        logical, intent(in), optional :: conjugate !! ignored for real Q
        real(dp), allocatable :: af(:, :), tau(:)
        integer :: m, n, k, i

        product = 0.0_dp
        r = 0.0_dp
        if (present(pivots)) pivots = 0
        status = linalg_status_invalid_shape
        m = size(a,1)
        n = size(a,2)
        k = min(m,n)
        if (m < 1 .or. n < 1) return
        if (size(r,1) /= k .or. size(r,2) /= n) return
        if (present(pivots)) then
            if (size(pivots) /= n) return
        end if
        allocate(af(m,n),tau(k))
        call qr_factor(a,af,tau,status,pivots)
        if (status /= linalg_status_success) return
        call apply_real_matrix(af,tau,c,product,status,mode,conjugate)
        if (status /= linalg_status_success) then
            if (present(pivots)) pivots = 0
            return
        end if
        do i = 1, k
            r(i,i:n) = af(i,i:n)
        end do
    end subroutine multiply_real_matrix

    subroutine multiply_complex_matrix(a, c, product, r, status, mode, pivots, conjugate)
        complex(dp), intent(in) :: a(:, :) !! m-by-n matrix to QR-factor
        complex(dp), intent(in) :: c(:, :) !! k-by-p (left) or p-by-m (right)
        complex(dp), intent(out) :: product(:, :) !! m-by-p (left) or p-by-k (right)
        complex(dp), intent(out) :: r(:, :) !! k-by-n economic upper trapezoidal factor
        integer, intent(out) :: status !! linalg_status_* outcome
        character(*), intent(in), optional :: mode !! 'right' (default) or 'left'
        integer, intent(out), optional :: pivots(:) !! one-based permutation; activates pivoting
        logical, intent(in), optional :: conjugate !! use elementwise conjugate Q
        complex(dp), allocatable :: af(:, :), tau(:)
        integer :: m, n, k, i

        product = cmplx(0.0_dp,0.0_dp,dp)
        r = cmplx(0.0_dp,0.0_dp,dp)
        if (present(pivots)) pivots = 0
        status = linalg_status_invalid_shape
        m = size(a,1)
        n = size(a,2)
        k = min(m,n)
        if (m < 1 .or. n < 1) return
        if (size(r,1) /= k .or. size(r,2) /= n) return
        if (present(pivots)) then
            if (size(pivots) /= n) return
        end if
        allocate(af(m,n),tau(k))
        call qr_factor(a,af,tau,status,pivots)
        if (status /= linalg_status_success) return
        call apply_complex_matrix(af,tau,c,product,status,mode,conjugate)
        if (status /= linalg_status_success) then
            if (present(pivots)) pivots = 0
            return
        end if
        do i = 1, k
            r(i,i:n) = af(i,i:n)
        end do
    end subroutine multiply_complex_matrix

    subroutine multiply_real_vector(a, c, product, r, status, mode, pivots, conjugate)
        real(dp), intent(in) :: a(:, :) !! real m-by-n coefficient matrix
        real(dp), intent(in) :: c(:) !! k (left) or m (right) elements
        real(dp), intent(out) :: product(:) !! m (left) or k (right) elements
        real(dp), intent(out) :: r(:, :) !! k-by-n R matrix
        integer, intent(out) :: status !! linalg_status_* outcome
        character(*), intent(in), optional :: mode !! 'right' (default) or 'left'
        integer, intent(out), optional :: pivots(:) !! one-based permutation; activates pivoting
        logical, intent(in), optional :: conjugate !! ignored for real Q
        real(dp), allocatable :: result(:, :)
        character(:), allocatable :: side

        product = 0.0_dp
        side = 'right'
        if (present(mode)) side = trim(adjustl(mode))
        if (side == 'left') then
            allocate(result(size(product),1))
            call multiply_real_matrix(a,reshape(c,[size(c),1]),result,r,status,side,pivots,conjugate)
            if (status == linalg_status_success) product = result(:,1)
        else
            allocate(result(1,size(product)))
            call multiply_real_matrix(a,reshape(c,[1,size(c)]),result,r,status,side,pivots,conjugate)
            if (status == linalg_status_success) product = result(1,:)
        end if
    end subroutine multiply_real_vector

    subroutine multiply_complex_vector(a, c, product, r, status, mode, pivots, conjugate)
        complex(dp), intent(in) :: a(:, :) !! complex m-by-n coefficient matrix
        complex(dp), intent(in) :: c(:) !! k (left) or m (right) elements
        complex(dp), intent(out) :: product(:) !! m (left) or k (right) elements
        complex(dp), intent(out) :: r(:, :) !! k-by-n R matrix
        integer, intent(out) :: status !! linalg_status_* outcome
        character(*), intent(in), optional :: mode !! 'right' (default) or 'left'
        integer, intent(out), optional :: pivots(:) !! one-based permutation; activates pivoting
        logical, intent(in), optional :: conjugate !! use elementwise conjugate Q
        complex(dp), allocatable :: result(:, :)
        character(:), allocatable :: side

        product = cmplx(0.0_dp,0.0_dp,dp)
        side = 'right'
        if (present(mode)) side = trim(adjustl(mode))
        if (side == 'left') then
            allocate(result(size(product),1))
            call multiply_complex_matrix(a,reshape(c,[size(c),1]),result,r,status,side,pivots,conjugate)
            if (status == linalg_status_success) product = result(:,1)
        else
            allocate(result(1,size(product)))
            call multiply_complex_matrix(a,reshape(c,[1,size(c)]),result,r,status,side,pivots,conjugate)
            if (status == linalg_status_success) product = result(1,:)
        end if
    end subroutine multiply_complex_vector

    subroutine multiply_real_complex_matrix(a, c, product, r, status, mode, pivots, conjugate)
        real(dp), intent(in) :: a(:, :) !! real m-by-n coefficient matrix
        complex(dp), intent(in) :: c(:, :) !! complex k-by-p (left) or p-by-m (right)
        complex(dp), intent(out) :: product(:, :) !! complex result, shape determined by mode
        complex(dp), intent(out) :: r(:, :) !! complex k-by-n upper factor
        integer, intent(out) :: status !! linalg_status_* outcome
        character(*), intent(in), optional :: mode !! 'right' (default) or 'left'
        integer, intent(out), optional :: pivots(:) !! one-based permutation; activates pivoting
        logical, intent(in), optional :: conjugate !! use elementwise conjugate Q
        call multiply_complex_matrix(cmplx(a,0.0_dp,dp),c,product,r,status,mode,pivots,conjugate)
    end subroutine multiply_real_complex_matrix

    subroutine multiply_complex_real_matrix(a, c, product, r, status, mode, pivots, conjugate)
        complex(dp), intent(in) :: a(:, :) !! complex m-by-n coefficient matrix
        real(dp), intent(in) :: c(:, :) !! real k-by-p (left) or p-by-m (right)
        complex(dp), intent(out) :: product(:, :) !! complex result, shape determined by mode
        complex(dp), intent(out) :: r(:, :) !! complex k-by-n upper factor
        integer, intent(out) :: status !! linalg_status_* outcome
        character(*), intent(in), optional :: mode !! 'right' (default) or 'left'
        integer, intent(out), optional :: pivots(:) !! one-based permutation; activates pivoting
        logical, intent(in), optional :: conjugate !! use elementwise conjugate Q
        call multiply_complex_matrix(a,cmplx(c,0.0_dp,dp),product,r,status,mode,pivots,conjugate)
    end subroutine multiply_complex_real_matrix

    subroutine multiply_real_complex_vector(a, c, product, r, status, mode, pivots, conjugate)
        real(dp), intent(in) :: a(:, :) !! real m-by-n coefficient matrix
        complex(dp), intent(in) :: c(:) !! complex k (left) or m (right) elements
        complex(dp), intent(out) :: product(:) !! complex result, shape determined by mode
        complex(dp), intent(out) :: r(:, :) !! complex k-by-n upper factor
        integer, intent(out) :: status !! linalg_status_* outcome
        character(*), intent(in), optional :: mode !! 'right' (default) or 'left'
        integer, intent(out), optional :: pivots(:) !! one-based permutation; activates pivoting
        logical, intent(in), optional :: conjugate !! use elementwise conjugate Q
        call multiply_complex_vector(cmplx(a,0.0_dp,dp),c,product,r,status,mode,pivots,conjugate)
    end subroutine multiply_real_complex_vector

    subroutine multiply_complex_real_vector(a, c, product, r, status, mode, pivots, conjugate)
        complex(dp), intent(in) :: a(:, :) !! complex m-by-n coefficient matrix
        real(dp), intent(in) :: c(:) !! real k (left) or m (right) elements
        complex(dp), intent(out) :: product(:) !! complex result, shape determined by mode
        complex(dp), intent(out) :: r(:, :) !! complex k-by-n upper factor
        integer, intent(out) :: status !! linalg_status_* outcome
        character(*), intent(in), optional :: mode !! 'right' (default) or 'left'
        integer, intent(out), optional :: pivots(:) !! one-based permutation; activates pivoting
        logical, intent(in), optional :: conjugate !! use elementwise conjugate Q
        call multiply_complex_vector(a,cmplx(c,0.0_dp,dp),product,r,status,mode,pivots,conjugate)
    end subroutine multiply_complex_real_vector
end module scifort_linalg_qr_multiply
