! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Toeplitz constructors and quadratic-time bordering/Levinson-type solves.
! Original implementation of the block-border identities for Toeplitz matrices;
! see docs/linalg_toeplitz.md for mathematical derivation and limitations.
module scifort_linalg_toeplitz
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_singular, linalg_status_no_convergence
    use scifort_linalg_dense, only : lu_factor, lu_solve
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    implicit none
    private

    public :: toeplitz, solve_toeplitz

    interface toeplitz
        module procedure toeplitz_real
        module procedure toeplitz_complex
    end interface toeplitz

    interface solve_toeplitz
        module procedure solve_toeplitz_real_vector
        module procedure solve_toeplitz_real_matrix
        module procedure solve_toeplitz_complex_vector
        module procedure solve_toeplitz_complex_matrix
        module procedure solve_toeplitz_real_complex_vector
        module procedure solve_toeplitz_real_complex_matrix
        module procedure solve_toeplitz_complex_real_vector
        module procedure solve_toeplitz_complex_real_matrix
    end interface solve_toeplitz

contains

    subroutine toeplitz_real(c, t, status, r)
        real(dp), intent(in) :: c(:) !! first column, length n (may be zero)
        real(dp), intent(out) :: t(:, :) !! Toeplitz matrix, shape (n,m)
        integer, intent(out) :: status !! linalg_status_* outcome
        real(dp), intent(in), optional :: r(:) !! first row, length m; r(1) ignored
        integer :: i, j, m, n

        t = 0.0_dp
        status = linalg_status_invalid_shape
        n = size(c)
        m = n
        if (present(r)) m = size(r)
        if (size(t,1) /= n .or. size(t,2) /= m) return
        status = linalg_status_no_convergence
        if (.not. all(ieee_is_finite(c))) return
        if (present(r)) then
            if (.not. all(ieee_is_finite(r))) return
        end if
        do j = 1, m
            do i = 1, n
                if (i >= j) then
                    t(i,j) = c(i-j+1)
                else if (present(r)) then
                    t(i,j) = r(j-i+1)
                else
                    t(i,j) = c(j-i+1)
                end if
            end do
        end do
        status = linalg_status_success
    end subroutine toeplitz_real

    subroutine toeplitz_complex(c, t, status, r)
        complex(dp), intent(in) :: c(:) !! first column, length n (may be zero)
        complex(dp), intent(out) :: t(:, :) !! Toeplitz matrix, shape (n,m)
        integer, intent(out) :: status !! linalg_status_* outcome
        complex(dp), intent(in), optional :: r(:) !! first row, length m; r(1) ignored
        integer :: i, j, m, n

        t = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_invalid_shape
        n = size(c)
        m = n
        if (present(r)) m = size(r)
        if (size(t,1) /= n .or. size(t,2) /= m) return
        status = linalg_status_no_convergence
        if (.not. finite_complex_vector(c)) return
        if (present(r)) then
            if (.not. finite_complex_vector(r)) return
        end if
        do j = 1, m
            do i = 1, n
                if (i >= j) then
                    t(i,j) = c(i-j+1)
                else if (present(r)) then
                    t(i,j) = r(j-i+1)
                else
                    t(i,j) = conjg(c(j-i+1))
                end if
            end do
        end do
        status = linalg_status_success
    end subroutine toeplitz_complex

    subroutine solve_toeplitz_real_vector(c, b, x, status, r, check_finite, used_dense)
        real(dp), intent(in) :: c(:) !! first column, length n
        real(dp), intent(in) :: b(:) !! right-hand side, length n
        real(dp), intent(out) :: x(:) !! solution vector, length n; zero on failure
        integer, intent(out) :: status !! linalg_status_* outcome
        real(dp), intent(in), optional :: r(:) !! first row, length n; r(1) ignored
        logical, intent(in), optional :: check_finite !! reject nonfinite inputs (default true)
        logical, intent(out), optional :: used_dense !! true if pivoted LU was used
        real(dp), allocatable :: rhs(:, :), sol(:, :)
        integer :: n

        x = 0.0_dp
        status = linalg_status_invalid_shape
        if (present(used_dense)) used_dense = .false.
        n = size(c)
        if (size(b) /= n .or. size(x) /= n) return
        allocate(rhs(n,1),sol(n,1))
        rhs(:,1) = b
        call solve_toeplitz_real_matrix(c,rhs,sol,status,r,check_finite,used_dense)
        if (status == linalg_status_success) x = sol(:,1)
    end subroutine solve_toeplitz_real_vector

    subroutine solve_toeplitz_complex_vector(c, b, x, status, r, check_finite, used_dense)
        complex(dp), intent(in) :: c(:) !! first column, length n
        complex(dp), intent(in) :: b(:) !! right-hand side, length n
        complex(dp), intent(out) :: x(:) !! solution vector, length n; zero on failure
        integer, intent(out) :: status !! linalg_status_* outcome
        complex(dp), intent(in), optional :: r(:) !! first row, length n; r(1) ignored
        logical, intent(in), optional :: check_finite !! reject nonfinite inputs (default true)
        logical, intent(out), optional :: used_dense !! true if pivoted LU was used
        complex(dp), allocatable :: rhs(:, :), sol(:, :)
        integer :: n

        x = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_invalid_shape
        if (present(used_dense)) used_dense = .false.
        n = size(c)
        if (size(b) /= n .or. size(x) /= n) return
        allocate(rhs(n,1),sol(n,1))
        rhs(:,1) = b
        call solve_toeplitz_complex_matrix(c,rhs,sol,status,r,check_finite,used_dense)
        if (status == linalg_status_success) x = sol(:,1)
    end subroutine solve_toeplitz_complex_vector

    subroutine solve_toeplitz_real_complex_vector(c, b, x, status, r, check_finite, used_dense)
        real(dp), intent(in) :: c(:) !! real first column, length n
        complex(dp), intent(in) :: b(:) !! complex right-hand side, length n
        complex(dp), intent(out) :: x(:) !! complex solution, length n; zero on failure
        integer, intent(out) :: status !! linalg_status_* outcome
        real(dp), intent(in), optional :: r(:) !! real first row, length n; r(1) ignored
        logical, intent(in), optional :: check_finite !! reject nonfinite inputs (default true)
        logical, intent(out), optional :: used_dense !! true if pivoted LU was used
        complex(dp), allocatable :: rhs(:, :), sol(:, :)
        integer :: n

        x = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_invalid_shape
        if (present(used_dense)) used_dense = .false.
        n = size(c)
        if (size(b) /= n .or. size(x) /= n) return
        allocate(rhs(n,1),sol(n,1))
        rhs(:,1) = b
        call solve_toeplitz_real_complex_matrix(c,rhs,sol,status,r,check_finite,used_dense)
        if (status == linalg_status_success) x = sol(:,1)
    end subroutine solve_toeplitz_real_complex_vector

    subroutine solve_toeplitz_complex_real_vector(c, b, x, status, r, check_finite, used_dense)
        complex(dp), intent(in) :: c(:) !! complex first column, length n
        real(dp), intent(in) :: b(:) !! real right-hand side, length n
        complex(dp), intent(out) :: x(:) !! complex solution, length n; zero on failure
        integer, intent(out) :: status !! linalg_status_* outcome
        complex(dp), intent(in), optional :: r(:) !! complex first row, length n; r(1) ignored
        logical, intent(in), optional :: check_finite !! reject nonfinite inputs (default true)
        logical, intent(out), optional :: used_dense !! true if pivoted LU was used
        complex(dp), allocatable :: rhs(:, :), sol(:, :)
        integer :: n

        x = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_invalid_shape
        if (present(used_dense)) used_dense = .false.
        n = size(c)
        if (size(b) /= n .or. size(x) /= n) return
        allocate(rhs(n,1),sol(n,1))
        rhs(:,1) = cmplx(b,kind=dp)
        call solve_toeplitz_complex_matrix(c,rhs,sol,status,r,check_finite,used_dense)
        if (status == linalg_status_success) x = sol(:,1)
    end subroutine solve_toeplitz_complex_real_vector

    subroutine solve_toeplitz_real_complex_matrix(c, b, x, status, r, check_finite, used_dense)
        real(dp), intent(in) :: c(:) !! real first column, length n
        complex(dp), intent(in) :: b(:, :) !! complex right-hand sides, shape (n,nrhs)
        complex(dp), intent(out) :: x(:, :) !! complex solutions, shape (n,nrhs)
        integer, intent(out) :: status !! linalg_status_* outcome
        real(dp), intent(in), optional :: r(:) !! real first row, length n; r(1) ignored
        logical, intent(in), optional :: check_finite !! reject nonfinite inputs (default true)
        logical, intent(out), optional :: used_dense !! true if pivoted LU was used

        if (present(r)) then
            call solve_toeplitz_complex_matrix(cmplx(c,kind=dp),b,x,status, &
                cmplx(r,kind=dp),check_finite,used_dense)
        else
            call solve_toeplitz_complex_matrix(cmplx(c,kind=dp),b,x,status, &
                check_finite=check_finite,used_dense=used_dense)
        end if
    end subroutine solve_toeplitz_real_complex_matrix

    subroutine solve_toeplitz_complex_real_matrix(c, b, x, status, r, check_finite, used_dense)
        complex(dp), intent(in) :: c(:) !! complex first column, length n
        real(dp), intent(in) :: b(:, :) !! real right-hand sides, shape (n,nrhs)
        complex(dp), intent(out) :: x(:, :) !! complex solutions, shape (n,nrhs)
        integer, intent(out) :: status !! linalg_status_* outcome
        complex(dp), intent(in), optional :: r(:) !! complex first row, length n; r(1) ignored
        logical, intent(in), optional :: check_finite !! reject nonfinite inputs (default true)
        logical, intent(out), optional :: used_dense !! true if pivoted LU was used

        call solve_toeplitz_complex_matrix(c,cmplx(b,kind=dp),x,status,r,check_finite,used_dense)
    end subroutine solve_toeplitz_complex_real_matrix

    subroutine solve_toeplitz_real_matrix(c, b, x, status, r, check_finite, used_dense)
        real(dp), intent(in) :: c(:) !! first column, length n
        real(dp), intent(in) :: b(:, :) !! right-hand sides, shape (n,nrhs)
        real(dp), intent(out) :: x(:, :) !! solutions, shape (n,nrhs); zero on failure
        integer, intent(out) :: status !! linalg_status_* outcome
        real(dp), intent(in), optional :: r(:) !! first row, length n; r(1) ignored
        logical, intent(in), optional :: check_finite !! reject nonfinite inputs (default true)
        logical, intent(out), optional :: used_dense !! true if pivoted LU was used
        real(dp), allocatable :: cc(:), rr(:), a(:), d(:), next_a(:), next_d(:)
        real(dp), allocatable :: work(:, :), bb(:, :)
        real(dp) :: scale, schur, delta_r, delta_c
        integer :: n, nrhs, k, j
        logical :: finite_check, fallback

        x = 0.0_dp
        status = linalg_status_invalid_shape
        if (present(used_dense)) used_dense = .false.
        n = size(c)
        nrhs = size(b,2)
        if (size(b,1) /= n .or. size(x,1) /= n .or. size(x,2) /= nrhs) return
        if (present(r)) then
            if (size(r) /= n) return
        end if
        finite_check = .true.
        if (present(check_finite)) finite_check = check_finite
        status = linalg_status_no_convergence
        if (finite_check) then
            if (.not. all(ieee_is_finite(c))) return
            if (.not. all(ieee_is_finite(b))) return
            if (present(r)) then
                if (.not. all(ieee_is_finite(r))) return
            end if
        end if
        if (n == 0 .or. nrhs == 0) then
            status = linalg_status_success
            return
        end if
        scale = maxval(abs(c))
        if (present(r)) then
            if (n > 1) scale = max(scale,maxval(abs(r(2:n))))
        end if
        status = linalg_status_singular
        if (scale == 0.0_dp) return
        status = linalg_status_no_convergence
        allocate(cc(n),rr(n),a(n),d(n),next_a(n),next_d(n),work(n,nrhs),bb(n,nrhs))
        cc = c/scale
        rr = cc
        if (present(r)) then
            if (n > 1) rr(2:n) = r(2:n)/scale
        end if
        bb = b/scale
        fallback = .false.
        if (abs(cc(1)) <= 32.0_dp*epsilon(1.0_dp)) fallback = .true.
        if (.not. fallback) then
            work(1,:) = bb(1,:)/cc(1)
            if (n > 1) then
                a(1) = rr(2)/cc(1)
                d(1) = cc(2)/cc(1)
            end if
            do k = 1, n-1
                schur = cc(1)-sum(cc(k+1:2:-1)*a(1:k))
                if (.not. ieee_is_finite(schur)) then
                    fallback = .true.
                    exit
                end if
                if (abs(schur) <= 32.0_dp*epsilon(1.0_dp)) then
                    fallback = .true.
                    exit
                end if
                do j = 1, nrhs
                    work(k+1,j) = (bb(k+1,j)- &
                        sum(cc(k+1:2:-1)*work(1:k,j)))/schur
                    work(1:k,j) = work(1:k,j)-a(1:k)*work(k+1,j)
                end do
                if (k == n-1) exit
                delta_r = rr(k+2)-sum(rr(2:k+1)*a(1:k))
                delta_c = cc(k+2)-sum(cc(2:k+1)*d(1:k))
                next_a(1) = delta_r/schur
                next_d(1) = delta_c/schur
                next_a(2:k+1) = a(1:k)-delta_r*d(k:1:-1)/schur
                next_d(2:k+1) = d(1:k)-delta_c*a(k:1:-1)/schur
                a(1:k+1) = next_a(1:k+1)
                d(1:k+1) = next_d(1:k+1)
            end do
        end if
        if (.not. fallback) then
            if (all(ieee_is_finite(work))) then
                x = work
                status = linalg_status_success
                return
            end if
            fallback = .true.
        end if
        if (present(used_dense)) used_dense = .true.
        call dense_toeplitz_real(cc,rr,bb,x,status)
    end subroutine solve_toeplitz_real_matrix

    subroutine solve_toeplitz_complex_matrix(c, b, x, status, r, check_finite, used_dense)
        complex(dp), intent(in) :: c(:) !! first column, length n
        complex(dp), intent(in) :: b(:, :) !! right-hand sides, shape (n,nrhs)
        complex(dp), intent(out) :: x(:, :) !! solutions, shape (n,nrhs); zero on failure
        integer, intent(out) :: status !! linalg_status_* outcome
        complex(dp), intent(in), optional :: r(:) !! first row, length n; r(1) ignored
        logical, intent(in), optional :: check_finite !! reject nonfinite inputs (default true)
        logical, intent(out), optional :: used_dense !! true if pivoted LU was used
        complex(dp), allocatable :: cc(:), rr(:), a(:), d(:), next_a(:), next_d(:)
        complex(dp), allocatable :: work(:, :), bb(:, :)
        complex(dp) :: schur, delta_r, delta_c
        real(dp) :: scale
        integer :: n, nrhs, k, j
        logical :: finite_check, fallback

        x = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_invalid_shape
        if (present(used_dense)) used_dense = .false.
        n = size(c)
        nrhs = size(b,2)
        if (size(b,1) /= n .or. size(x,1) /= n .or. size(x,2) /= nrhs) return
        if (present(r)) then
            if (size(r) /= n) return
        end if
        finite_check = .true.
        if (present(check_finite)) finite_check = check_finite
        status = linalg_status_no_convergence
        if (finite_check) then
            if (.not. finite_complex_vector(c)) return
            if (.not. finite_complex_matrix(b)) return
            if (present(r)) then
                if (.not. finite_complex_vector(r)) return
            end if
        end if
        if (n == 0 .or. nrhs == 0) then
            status = linalg_status_success
            return
        end if
        scale = maxval(abs(c))
        if (present(r)) then
            if (n > 1) scale = max(scale,maxval(abs(r(2:n))))
        end if
        status = linalg_status_singular
        if (scale == 0.0_dp) return
        status = linalg_status_no_convergence
        allocate(cc(n),rr(n),a(n),d(n),next_a(n),next_d(n),work(n,nrhs),bb(n,nrhs))
        cc = c/scale
        rr = conjg(cc)
        rr(1) = cc(1)
        if (present(r)) then
            if (n > 1) rr(2:n) = r(2:n)/scale
        end if
        bb = b/scale
        fallback = .false.
        if (abs(cc(1)) <= 32.0_dp*epsilon(1.0_dp)) fallback = .true.
        if (.not. fallback) then
            work(1,:) = bb(1,:)/cc(1)
            if (n > 1) then
                a(1) = rr(2)/cc(1)
                d(1) = cc(2)/cc(1)
            end if
            do k = 1, n-1
                ! These products are non-conjugating: this also handles
                ! general non-Hermitian and complex-symmetric Toeplitz matrices.
                schur = cc(1)-sum(cc(k+1:2:-1)*a(1:k))
                if (.not. finite_complex_scalar(schur)) then
                    fallback = .true.
                    exit
                end if
                if (abs(schur) <= 32.0_dp*epsilon(1.0_dp)) then
                    fallback = .true.
                    exit
                end if
                do j = 1, nrhs
                    work(k+1,j) = (bb(k+1,j)- &
                        sum(cc(k+1:2:-1)*work(1:k,j)))/schur
                    work(1:k,j) = work(1:k,j)-a(1:k)*work(k+1,j)
                end do
                if (k == n-1) exit
                delta_r = rr(k+2)-sum(rr(2:k+1)*a(1:k))
                delta_c = cc(k+2)-sum(cc(2:k+1)*d(1:k))
                next_a(1) = delta_r/schur
                next_d(1) = delta_c/schur
                next_a(2:k+1) = a(1:k)-delta_r*d(k:1:-1)/schur
                next_d(2:k+1) = d(1:k)-delta_c*a(k:1:-1)/schur
                a(1:k+1) = next_a(1:k+1)
                d(1:k+1) = next_d(1:k+1)
            end do
        end if
        if (.not. fallback) then
            if (finite_complex_matrix(work)) then
                x = work
                status = linalg_status_success
                return
            end if
            fallback = .true.
        end if
        if (present(used_dense)) used_dense = .true.
        call dense_toeplitz_complex(cc,rr,bb,x,status)
    end subroutine solve_toeplitz_complex_matrix

    subroutine dense_toeplitz_real(c, r, b, x, status)
        real(dp), intent(in) :: c(:) !! normalized Toeplitz first column, length n
        real(dp), intent(in) :: r(:) !! normalized Toeplitz first row, length n
        real(dp), intent(in) :: b(:, :) !! normalized right-hand sides, shape (n,nrhs)
        real(dp), intent(out) :: x(:, :) !! solutions or zeros, shape (n,nrhs)
        integer, intent(out) :: status !! linalg_status_* outcome
        real(dp), allocatable :: a(:, :), lu(:, :)
        integer, allocatable :: piv(:)
        integer :: n

        x = 0.0_dp
        n = size(c)
        allocate(a(n,n),lu(n,n),piv(n))
        call toeplitz_real(c,a,status,r)
        if (status /= linalg_status_success) return
        call lu_factor(a,lu,piv,status)
        if (status /= linalg_status_success) return
        call lu_solve(lu,piv,b,x,status)
    end subroutine dense_toeplitz_real

    subroutine dense_toeplitz_complex(c, r, b, x, status)
        complex(dp), intent(in) :: c(:) !! normalized Toeplitz first column, length n
        complex(dp), intent(in) :: r(:) !! normalized Toeplitz first row, length n
        complex(dp), intent(in) :: b(:, :) !! normalized right-hand sides, shape (n,nrhs)
        complex(dp), intent(out) :: x(:, :) !! solutions or zeros, shape (n,nrhs)
        integer, intent(out) :: status !! linalg_status_* outcome
        complex(dp), allocatable :: a(:, :), lu(:, :)
        integer, allocatable :: piv(:)
        integer :: n

        x = cmplx(0.0_dp,0.0_dp,dp)
        n = size(c)
        allocate(a(n,n),lu(n,n),piv(n))
        call toeplitz_complex(c,a,status,r)
        if (status /= linalg_status_success) return
        call lu_factor(a,lu,piv,status)
        if (status /= linalg_status_success) return
        call lu_solve(lu,piv,b,x,status)
    end subroutine dense_toeplitz_complex

    pure logical function finite_complex_scalar(z) result(ok)
        complex(dp), intent(in) :: z !! scalar tested for finite real and imaginary parts
        ok = ieee_is_finite(real(z,dp)) .and. ieee_is_finite(aimag(z))
    end function finite_complex_scalar

    pure logical function finite_complex_vector(z) result(ok)
        complex(dp), intent(in) :: z(:) !! vector tested for finite values
        ok = all(ieee_is_finite(real(z,dp))) .and. all(ieee_is_finite(aimag(z)))
    end function finite_complex_vector

    pure logical function finite_complex_matrix(z) result(ok)
        complex(dp), intent(in) :: z(:, :) !! matrix tested for finite values
        ok = all(ieee_is_finite(real(z,dp))) .and. all(ieee_is_finite(aimag(z)))
    end function finite_complex_matrix

end module scifort_linalg_toeplitz
