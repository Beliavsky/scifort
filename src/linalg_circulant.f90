! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Circulant constructors and FFT-diagonalized linear systems.
! Independent implementation of the discrete Fourier eigenbasis identity.
! Reuses SciFort's existing FFT plan; see docs/linalg_circulant.md.
module scifort_linalg_circulant
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_singular, linalg_status_no_convergence
    use scifort_signal_spectral, only : fft_plan, spectral_status_success
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    implicit none
    private

    public :: circulant, solve_circulant

    interface circulant
        module procedure circulant_real
        module procedure circulant_complex
    end interface circulant

    interface solve_circulant
        module procedure solve_circulant_real_vector
        module procedure solve_circulant_real_matrix
        module procedure solve_circulant_complex_vector
        module procedure solve_circulant_complex_matrix
        module procedure solve_circulant_real_complex_vector
        module procedure solve_circulant_real_complex_matrix
        module procedure solve_circulant_complex_real_vector
        module procedure solve_circulant_complex_real_matrix
    end interface solve_circulant

contains

    subroutine circulant_real(c, a, status)
        real(dp), intent(in) :: c(:) !! first column, length n (can be empty)
        real(dp), intent(out) :: a(:, :) !! circulant output, shape (n,n)
        integer, intent(out) :: status !! linalg_status_* result
        integer :: i, j, n

        a = 0.0_dp
        status = linalg_status_invalid_shape
        n = size(c)
        if (size(a,1) /= n .or. size(a,2) /= n) return
        status = linalg_status_no_convergence
        if (.not. all(ieee_is_finite(c))) return
        do j = 1, n
            do i = 1, n
                a(i,j) = c(modulo(i-j,n)+1)
            end do
        end do
        status = linalg_status_success
    end subroutine circulant_real

    subroutine circulant_complex(c, a, status)
        complex(dp), intent(in) :: c(:) !! first column, length n (can be empty)
        complex(dp), intent(out) :: a(:, :) !! circulant output, shape (n,n)
        integer, intent(out) :: status !! linalg_status_* result
        integer :: i, j, n

        a = (0.0_dp,0.0_dp)
        status = linalg_status_invalid_shape
        n = size(c)
        if (size(a,1) /= n .or. size(a,2) /= n) return
        status = linalg_status_no_convergence
        if (.not. finite_vector(c)) return
        do j = 1, n
            do i = 1, n
                a(i,j) = c(modulo(i-j,n)+1)
            end do
        end do
        status = linalg_status_success
    end subroutine circulant_complex

    subroutine solve_circulant_real_vector(c, b, x, status, singular, tol, check_finite)
        real(dp), intent(in) :: c(:) !! first column, length n
        real(dp), intent(in) :: b(:) !! right-hand side, length n
        real(dp), intent(out) :: x(:) !! real solution, zero on failure
        integer, intent(out) :: status !! linalg_status_* result
        character(*), intent(in), optional :: singular !! 'raise' (default) or 'lstsq'
        real(dp), intent(in), optional :: tol !! cutoff >= 0; default n*eps*max(abs(fft(c)))
        logical, intent(in), optional :: check_finite !! reject nonfinite operands, default true
        complex(dp), allocatable :: sol(:,:)
        integer :: n

        x = 0.0_dp
        status = linalg_status_invalid_shape
        n = size(c)
        if (size(b) /= n .or. size(x) /= n) return
        allocate(sol(n,1))
        call fft_solve(cmplx(c,kind=dp),reshape(cmplx(b,kind=dp),[n,1]), &
            sol,status,singular,tol,check_finite)
        if (status == linalg_status_success) x = real(sol(:,1),dp)
    end subroutine solve_circulant_real_vector

    subroutine solve_circulant_complex_vector(c, b, x, status, singular, tol, check_finite)
        complex(dp), intent(in) :: c(:) !! first column, length n
        complex(dp), intent(in) :: b(:) !! right-hand side, length n
        complex(dp), intent(out) :: x(:) !! complex solution, zero on failure
        integer, intent(out) :: status !! linalg_status_* result
        character(*), intent(in), optional :: singular !! 'raise' (default) or 'lstsq'
        real(dp), intent(in), optional :: tol !! cutoff >= 0; default n*eps*max(abs(fft(c)))
        logical, intent(in), optional :: check_finite !! reject nonfinite operands, default true
        complex(dp), allocatable :: sol(:,:)
        integer :: n

        x = (0.0_dp,0.0_dp)
        status = linalg_status_invalid_shape
        n = size(c)
        if (size(b) /= n .or. size(x) /= n) return
        allocate(sol(n,1))
        call fft_solve(c,reshape(b,[n,1]),sol,status,singular,tol,check_finite)
        if (status == linalg_status_success) x = sol(:,1)
    end subroutine solve_circulant_complex_vector

    subroutine solve_circulant_real_complex_vector(c, b, x, status, singular, tol, check_finite)
        real(dp), intent(in) :: c(:) !! real first column, length n
        complex(dp), intent(in) :: b(:) !! complex right-hand side, length n
        complex(dp), intent(out) :: x(:) !! complex solution, zero on failure
        integer, intent(out) :: status !! linalg_status_* result
        character(*), intent(in), optional :: singular !! 'raise' (default) or 'lstsq'
        real(dp), intent(in), optional :: tol !! cutoff >= 0; default n*eps*max(abs(fft(c)))
        logical, intent(in), optional :: check_finite !! reject nonfinite operands, default true
        complex(dp), allocatable :: sol(:,:)
        integer :: n

        x = (0.0_dp,0.0_dp)
        status = linalg_status_invalid_shape
        n = size(c)
        if (size(b) /= n .or. size(x) /= n) return
        allocate(sol(n,1))
        call fft_solve(cmplx(c,kind=dp),reshape(b,[n,1]),sol,status,singular,tol,check_finite)
        if (status == linalg_status_success) x = sol(:,1)
    end subroutine solve_circulant_real_complex_vector

    subroutine solve_circulant_complex_real_vector(c, b, x, status, singular, tol, check_finite)
        complex(dp), intent(in) :: c(:) !! complex first column, length n
        real(dp), intent(in) :: b(:) !! real right-hand side, length n
        complex(dp), intent(out) :: x(:) !! complex solution, zero on failure
        integer, intent(out) :: status !! linalg_status_* result
        character(*), intent(in), optional :: singular !! 'raise' (default) or 'lstsq'
        real(dp), intent(in), optional :: tol !! cutoff >= 0; default n*eps*max(abs(fft(c)))
        logical, intent(in), optional :: check_finite !! reject nonfinite operands, default true
        complex(dp), allocatable :: sol(:,:)
        integer :: n

        x = (0.0_dp,0.0_dp)
        status = linalg_status_invalid_shape
        n = size(c)
        if (size(b) /= n .or. size(x) /= n) return
        allocate(sol(n,1))
        call fft_solve(c,reshape(cmplx(b,kind=dp),[n,1]), &
            sol,status,singular,tol,check_finite)
        if (status == linalg_status_success) x = sol(:,1)
    end subroutine solve_circulant_complex_real_vector

    subroutine solve_circulant_real_matrix(c, b, x, status, singular, tol, &
            check_finite, baxis, outaxis)
        real(dp), intent(in) :: c(:) !! first column, length n
        real(dp), intent(in) :: b(:,:) !! RHS; length n on selected baxis
        real(dp), intent(out) :: x(:,:) !! solutions; length n on selected outaxis
        integer, intent(out) :: status !! linalg_status_* result
        character(*), intent(in), optional :: singular !! 'raise' (default) or 'lstsq'
        real(dp), intent(in), optional :: tol !! cutoff >= 0; default n*eps*max(abs(fft(c)))
        logical, intent(in), optional :: check_finite !! reject nonfinite operands, default true
        integer, intent(in), optional :: baxis !! RHS sample axis (1/-2 or 2/-1); default 1
        integer, intent(in), optional :: outaxis !! output sample axis (1/-2 or 2/-1); default 1
        complex(dp), allocatable :: sol(:,:)

        allocate(sol(size(x,1),size(x,2)))
        call axis_solve(cmplx(c,kind=dp),cmplx(b,kind=dp),sol,status, &
            singular,tol,check_finite,baxis,outaxis)
        x = 0.0_dp
        if (status == linalg_status_success) x = real(sol,dp)
    end subroutine solve_circulant_real_matrix

    subroutine solve_circulant_complex_matrix(c, b, x, status, singular, tol, &
            check_finite, baxis, outaxis)
        complex(dp), intent(in) :: c(:) !! first column, length n
        complex(dp), intent(in) :: b(:,:) !! RHS; length n on selected baxis
        complex(dp), intent(out) :: x(:,:) !! solutions; length n on selected outaxis
        integer, intent(out) :: status !! linalg_status_* result
        character(*), intent(in), optional :: singular !! 'raise' (default) or 'lstsq'
        real(dp), intent(in), optional :: tol !! cutoff >= 0; default n*eps*max(abs(fft(c)))
        logical, intent(in), optional :: check_finite !! reject nonfinite operands, default true
        integer, intent(in), optional :: baxis !! RHS sample axis (1/-2 or 2/-1); default 1
        integer, intent(in), optional :: outaxis !! output sample axis (1/-2 or 2/-1); default 1

        call axis_solve(c,b,x,status,singular,tol,check_finite,baxis,outaxis)
    end subroutine solve_circulant_complex_matrix

    subroutine solve_circulant_real_complex_matrix(c, b, x, status, singular, tol, &
            check_finite, baxis, outaxis)
        real(dp), intent(in) :: c(:) !! real first column, length n
        complex(dp), intent(in) :: b(:,:) !! complex RHS; length n on selected baxis
        complex(dp), intent(out) :: x(:,:) !! complex solution; length n on outaxis
        integer, intent(out) :: status !! linalg_status_* result
        character(*), intent(in), optional :: singular !! 'raise' (default) or 'lstsq'
        real(dp), intent(in), optional :: tol !! cutoff >= 0; default n*eps*max(abs(fft(c)))
        logical, intent(in), optional :: check_finite !! reject nonfinite operands, default true
        integer, intent(in), optional :: baxis !! RHS sample axis (1/-2 or 2/-1); default 1
        integer, intent(in), optional :: outaxis !! output sample axis (1/-2 or 2/-1); default 1

        call axis_solve(cmplx(c,kind=dp),b,x,status,singular,tol,check_finite,baxis,outaxis)
    end subroutine solve_circulant_real_complex_matrix

    subroutine solve_circulant_complex_real_matrix(c, b, x, status, singular, tol, &
            check_finite, baxis, outaxis)
        complex(dp), intent(in) :: c(:) !! complex first column, length n
        real(dp), intent(in) :: b(:,:) !! real RHS; length n on selected baxis
        complex(dp), intent(out) :: x(:,:) !! complex solution; length n on outaxis
        integer, intent(out) :: status !! linalg_status_* result
        character(*), intent(in), optional :: singular !! 'raise' (default) or 'lstsq'
        real(dp), intent(in), optional :: tol !! cutoff >= 0; default n*eps*max(abs(fft(c)))
        logical, intent(in), optional :: check_finite !! reject nonfinite operands, default true
        integer, intent(in), optional :: baxis !! RHS sample axis (1/-2 or 2/-1); default 1
        integer, intent(in), optional :: outaxis !! output sample axis (1/-2 or 2/-1); default 1

        call axis_solve(c,cmplx(b,kind=dp),x,status,singular,tol,check_finite,baxis,outaxis)
    end subroutine solve_circulant_complex_real_matrix

    ! 1-based axis convention. Negative axes are accepted as in NumPy:
    ! axis=-1 selects the last axis, and -2 the penultimate axis.
    pure integer function normalized_axis(axis) result(result_axis)
        integer, intent(in) :: axis !! specified axis; 1,2,-1,-2 valid
        result_axis = axis
        if (result_axis < 0) result_axis = 3+result_axis
        if (result_axis < 1 .or. result_axis > 2) result_axis = 0
    end function normalized_axis

    subroutine axis_solve(c, b, x, status, singular, tol, check_finite, baxis, outaxis)
        complex(dp), intent(in) :: c(:) !! first column, length n
        complex(dp), intent(in) :: b(:,:) !! RHS with vector axis selected by baxis
        complex(dp), intent(out) :: x(:,:) !! output with vector axis selected by outaxis
        integer, intent(out) :: status !! linalg_status_* result
        character(*), intent(in), optional :: singular !! near singular handling
        real(dp), intent(in), optional :: tol !! FFT eigenvalue cutoff
        logical, intent(in), optional :: check_finite !! reject nonfinite inputs, default true
        integer, intent(in), optional :: baxis !! RHS sample axis, default 1
        integer, intent(in), optional :: outaxis !! output sample axis, default 1
        integer :: ba, oa, n, nrhs
        complex(dp), allocatable :: rhs(:,:), sol(:,:)

        x = (0.0_dp,0.0_dp)
        status = linalg_status_invalid_shape
        ba = 1
        if (present(baxis)) ba = normalized_axis(baxis)
        oa = 1
        if (present(outaxis)) oa = normalized_axis(outaxis)
        if (ba == 0 .or. oa == 0) return
        n = size(c)
        if (size(b,ba) /= n) return
        nrhs = size(b,3-ba)
        if (size(x,oa) /= n .or. size(x,3-oa) /= nrhs) return
        allocate(rhs(n,nrhs),sol(n,nrhs))
        if (ba == 1) then
            rhs = b
        else
            rhs = transpose(b)
        end if
        call fft_solve(c,rhs,sol,status,singular,tol,check_finite)
        if (status /= linalg_status_success) return
        if (oa == 1) then
            x = sol
        else
            x = transpose(sol)
        end if
    end subroutine axis_solve

    subroutine fft_solve(c, b, x, status, singular, tol, check_finite)
        complex(dp), intent(in) :: c(:) !! first column, length n
        complex(dp), intent(in) :: b(:,:) !! n-by-nrhs right-hand sides
        complex(dp), intent(out) :: x(:,:) !! n-by-nrhs solutions; zero on failure
        integer, intent(out) :: status !! linalg_status_* result
        character(*), intent(in), optional :: singular !! 'raise' (default) or 'lstsq'
        real(dp), intent(in), optional :: tol !! explicit absolute eigenvalue cutoff, >= 0
        logical, intent(in), optional :: check_finite !! reject nonfinite inputs, default true
        type(fft_plan) :: plan
        complex(dp), allocatable :: eig(:), spectrum(:), solution(:,:)
        logical, allocatable :: small(:)
        real(dp) :: cutoff
        integer :: n, nrhs, j, i, fft_status
        logical :: require_finite, least_squares
        character(:), allocatable :: policy

        x = (0.0_dp,0.0_dp)
        status = linalg_status_invalid_shape
        n = size(c)
        nrhs = size(b,2)
        if (size(b,1) /= n .or. size(x,1) /= n .or. size(x,2) /= nrhs) return
        policy = 'raise'
        if (present(singular)) policy = trim(singular)
        if (policy /= 'raise' .and. policy /= 'lstsq') return
        least_squares = policy == 'lstsq'
        if (present(tol)) then
            if (.not. ieee_is_finite(tol)) return
            if (tol < 0.0_dp) return
        end if
        require_finite = .true.
        if (present(check_finite)) require_finite = check_finite
        status = linalg_status_no_convergence
        if (require_finite) then
            if (.not. finite_vector(c) .or. .not. finite_matrix(b)) return
        end if
        if (n == 0 .or. nrhs == 0) then
            status = linalg_status_success
            return
        end if
        call plan%initialize(n,fft_status)
        if (fft_status /= spectral_status_success) return
        allocate(eig(n),small(n),spectrum(n),solution(n,nrhs))
        eig = c
        call plan%execute(eig,fft_status)
        if (fft_status /= spectral_status_success) return
        if (.not. finite_vector(eig)) return
        cutoff = maxval(abs(eig))*(real(n,dp)*epsilon(1.0_dp))
        if (present(tol)) cutoff = tol
        small = abs(eig) <= cutoff
        if (any(small) .and. .not. least_squares) then
            status = linalg_status_singular
            return
        end if
        do j = 1, nrhs
            spectrum = b(:,j)
            call plan%execute(spectrum,fft_status)
            if (fft_status /= spectral_status_success) return
            do i = 1, n
                if (small(i)) then
                    spectrum(i) = (0.0_dp,0.0_dp)
                else
                    spectrum(i) = spectrum(i)/eig(i)
                end if
            end do
            call plan%inverse(spectrum,fft_status)
            if (fft_status /= spectral_status_success) return
            if (.not. finite_vector(spectrum)) return
            solution(:,j) = spectrum
        end do
        x = solution
        status = linalg_status_success
    end subroutine fft_solve

    pure logical function finite_vector(z) result(ok)
        complex(dp), intent(in) :: z(:) !! vector to check for nonfinite entries
        ok = all(ieee_is_finite(real(z,dp))) .and. all(ieee_is_finite(aimag(z)))
    end function finite_vector

    pure logical function finite_matrix(z) result(ok)
        complex(dp), intent(in) :: z(:,:) !! matrix to check for nonfinite entries
        ok = all(ieee_is_finite(real(z,dp))) .and. all(ieee_is_finite(aimag(z)))
    end function finite_matrix

end module scifort_linalg_circulant
