! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Generalized nonsymmetric eigenvalues and left/right eigenvectors from
! a complex two-sided QZ decomposition of a regular matrix pencil.
! The homogeneous pair (alpha,beta) represents beta*A*v=alpha*B*v.
module scifort_linalg_generalized_eigen
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence, linalg_status_singular
    use scifort_linalg_qz, only : qz
    use scifort_linalg_matrix_functions, only : matrix_function_max_order
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_positive_inf
    implicit none
    private
    public :: eig, eigvals

    interface eig
        module procedure eig_rr
        module procedure eig_rc
        module procedure eig_cr
        module procedure eig_cc
    end interface eig

    interface eigvals
        module procedure eigvals_rr
        module procedure eigvals_rc
        module procedure eigvals_cr
        module procedure eigvals_cc
    end interface eigvals

contains

    subroutine eig_rr(a, b, w, status, vr, vl, alpha, beta, method)
        real(dp), intent(in) :: a(:, :) !! first square matrix A, order 0..40
        real(dp), intent(in) :: b(:, :) !! second square matrix B, order 0..40
        complex(dp), intent(out) :: w(:) !! generalized eigenvalues, infinite when beta=0
        integer, intent(out) :: status !! linalg_status_* outcome
        complex(dp), intent(out), optional :: vr(:, :) !! normalized right eigenvectors
        complex(dp), intent(out), optional :: vl(:, :) !! normalized left eigenvectors
        complex(dp), intent(out), optional :: alpha(:) !! homogeneous eigenvalue numerators
        complex(dp), intent(out), optional :: beta(:) !! homogeneous eigenvalue denominators
        character(len=*), intent(in), optional :: method !! QZ algorithm: auto, iterative, inverse

        call generalized_eig_core(cmplx(a,kind=dp),cmplx(b,kind=dp),w,status,vr,vl,alpha,beta,method)
    end subroutine eig_rr

    subroutine eigvals_rr(a, b, w, status, alpha, beta, method)
        real(dp), intent(in) :: a(:, :) !! first square matrix A, order 0..40
        real(dp), intent(in) :: b(:, :) !! second square matrix B, order 0..40
        complex(dp), intent(out) :: w(:) !! generalized eigenvalues
        integer, intent(out) :: status !! linalg_status_* outcome
        complex(dp), intent(out), optional :: alpha(:) !! homogeneous numerators
        complex(dp), intent(out), optional :: beta(:) !! homogeneous denominators
        character(len=*), intent(in), optional :: method !! QZ algorithm selection

        call generalized_eig_core(cmplx(a,kind=dp),cmplx(b,kind=dp),w,status,alpha=alpha,beta=beta,method=method)
    end subroutine eigvals_rr

    subroutine eig_rc(a, b, w, status, vr, vl, alpha, beta, method)
        real(dp), intent(in) :: a(:, :) !! first square matrix A, order 0..40
        complex(dp), intent(in) :: b(:, :) !! second square matrix B, order 0..40
        complex(dp), intent(out) :: w(:) !! generalized eigenvalues, infinite when beta=0
        integer, intent(out) :: status !! linalg_status_* outcome
        complex(dp), intent(out), optional :: vr(:, :) !! normalized right eigenvectors
        complex(dp), intent(out), optional :: vl(:, :) !! normalized left eigenvectors
        complex(dp), intent(out), optional :: alpha(:) !! homogeneous eigenvalue numerators
        complex(dp), intent(out), optional :: beta(:) !! homogeneous eigenvalue denominators
        character(len=*), intent(in), optional :: method !! QZ algorithm: auto, iterative, inverse

        call generalized_eig_core(cmplx(a,kind=dp),b,w,status,vr,vl,alpha,beta,method)
    end subroutine eig_rc

    subroutine eigvals_rc(a, b, w, status, alpha, beta, method)
        real(dp), intent(in) :: a(:, :) !! first square matrix A, order 0..40
        complex(dp), intent(in) :: b(:, :) !! second square matrix B, order 0..40
        complex(dp), intent(out) :: w(:) !! generalized eigenvalues
        integer, intent(out) :: status !! linalg_status_* outcome
        complex(dp), intent(out), optional :: alpha(:) !! homogeneous numerators
        complex(dp), intent(out), optional :: beta(:) !! homogeneous denominators
        character(len=*), intent(in), optional :: method !! QZ algorithm selection

        call generalized_eig_core(cmplx(a,kind=dp),b,w,status,alpha=alpha,beta=beta,method=method)
    end subroutine eigvals_rc

    subroutine eig_cr(a, b, w, status, vr, vl, alpha, beta, method)
        complex(dp), intent(in) :: a(:, :) !! first square matrix A, order 0..40
        real(dp), intent(in) :: b(:, :) !! second square matrix B, order 0..40
        complex(dp), intent(out) :: w(:) !! generalized eigenvalues, infinite when beta=0
        integer, intent(out) :: status !! linalg_status_* outcome
        complex(dp), intent(out), optional :: vr(:, :) !! normalized right eigenvectors
        complex(dp), intent(out), optional :: vl(:, :) !! normalized left eigenvectors
        complex(dp), intent(out), optional :: alpha(:) !! homogeneous eigenvalue numerators
        complex(dp), intent(out), optional :: beta(:) !! homogeneous eigenvalue denominators
        character(len=*), intent(in), optional :: method !! QZ algorithm: auto, iterative, inverse

        call generalized_eig_core(a,cmplx(b,kind=dp),w,status,vr,vl,alpha,beta,method)
    end subroutine eig_cr

    subroutine eigvals_cr(a, b, w, status, alpha, beta, method)
        complex(dp), intent(in) :: a(:, :) !! first square matrix A, order 0..40
        real(dp), intent(in) :: b(:, :) !! second square matrix B, order 0..40
        complex(dp), intent(out) :: w(:) !! generalized eigenvalues
        integer, intent(out) :: status !! linalg_status_* outcome
        complex(dp), intent(out), optional :: alpha(:) !! homogeneous numerators
        complex(dp), intent(out), optional :: beta(:) !! homogeneous denominators
        character(len=*), intent(in), optional :: method !! QZ algorithm selection

        call generalized_eig_core(a,cmplx(b,kind=dp),w,status,alpha=alpha,beta=beta,method=method)
    end subroutine eigvals_cr

    subroutine eig_cc(a, b, w, status, vr, vl, alpha, beta, method)
        complex(dp), intent(in) :: a(:, :) !! first square matrix A, order 0..40
        complex(dp), intent(in) :: b(:, :) !! second square matrix B, order 0..40
        complex(dp), intent(out) :: w(:) !! generalized eigenvalues, infinite when beta=0
        integer, intent(out) :: status !! linalg_status_* outcome
        complex(dp), intent(out), optional :: vr(:, :) !! normalized right eigenvectors
        complex(dp), intent(out), optional :: vl(:, :) !! normalized left eigenvectors
        complex(dp), intent(out), optional :: alpha(:) !! homogeneous eigenvalue numerators
        complex(dp), intent(out), optional :: beta(:) !! homogeneous eigenvalue denominators
        character(len=*), intent(in), optional :: method !! QZ algorithm: auto, iterative, inverse

        call generalized_eig_core(a,b,w,status,vr,vl,alpha,beta,method)
    end subroutine eig_cc

    subroutine eigvals_cc(a, b, w, status, alpha, beta, method)
        complex(dp), intent(in) :: a(:, :) !! first square matrix A, order 0..40
        complex(dp), intent(in) :: b(:, :) !! second square matrix B, order 0..40
        complex(dp), intent(out) :: w(:) !! generalized eigenvalues
        integer, intent(out) :: status !! linalg_status_* outcome
        complex(dp), intent(out), optional :: alpha(:) !! homogeneous numerators
        complex(dp), intent(out), optional :: beta(:) !! homogeneous denominators
        character(len=*), intent(in), optional :: method !! QZ algorithm selection

        call generalized_eig_core(a,b,w,status,alpha=alpha,beta=beta,method=method)
    end subroutine eigvals_cc

    subroutine generalized_eig_core(a,b,w,status,vr,vl,alpha,beta,method)
        complex(dp), intent(in) :: a(:, :) !! first matrix A, square
        complex(dp), intent(in) :: b(:, :) !! second matrix B, square
        complex(dp), intent(out) :: w(:) !! eigenvalues of A v = lambda B v
        integer, intent(out) :: status !! linalg_status_* outcome
        complex(dp), intent(out), optional :: vr(:, :) !! normalized right eigenvectors
        complex(dp), intent(out), optional :: vl(:, :) !! normalized left eigenvectors
        complex(dp), intent(out), optional :: alpha(:) !! diagonal of first QZ factor
        complex(dp), intent(out), optional :: beta(:) !! diagonal of second QZ factor
        character(len=*), intent(in), optional :: method !! auto, iterative, inverse

        complex(dp), allocatable :: aa(:, :), bb(:, :), q(:, :), z(:, :), h(:, :)
        complex(dp), allocatable :: y(:), v(:)
        real(dp) :: sa, sb, scale_ab, residual, residual_bound
        complex(dp) :: scaled_alpha, scaled_beta
        integer :: n, j, code
        logical :: upper_pair, lower_pair

        w = cmplx(0.0_dp,0.0_dp,dp)
        if (present(vr)) vr = cmplx(0.0_dp,0.0_dp,dp)
        if (present(vl)) vl = cmplx(0.0_dp,0.0_dp,dp)
        if (present(alpha)) alpha = cmplx(0.0_dp,0.0_dp,dp)
        if (present(beta)) beta = cmplx(0.0_dp,0.0_dp,dp)
        status = linalg_status_invalid_shape
        n = size(a,1)
        if (size(a,2) /= n .or. n > matrix_function_max_order) return
        if (size(b,1) /= n .or. size(b,2) /= n) return
        if (size(w) /= n) return
        if (present(vr)) then
            if (size(vr,1) /= n .or. size(vr,2) /= n) return
        end if
        if (present(vl)) then
            if (size(vl,1) /= n .or. size(vl,2) /= n) return
        end if
        if (present(alpha)) then
            if (size(alpha) /= n) return
        end if
        if (present(beta)) then
            if (size(beta) /= n) return
        end if
        status = linalg_status_no_convergence
        if (.not. finite_complex(a) .or. .not. finite_complex(b)) return
        if (n == 0) then
            status = linalg_status_success
            return
        end if

        allocate(aa(n,n),bb(n,n),q(n,n),z(n,n))
        ! For a jointly triangular pair the generalized spectrum is already
        ! exactly on the diagonal. Do not introduce Schur roundoff into a
        ! potentially very nonnormal, large triangular pencil.
        if (present(method)) then
            if (trim(adjustl(method)) /= 'auto' .and. &
                trim(adjustl(method)) /= 'iterative' .and. &
                trim(adjustl(method)) /= 'inverse') then
                status = linalg_status_invalid_shape
                return
            end if
        end if
        upper_pair = .true.
        lower_pair = .true.
        do j=1,n-1
            if (any(abs(a(j+1:n,j)) > 0.0_dp) .or. &
                any(abs(b(j+1:n,j)) > 0.0_dp)) upper_pair = .false.
            if (any(abs(a(j,j+1:n)) > 0.0_dp) .or. &
                any(abs(b(j,j+1:n)) > 0.0_dp)) lower_pair = .false.
        end do
        if (upper_pair .or. lower_pair) then
            q = cmplx(0.0_dp,0.0_dp,dp)
            z = cmplx(0.0_dp,0.0_dp,dp)
            if (upper_pair) then
                aa = a
                bb = b
                do j=1,n
                    q(j,j) = cmplx(1.0_dp,0.0_dp,dp)
                    z(j,j) = cmplx(1.0_dp,0.0_dp,dp)
                end do
            else
                aa = a(n:1:-1,n:1:-1)
                bb = b(n:1:-1,n:1:-1)
                do j=1,n
                    q(j,n+1-j) = cmplx(1.0_dp,0.0_dp,dp)
                    z(j,n+1-j) = cmplx(1.0_dp,0.0_dp,dp)
                end do
            end if
        else
            call qz(a,b,aa,bb,q,z,code,method=method)
            status = code
            if (code /= linalg_status_success) return
        end if

        do j=1,n
            if (abs(aa(j,j)) <= 0.0_dp .and. abs(bb(j,j)) <= 0.0_dp) then
                status = linalg_status_singular
                return
            end if
            call safe_eigenvalue(aa(j,j),bb(j,j),w(j))
            if (present(alpha)) alpha(j) = aa(j,j)
            if (present(beta)) beta(j) = bb(j,j)
        end do

        if (.not. present(vr) .and. .not. present(vl)) then
            status = linalg_status_success
            return
        end if

        ! Divide A and B independently by their physical scales. The
        ! homogeneous eigenvector equation then remains well-scaled even
        ! for pairs whose coefficient magnitudes differ by > 1e200.
        sa = maxval(abs(aa))
        sb = maxval(abs(bb))
        if (sa <= 0.0_dp) sa = 1.0_dp
        if (sb <= 0.0_dp) sb = 1.0_dp
        allocate(h(n,n),y(n),v(n))
        residual_bound = 5.0e-8_dp
        do j=1,n
            scaled_alpha = aa(j,j)/sa
            scaled_beta = bb(j,j)/sb
            scale_ab = max(abs(scaled_alpha),abs(scaled_beta))
            if (scale_ab <= 0.0_dp) then
                status = linalg_status_singular
                return
            end if
            scaled_alpha = scaled_alpha/scale_ab
            scaled_beta = scaled_beta/scale_ab
            h = scaled_beta*(aa/sa)-scaled_alpha*(bb/sb)
            h(j,j) = cmplx(0.0_dp,0.0_dp,dp)
            if (present(vr)) then
                call triangular_pencil_vector(h,j,.false.,y)
                residual = scaled_residual(h,y,.false.)
                if (.not. ieee_is_finite(residual) .or. residual > residual_bound) then
                    status = linalg_status_no_convergence
                    return
                end if
                v = matmul(z,y)
                call normalize_vector(v,code)
                if (code /= linalg_status_success) then
                    status = code
                    return
                end if
                vr(:,j) = v
            end if
            if (present(vl)) then
                call triangular_pencil_vector(h,j,.true.,y)
                residual = scaled_residual(h,y,.true.)
                if (.not. ieee_is_finite(residual) .or. residual > residual_bound) then
                    status = linalg_status_no_convergence
                    return
                end if
                v = matmul(q,y)
                call normalize_vector(v,code)
                if (code /= linalg_status_success) then
                    status = code
                    return
                end if
                vl(:,j) = v
            end if
        end do
        status = linalg_status_success
    end subroutine generalized_eig_core

    pure logical function finite_complex(a) result(ok)
        complex(dp), intent(in) :: a(:, :) !! matrix with complex entries
        ok = all(ieee_is_finite(real(a,dp))) .and. all(ieee_is_finite(aimag(a)))
    end function finite_complex

    subroutine safe_eigenvalue(numerator,denominator,value)
        complex(dp), intent(in) :: numerator !! homogeneous alpha
        complex(dp), intent(in) :: denominator !! homogeneous beta
        complex(dp), intent(out) :: value !! finite or infinite generalized eigenvalue
        real(dp) :: an, bn, logmag, mag, inf
        complex(dp) :: phase

        an = abs(numerator)
        bn = abs(denominator)
        if (an <= 0.0_dp) then
            value = cmplx(0.0_dp,0.0_dp,dp)
        else if (bn <= 0.0_dp) then
            inf = ieee_value(0.0_dp,ieee_positive_inf)
            value = cmplx(inf,0.0_dp,dp)
        else
            ! Avoid overflow from a direct complex division on extreme scales.
            logmag = log(an)-log(bn)
            if (logmag > log(huge(1.0_dp))) then
                inf = ieee_value(0.0_dp,ieee_positive_inf)
                value = cmplx(inf,0.0_dp,dp)
            else if (logmag < log(tiny(1.0_dp))-35.0_dp) then
                value = cmplx(0.0_dp,0.0_dp,dp)
            else
                mag = exp(logmag)
                phase = (numerator/an)*conjg(denominator/bn)
                value = mag*phase
            end if
        end if
    end subroutine safe_eigenvalue

    subroutine triangular_pencil_vector(h,j,left,y)
        complex(dp), intent(in) :: h(:, :) !! upper triangular pencil factor, H(j,j)=0
        integer, intent(in) :: j !! target diagonal index, 1..n
        logical, intent(in) :: left !! if true, solve H^H y=0; else H y=0
        complex(dp), intent(out) :: y(:) !! (possibly scaled) eigenvector
        complex(dp) :: delta, subtotal
        real(dp) :: sep, magnitude
        integer :: i, n

        n = size(h,1)
        sep = 32.0_dp*epsilon(1.0_dp)*max(1.0_dp,maxval(abs(h)))
        y = cmplx(0.0_dp,0.0_dp,dp)
        y(j) = cmplx(1.0_dp,0.0_dp,dp)
        if (left) then
            do i=j+1,n
                subtotal = sum(conjg(h(j:i-1,i))*y(j:i-1))
                delta = conjg(h(i,i))
                if (abs(delta) < sep) then
                    if (abs(delta) > 0.0_dp) then
                        delta = sep*delta/abs(delta)
                    else
                        delta = cmplx(sep,0.0_dp,dp)
                    end if
                end if
                y(i) = -subtotal/delta
                magnitude = maxval(abs(y))
                if (magnitude > 1.0e60_dp) y = y/magnitude
            end do
        else
            do i=j-1,1,-1
                subtotal = sum(h(i,i+1:j)*y(i+1:j))
                delta = h(i,i)
                if (abs(delta) < sep) then
                    if (abs(delta) > 0.0_dp) then
                        delta = sep*delta/abs(delta)
                    else
                        delta = cmplx(sep,0.0_dp,dp)
                    end if
                end if
                y(i) = -subtotal/delta
                magnitude = maxval(abs(y))
                if (magnitude > 1.0e60_dp) y = y/magnitude
            end do
        end if
        magnitude = maxval(abs(y))
        if (magnitude > 0.0_dp) y = y/magnitude
    end subroutine triangular_pencil_vector

    pure real(dp) function scaled_residual(h,y,left) result(residual)
        complex(dp), intent(in) :: h(:, :) !! homogeneous triangular pencil
        complex(dp), intent(in) :: y(:) !! Schur-coordinate eigenvector
        logical, intent(in) :: left !! use conjugate transpose when true
        real(dp) :: scale

        scale = max(1.0_dp,maxval(abs(h)))
        if (left) then
            residual = maxval(abs(matmul(transpose(conjg(h/scale)),y))) / &
                max(1.0e-300_dp,maxval(abs(y)))
        else
            residual = maxval(abs(matmul(h/scale,y))) / &
                max(1.0e-300_dp,maxval(abs(y)))
        end if
    end function scaled_residual

    subroutine normalize_vector(v,status)
        complex(dp), intent(inout) :: v(:) !! eigenvector, returned with Euclidean norm 1
        integer, intent(out) :: status !! success or no convergence
        real(dp) :: scale, norm

        status = linalg_status_no_convergence
        scale = maxval(abs(v))
        if (.not. ieee_is_finite(scale) .or. scale <= 0.0_dp) return
        v = v/scale
        norm = sqrt(sum(real(conjg(v)*v,dp)))
        if (.not. ieee_is_finite(norm) .or. norm <= 0.0_dp) return
        v = v/norm
        status = linalg_status_success
    end subroutine normalize_vector
end module scifort_linalg_generalized_eigen
